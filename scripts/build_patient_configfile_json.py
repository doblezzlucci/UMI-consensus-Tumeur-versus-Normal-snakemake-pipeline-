#!/usr/bin/env python3


import argparse
import calendar
import csv
import json
import re
import sys
from datetime import date
from pathlib import Path


DEFAULT_LIBRARY = "lib1"
DEFAULT_PLATFORM = "illumina"

# Extrait l'identifiant d'échantillon juste avant "_R1.fastq.gz" / "_R2.fastq.gz"
SAMPLE_ID_RE = re.compile(r"_(\w+?)_R[12]\.fastq\.gz$")

# Détecte les colonnes circulant_N_R1 / circulant_N_R2 / circulant_N_date
CIRCULANT_COL_RE = re.compile(r"^circulant_(\d+)_R1$")


def parse_date_obj(raw_date: str) -> date:
   
    raw_date = raw_date.strip()
    parts = raw_date.split("/")
    if len(parts) != 3:
        raise ValueError(f"Format de date non reconnu : {raw_date!r}")
    day, month, year = parts
    if len(year) == 2:
        # Hypothèse : toutes les années à 2 chiffres sont en 20xx
        year = "20" + year
    return date(int(year), int(month), int(day))


def add_months(d: date, months: int) -> date:
    
    total_month_index = d.month - 1 + months
    year = d.year + total_month_index // 12
    month = total_month_index % 12 + 1
    day = min(d.day, calendar.monthrange(year, month)[1])
    return date(year, month, day)


def format_yyyymmdd(d: date) -> str:
    return d.strftime("%Y%m%d")


def format_ddmmyyyy(d: date) -> str:
    return d.strftime("%d/%m/%Y")


def extract_sample_id(fastq_path: str) -> str:
    match = SAMPLE_ID_RE.search(fastq_path)
    if not match:
        raise ValueError(f"Impossible d'extraire l'ID d'échantillon depuis : {fastq_path!r}")
    return match.group(1)


def insert_date_in_path(fastq_path: str, date_str: str) -> str:
    
    return re.sub(r"_(R[12])\.fastq\.gz$", rf"_{date_str}_\1.fastq.gz", fastq_path)


def build_sample_entry(role: str, sample_id: str, r1: str, r2: str, umi_structure: str,
                        date_str: str, date_display: str) -> dict:
    fastq_r1 = insert_date_in_path(r1, date_str)
    fastq_r2 = insert_date_in_path(r2, date_str)
    return {
        role: sample_id,
        "fastq_r1": fastq_r1,
        "fastq_r2": fastq_r2,
        "read_structures": [umi_structure, umi_structure],
        "read_group_id": sample_id,
        "sample_name": sample_id,
        "library": DEFAULT_LIBRARY,
        "platform": DEFAULT_PLATFORM,
        "date": date_display,
    }


def find_circulant_columns(fieldnames):
    
    indices = set()
    for name in fieldnames:
        m = CIRCULANT_COL_RE.match(name)
        if m:
            indices.add(int(m.group(1)))
    return sorted(indices)


def build_patient_from_row(row: dict, circulant_indices) -> tuple[str, dict]:
    patient_id = row["patient_id"].strip()
    umi_structure = row["UMI_structure"].strip()
    main_date_obj = parse_date_obj(row["Date"])
    main_date = format_yyyymmdd(main_date_obj)
    main_date_display = format_ddmmyyyy(main_date_obj)

    normal_id = extract_sample_id(row["normal_R1"])
    tumeur_id = extract_sample_id(row["tumeur_solide_R1"])

    normal_entry = build_sample_entry(
        "normal", normal_id, row["normal_R1"], row["normal_R2"], umi_structure,
        main_date, main_date_display,
    )
    tumeur_entry = build_sample_entry(
        "tumor", tumeur_id, row["tumeur_solide_R1"], row["tumeur_solide_R2"], umi_structure,
        main_date, main_date_display,
    )

    cfdna_samples = {}
    for idx in circulant_indices:
        r1_col = f"circulant_{idx}_R1"
        r2_col = f"circulant_{idx}_R2"
        month_col = f"circulant_{idx}_date"  # nombre de mois de suivi depuis la date principale
        if not row.get(r1_col) or not row.get(r2_col):
            continue
        months_val = int(row[month_col].strip()) if row.get(month_col, "").strip() else 0
        circ_date_obj = add_months(main_date_obj, months_val)
        circ_date = format_yyyymmdd(circ_date_obj)
        circ_date_display = format_ddmmyyyy(circ_date_obj)
        circ_id = extract_sample_id(row[r1_col])
        cfdna_samples[circ_id] = build_sample_entry(
            "tumor", circ_id, row[r1_col], row[r2_col], umi_structure,
            circ_date, circ_date_display,
        )

    patient_data = {
        "normal": {
            "gDNA_samples": {
                normal_id: normal_entry,
            }
        },
        "tumor": {
            "gDNA_samples": {
                tumeur_id: tumeur_entry,
            }
        },
    }
    if cfdna_samples:
        patient_data["tumor"]["cfDNA_samples"] = cfdna_samples

    return patient_id, patient_data


def deep_merge(base: dict, addition: dict) -> dict:
    """Fusionne 'addition' dans 'base', récursivement, sans écraser les clés voisines."""
    for key, value in addition.items():
        if key in base and isinstance(base[key], dict) and isinstance(value, dict):
            deep_merge(base[key], value)
        else:
            base[key] = value
    return base


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samplesheet", required=True, help="Chemin du sample sheet TSV")
    parser.add_argument("--json", default=None, help="Fichier JSON existant à mettre à jour (optionnel)")
    parser.add_argument("--output", required=True, help="Chemin du fichier JSON de sortie")
    args = parser.parse_args()

    # Charger le JSON existant, ou démarrer avec une structure vide
    if args.json and Path(args.json).exists():
        with open(args.json, "r", encoding="utf-8") as f:
            content = f.read()
        try:
            data = json.loads(content)
        except json.JSONDecodeError as e:
            lines = content.splitlines()
            start = max(0, e.lineno - 3)
            end = min(len(lines), e.lineno + 2)
            context = "\n".join(
                f"{i + 1:>4}: {lines[i]}" for i in range(start, end)
            )
            print(
                f"Erreur : {args.json} n'est pas un JSON valide "
                f"(ligne {e.lineno}, colonne {e.colno}) : {e.msg}\n"
                f"Contexte :\n{context}\n\n"
                f"Vérifiez notamment une virgule manquante/en trop entre deux "
                f"entrées, ou une accolade/crochet non fermé autour de cette ligne.",
                file=sys.stderr,
            )
            sys.exit(1)
    else:
        data = {}
    data.setdefault("patients", {})

    with open(args.samplesheet, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        circulant_indices = find_circulant_columns(reader.fieldnames)
        for row in reader:
            if not row.get("patient_id"):
                continue
            patient_id, patient_data = build_patient_from_row(row, circulant_indices)
            data["patients"].setdefault(patient_id, {})
            deep_merge(data["patients"][patient_id], patient_data)

    with open(args.output, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)
        f.write("\n")

    print(f"OK : {len(data['patients'])} patient(s) dans {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
