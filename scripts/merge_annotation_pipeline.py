#!/usr/bin/env python3
"""
Fusionne les annotations SNP/indel, CNV/SV et microsatellites d'un pipeline
somatique (ANNOVAR, snpEff/SnpSift, AnnotSV, VisCap, MSIsensor...) en un seul
tableau final, à partir de fichiers passés en argument.

Chaque type de fichier est fourni via une option répétable :
    --vcf SAMPLE PATTERN        VCF annoté (ANNOVAR, snpEff/SnpSift...)
    --cnv SAMPLE PATTERN        Tableau CNV/SV annoté (AnnotSV, VisCap, Manta...)
    --microsat SAMPLE PATTERN   Sortie MSIsensor (microsatellites)

PATTERN peut être un chemin exact ou un motif glob (ex: "*.vcf") : tous les
fichiers correspondants sont lus et empilés pour le sample donné.

Exemple (reprend les fichiers du script d'origine pour l'échantillon 846352-kkl) :

    python merge_annotation_pipeline.py \\
        --cnv 846352-kkl "path/Annot_SV_CNV_annotation_21/846352-kkl.cns.vcf.annotated.tsv" \\
        --cnv 846352-kkl "path/Annot_SV_CNV_viscap_annotation_13/846352-kkl.viscapCNV.annotated.tsv" \\
        --cnv 846352-kkl "path/AnnotSV_SV_annotation_19/846352-kkl.somaticSV.annotated.tsv" \\
        --vcf 87334-M6   "path/snp_sift_annotation_18/846352-kkl.ann_snpEff_snpSift.vcf" \\
        --vcf 87334-M6   "path/Annovar_Snp_and_indel_annotation_16/846352-kkl.hg38_multianno.vcf" \\
        --microsat 846352-kkl "/path/Msisensor_08/846352-kkl_somatic.prefix_somatic" \\
        --outdir path/final_dataset_merge/846352-kkl.tsv
"""

import argparse
import glob
import os
import sys

import pandas as pd
from cyvcf2 import VCF


MICROSAT_COLUMNS = [
    "CHROM", "POS_start", "left_flank", "repeat_times",
    "repeat_unit_bases", "right_flank", "difference", "p-value", "FDR",
]


def _rule_from_path(filepath: str):

    parts = filepath.split(os.sep)
    return parts[-2] if len(parts) >= 2 else None


def _add_metadata(df: pd.DataFrame, filepath: str, sample: str, variant_caller: str) -> pd.DataFrame:
    df.insert(4, "output_file", os.path.basename(filepath))
    df.insert(5, "sample", sample)
    df.insert(6, "rule", _rule_from_path(filepath))
    df.insert(7,"variant calling tool ", variant_caller)
    return df


def _expand(pattern: str, kind: str) -> list:
    filepaths = sorted(glob.glob(pattern))
    if not filepaths:
        raise ValueError(f"Aucun fichier {kind} trouvé pour le pattern : {pattern}")
    return filepaths


def read_vcf(pattern: str, sample: str, variant_caller: str) -> pd.DataFrame:

    frames = []
    for filepath in _expand(pattern, "VCF"):
        vcf = VCF(filepath, samples=[sample])
        if sample not in vcf.samples:
            raise ValueError(
                f"Sample '{sample}' introuvable dans {filepath}. Disponibles : {vcf.samples}"
            )

        file_rows = []
        sample_idx = vcf.samples.index(sample)
        for variant in vcf:
            alt_alleles = variant.ALT or []
            first_alt = alt_alleles[0] if alt_alleles else ""
            row = {
                "CHROM": variant.CHROM,
                "POS": variant.POS,
                "REF": variant.REF,
                "ALT": ",".join(alt_alleles) if alt_alleles else None,
                # comparaison des longueurs d'allèles (et non nb d'ALT vs longueur REF)
                "VARIANT_TYPE": variant.var_type,
                "QUAL": variant.QUAL,
                "FILTERS": variant.FILTERS,
            }
            for key, value in variant.INFO:
                row[key] = value

            for fmt_field in variant.FORMAT:
              try:
                  values = variant.format(fmt_field)
                  if values is not None:
                      row[fmt_field] = values[sample_idx]
              except Exception:
                  pass

            
            file_rows.append(row)

        if not file_rows:
            raise ValueError(f"Aucun variant trouvé dans : {filepath}")

        df = pd.DataFrame(file_rows)
        frames.append(_add_metadata(df, filepath, sample,variant_caller))

    return pd.concat(frames, ignore_index=True)


def read_cnv(pattern: str, sample: str, variant_caller: str) -> pd.DataFrame:

    frames = []
    for filepath in _expand(pattern, "CNV/SV"):
        df = pd.read_csv(filepath, sep="\t")
        frames.append(_add_metadata(df, filepath, sample, variant_caller))

    return pd.concat(frames, ignore_index=True)


def read_microsatellite(pattern: str, sample: str, variant_caller: str) -> pd.DataFrame:

    frames = []
    for filepath in _expand(pattern, "microsatellite"):
        df = pd.read_csv(filepath, sep="\t", header=None, names=MICROSAT_COLUMNS, index_col=False)
        frames.append(_add_metadata(df, filepath, sample, variant_caller))

    return pd.concat(frames, ignore_index=True)



def read_manta_vcf(pattern: str, sample: str, caller: str) -> pd.DataFrame:

    frames = []
    for filepath in _expand(pattern, "Manta VCF"):
        vcf = VCF(filepath)
 
        file_rows = []
        for variant in vcf:
            alt_alleles = variant.ALT or []
            row = {
                "CHROM": variant.CHROM,
                "POS": variant.POS,
                "REF": variant.REF,
                "ALT": alt_alleles,
                "QUAL": variant.QUAL,
                "FILTER": variant.FILTERS,
            }
            for key, value in variant.INFO:
                row[key] = value
            file_rows.append(row)
 
        if not file_rows:
            raise ValueError(f"Aucun variant trouvé dans : {filepath}")
 
        df = pd.DataFrame(file_rows)
        frames.append(_add_metadata(df, filepath, sample, caller))
 
    return pd.concat(frames, ignore_index=True)


def format_vcf(df: pd.DataFrame) -> pd.DataFrame:
    df = df.rename(columns={"POS": "POS_start"})
    df.insert(2, "POS_end", df["POS_start"])
    return df

def format_manta_vcf(df: pd.DataFrame) -> pd.DataFrame:
    df = df.rename(columns={"POS": "POS_start"})
 
    pos_end = df["END"] if "END" in df.columns else pd.Series(index=df.index, dtype="float64")
    pos_end = pos_end.fillna(df["POS_start"])
    if "END" in df.columns:
        df = df.drop(columns=["END"])
    df.insert(2, "POS_end", pos_end)
 
    variant_type = df["SVTYPE"] if "SVTYPE" in df.columns else pd.Series(index=df.index, dtype="object")
    variant_type = variant_type.replace({"DEL": "SV_DEL", "DUP": "SV_DUP"})
    if "SVTYPE" in df.columns:
        df = df.drop(columns=["SVTYPE"])
    df.insert(3, "VARIANT_TYPE", variant_type)
 
    return df


def format_cnv(df: pd.DataFrame) -> pd.DataFrame:
    df = df.rename(columns={
        "SV_chrom": "CHROM",
        "SV_start": "POS_start",
        "SV_end": "POS_end",
        "SV_type": "VARIANT_TYPE",
    })

    if "VARIANT_TYPE" in df.columns:

        df["VARIANT_TYPE"] = df["VARIANT_TYPE"].replace({"DEL": "CNV_DEL", "DUP": "CNV_DUP"})

    df = df.drop(columns=["ID"], errors="ignore")

    for col, pos in (("AnnotSV_ID", 7), ("REF", 3), ("ALT", 4)):
        if col in df.columns:
            values = df.pop(col)
            df.insert(min(pos, len(df.columns)), col, values)

    if "CHROM" in df.columns:
        df["CHROM"] = df["CHROM"].astype(str).apply(
            lambda x: x if x.startswith("chr") else "chr" + x
        )

    return df


def format_microsatellite(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df.insert(2, "POS_end", df["POS_start"])
    df.insert(3, "VARIANT_TYPE", "microsatellite_instability")
    return df


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Fusionne les annotations SNP, CNV/SV et microsatellites en un seul tableau.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "--vcf", nargs=3, action="append", default=[], metavar=("SAMPLE", "PATTERN","VARIANT CALLING TOOL"),
        help="VCF annoté (ANNOVAR, snpEff/SnpSift...). + outil utilisé Répétable.",
    )
    
    parser.add_argument(
        "--cnv", nargs=3, action="append", default=[], metavar=("SAMPLE", "PATTERN", "VARIANT CALLING TOOL"),
        help="Tableau CNV/SV annoté (AnnotSV, VisCap, Manta...). + outil utilisé Répétable.",
    )

    parser.add_argument(
        "--manta_vcf", nargs=3, action="append", default=[], metavar=("SAMPLE", "PATTERN", "CALLER"),
        help="VCF Manta brut (SVTYPE/END en INFO, allèles symboliques/BND) + outil de calling + outil d'annotation. Répétable.",
    )

    parser.add_argument(
        "--microsat", nargs=3, action="append", default=[], metavar=("SAMPLE", "PATTERN", "VARIANT CALLING TOOL"),
        help="Sortie MSIsensor (microsatellites). + outil utilisé Répétable.",
    )
    parser.add_argument("--outdir", required=True, help="Chemin du fichier TSV final à écrire.")
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)

    if not (args.vcf or args.cnv or args.microsat):
        sys.exit("Aucun fichier fourni : utilisez --vcf, --cnv  et/ou --microsat (voir --help).")

    datasets = []

    for sample, pattern, variant_caller in args.vcf:
        datasets.append(format_vcf(read_vcf(pattern, sample, variant_caller )))

    for sample, pattern, variant_caller  in args.cnv :
        datasets.append(format_cnv(read_cnv(pattern, sample, variant_caller )))

    for sample, pattern, variant_caller  in args.manta_vcf :
        datasets.append(format_manta_vcf(read_manta_vcf(pattern, sample, variant_caller )))
    
    for sample, pattern, variant_caller  in args.microsat:
        datasets.append(format_microsatellite(read_microsatellite(pattern, sample, variant_caller )))

    annotation_dataset = pd.concat(datasets, ignore_index=True)

    outdir_parent = os.path.dirname(args.outdir)
    if outdir_parent:
        os.makedirs(outdir_parent, exist_ok=True)
    annotation_dataset.to_csv(args.outdir, sep="\t", index=False)
    print(f"{len(annotation_dataset)} lignes écrites dans {args.outdir}")


if __name__ == "__main__":
    main()
