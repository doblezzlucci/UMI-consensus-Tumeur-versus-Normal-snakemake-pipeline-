import os
import gzip


def build_samples_and_comparisons(config):
    samples = {}
    comparisons = []

    for patient_id, pdata in config["patients"].items():

        normal_samples = pdata["normal"]["gDNA_samples"]

        if len(normal_samples) != 1:
            raise ValueError(
                f"{patient_id}: attendu exactement 1 normal, trouvé {len(normal_samples)}"
            )

        normal_id, normal_cfg = next(iter(normal_samples.items()))

        samples[(patient_id, normal_id)] = {
            **normal_cfg,
            "role": "normal",
            "sample_type": "gDNA"
        }

        tumor_groups = pdata["tumor"]

        for sample_type, tumor_samples in tumor_groups.items():
            for tumor_id, tumor_cfg in tumor_samples.items():

                samples[(patient_id, tumor_id)] = {
                    **tumor_cfg,
                    "role": "tumor",
                    "sample_type": sample_type
                }

                comparisons.append({
                    "patient": patient_id,
                    "normal": normal_id,
                    "tumor": tumor_id,
                    "tumor_type": sample_type
                })

    return samples, comparisons

SAMPLES, COMPARISONS = build_samples_and_comparisons(config)


def vcf_has_variants(vcf_path: str) -> bool:

    opener = gzip.open if vcf_path.endswith(".gz") else open
    with opener(vcf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            if line.strip():
                return True
    return False





def valid_comparaison(comparaisons: list) -> list:
    cfdna_count = sum(1 for c in comparaisons if c.get("tumor_type") == "cfDNA_samples")
    return comparaisons if cfdna_count >= 2 else []

VALID_COMPARISONS = valid_comparaison(COMPARISONS)

def ref():
    return config["reference_fasta"]

def ref_bed():
    return config["refrence_bed"]

def ref_gnomad():
    return config["reference_annotation_gnomad"]

def ref_satellite_list():
    return config["microsatellite_list"]

def ref_accessible_bed():
    return config["accessible_bed"]

def ref_gene_labels():
    return config["reference_gene_labels"]

def ref_annotated_bed():
    return config["annotated_bed"]


def recherche_databse(a, b, config, c, d):
    result = config.get(a, {}).get(b)
    if d == "name" :
        return c.upper()+"_"

    elif d == "path":
        return result.get(c)
    else:
        raise ValueError(f"Invalid value for d: {d!r}")



def sample_cfig(patient, sample):
    return SAMPLES[(patient, sample)]

def sample_cfg(wc):
    return sample_cfig(wc.patient, wc.sample)

def comparison_cfg(patient, tumor):
    for comp in COMPARISONS:
        if comp["patient"] == patient and comp["tumor"] == tumor:
            return comp
    raise KeyError(f"No comparison found for {patient}, {tumor}")


def tumor_cfg(wc):
    return sample_cfig(wc.patient, wc.tumor)

def normal_cfg(wc):
    return sample_cfig(wc.patient, wc.normal)


def sample_name(wc):
    return sample_cfg(wc)["sample_name"]
    
def fastq_r1(wc):
    return sample_cfg(wc)["fastq_r1"]

def fastq_r2(wc):
    return sample_cfg(wc)["fastq_r2"]

def consensus_bam(patient, sample):
    return os.path.join(
        config["output_dir"],
        patient,
        sample,
        "19_consensus_final",
        f"{sample}.consensus.zipped.coord.bam"
    )

def consensus_bai(patient, sample):
    return consensus_bam(patient, sample) + ".bai"



PAIR_FINAL_OUTPUTS = [
    os.path.join(
        config["output_dir"],
        comp["patient"],
        f"{comp['normal']}_vs_{comp['tumor']}",
        "Annot_SV_CNV_annotation_21",
        "done.py")
    for comp in COMPARISONS
]

PAIR_FINAL_OUTPUTS_cfDNA = [
    os.path.join(
        config["output_dir"],
        comp["patient"],
        comp['tumor'],
        "Annot_SV_CNV_viscap_annotation_13",
        "done.py")
    for comp in VALID_COMPARISONS
]

CONSENSUS_BAMS = [
    consensus_bam(patient, sample)
    for patient, sample in SAMPLES.keys()
]

CNV_DEPTH_OUT=[
    os.path.join(
        config["output_dir"],
        v["patient"],
        v["tumor"],
        "CNV_depth_coverage_10",
        "done.py")
    for v in VALID_COMPARISONS
]

rule all:
    input:
        CONSENSUS_BAMS,
        PAIR_FINAL_OUTPUTS,
        PAIR_FINAL_OUTPUTS_cfDNA,
        CNV_DEPTH_OUT
        
        
        
        

# 01 - FASTQ -> uBAM brut
rule fastq_to_ubam_raw_01:
    input:
        fastq1=lambda wc: sample_cfg(wc)["fastq_r1"],
        fastq2=lambda wc: sample_cfg(wc)["fastq_r2"]
    output:
        bam=os.path.join(config["output_dir"],"{patient}","{sample}", "01_fastq_to_ubam","{sample}.unmapped.bam")
    threads: 1
    params:
        read_structure1=lambda wc: sample_cfg(wc).get("read_structures", ["5M+T", "5M+T"])[0],
        read_structure2=lambda wc: sample_cfg(wc).get("read_structures", ["5M+T", "5M+T"])[1],
        rgid=lambda wc: sample_cfg(wc).get("read_group_id", wc.sample),
        sample_name=lambda wc: sample_cfg(wc).get("sample_name", wc.sample),
        library=lambda wc: sample_cfg(wc).get("library", "lib1"),
        platform=lambda wc: sample_cfg(wc).get("platform", "illumina")
    shell:
        r"""
        mkdir -p $(dirname {output.bam})
        java -Xmx4g -jar /usr/share/java/fgbio-2.2.1.jar FastqToBam \
          --input {input.fastq1} \
          --input {input.fastq2} \
          --read-structures {params.read_structure1} {params.read_structure2} \
          --output {output.bam} \
          --read-group-id {params.rgid} \
          --sample {params.sample_name} \
          --library {params.library} \
          --platform {params.platform}
        """


# 02 - uBAM -> FASTQ brut
rule ubam_to_fastq_raw_02:
    input:
        bam=rules.fastq_to_ubam_raw_01.output.bam
    output:
        fastq1=os.path.join(config["output_dir"],"{patient}",  "{sample}", "02_ubam_fastq", "R1", "{sample}.R1.fastq.gz"),
        fastq2=os.path.join(config["output_dir"],"{patient}",  "{sample}", "02_ubam_fastq", "R2", "{sample}.R2.fastq.gz")
    threads: 4
    shell:
        r"""
        mkdir -p $(dirname {output.fastq1}) $(dirname {output.fastq2})
        java -XX:ParallelGCThreads={threads} -Xmx8g -Xms8g -jar /usr/share/java/picard-2.8.1.jar SamToFastq \
          I={input.bam} \
          F={output.fastq1} \
          F2={output.fastq2}
        """


# 03 - Alignement brut
rule alignment_raw_03:
    input:
        fastq1=rules.ubam_to_fastq_raw_02.output.fastq1,
        fastq2=rules.ubam_to_fastq_raw_02.output.fastq2
    output:
        sam=os.path.join(config["output_dir"],"{patient}",  "{sample}", "03_align_raw", "{sample}.aligned.sam")
    threads: 8
    params:
        rg=lambda wc: '@RG\\tID:{0}\\tSM:{1}\\tPL:{2}'.format(
            sample_cfg(wc).get("read_group_id", wc.sample),
            sample_cfg(wc).get("sample_name", wc.sample),
            sample_cfg(wc).get("platform", "ILLUMINA").upper(),
        ),
        reference=lambda wc: ref()
    shell:
        r"""
        mkdir -p $(dirname {output.sam})
        bwa-mem2 mem -t {threads} \
          -R "{params.rg}" \
          {params.reference} \
          {input.fastq1} \
          {input.fastq2} \
          > {output.sam}
        """


# 04 - SAM -> BAM brut
rule sam_to_bam_raw_04:
    input:
        sam=rules.alignment_raw_03.output.sam
    output:
        bam=os.path.join(config["output_dir"],"{patient}",  "{sample}", "04_sam_to_bam_raw", "{sample}.aligned.bam")
    threads: 8
    shell:
        r"""
        mkdir -p $(dirname {output.bam})
        samtools view -@ {threads} -bS {input.sam} > {output.bam}
        """


# 05 - Supprime secondary/supplementary
rule remove_secondary_supplementary_raw_05:
    input:
        bam=rules.sam_to_bam_raw_04.output.bam
    output:
        bam=os.path.join(config["output_dir"],"{patient}",  "{sample}", "05_primary_only", "{sample}.aligned.primary.bam")
    threads: 4
    shell:
        r"""
        mkdir -p $(dirname {output.bam})
        samtools view -@ {threads} -h -F 0x900 {input.bam} | samtools view -@ {threads} -b -o {output.bam}
        """


# 06 - Tri queryname
rule sort_unmapped_queryname_raw_06:
    input:
        bam=rules.fastq_to_ubam_raw_01.output.bam
    output:
        bam=os.path.join(config["output_dir"],"{patient}",  "{sample}", "06_queryname_sort", "{sample}.unmapped.queryname.bam")
    threads: 4
    shell:
        r"""
        mkdir -p $(dirname {output.bam})
        java -Xmx4g -jar /usr/share/java/picard.jar SortSam \
          I={input.bam} \
          O={output.bam} \
          SORT_ORDER=queryname
        """


rule sort_primary_queryname_raw_06:
    input:
        bam=rules.remove_secondary_supplementary_raw_05.output.bam
    output:
        bam=os.path.join(config["output_dir"],"{patient}",  "{sample}", "06_queryname_sort", "{sample}.aligned.primary.queryname.bam")
    threads: 4
    shell:
        r"""
        mkdir -p $(dirname {output.bam})
        java -Xmx4g -jar /usr/share/java/picard.jar SortSam \
          I={input.bam} \
          O={output.bam} \
          SORT_ORDER=queryname
        """


# 07 - Retire @RG de l'en-tête du BAM aligné avant zipper
rule remove_rg_header_raw_07:
    input:
        bam=rules.sort_primary_queryname_raw_06.output.bam
    output:
        header=os.path.join(config["output_dir"],"{patient}",  "{sample}", "07_reheader", "{sample}.aligned.noRG.header.sam")
    threads: 1
    shell:
        r"""
        mkdir -p $(dirname {output.header})
        samtools view -H {input.bam} | grep -v '^@RG' > {output.header}
        """


rule reheader_primary_queryname_raw_07:
    input:
        header=rules.remove_rg_header_raw_07.output.header,
        bam=rules.sort_primary_queryname_raw_06.output.bam
    output:
        bam=os.path.join(config["output_dir"],"{patient}",  "{sample}", "07_reheader", "{sample}.aligned.primary.queryname.noRG.bam")
    threads: 1
    shell:
        r"""
        mkdir -p $(dirname {output.bam})
        samtools reheader {input.header} {input.bam} > {output.bam}
        """


# 08 - Zipper raw
rule zipper_bams_raw_08:
    input:
        unmapped=rules.sort_unmapped_queryname_raw_06.output.bam,
        aligned=rules.reheader_primary_queryname_raw_07.output.bam
    output:
        bam=os.path.join(config["output_dir"],"{patient}",  "{sample}", "08_zipper_raw", "{sample}.zipped.bam")
    params:
        reference=lambda wc: ref()
    threads: 4
    shell:
        r"""
        mkdir -p $(dirname {output.bam})
        java -Xmx4g -jar /usr/share/java/fgbio-2.2.1.jar ZipperBams \
          --unmapped {input.unmapped} \
          --input {input.aligned} \
          --output {output.bam} \
          --ref {params.reference}
        """


# 09 - Tri/index du BAM zipped raw
rule sort_zipped_raw_09:
    input:
        bam=rules.zipper_bams_raw_08.output.bam
    output:
        bam=os.path.join(config["output_dir"],"{patient}",  "{sample}", "09_sort_index_raw", "{sample}.zipped.coord.bam")
    threads: 4
    shell:
        r"""
        mkdir -p $(dirname {output.bam})
        samtools sort --threads {threads} -o {output.bam} {input.bam}
        """


rule index_zipped_raw_09:
    input:
        bam=rules.sort_zipped_raw_09.output.bam
    output:
        bai=os.path.join(config["output_dir"],"{patient}",  "{sample}", "09_sort_index_raw", "{sample}.zipped.coord.bam.bai")
    threads: 4
    shell:
        r"""
        samtools index -@ {threads} {input.bam}
        """


# 10 - GroupReadsByUmi
rule group_reads_by_umi_10:
    input:
        bam=rules.sort_zipped_raw_09.output.bam
    output:
        bam=os.path.join(config["output_dir"],"{patient}",  "{sample}", "10_grouped_by_umi", "{sample}.zipped.grouped.bam"),
        hist=os.path.join(config["output_dir"],"{patient}",  "{sample}", "10_grouped_by_umi", "{sample}.tag-family-sizes.txt")
    params:
        strategy=lambda wc: config.get("group_reads_by_umi", {}).get("strategy", "Adjacency"),
        edits=lambda wc: config.get("group_reads_by_umi", {}).get("edits", 1)
    threads: 8
    shell:
        r"""
        mkdir -p $(dirname {output.bam})
        java -Xmx8g -jar /usr/share/java/fgbio-2.2.1.jar GroupReadsByUmi \
          --input {input.bam} \
          --strategy {params.strategy} \
          --edits {params.edits} \
          --output {output.bam} \
          --family-size-histogram {output.hist}
        """


# 11 - Consensus non filtré
rule call_molecular_consensus_reads_11:
    input:
        bam=rules.group_reads_by_umi_10.output.bam
    output:
        bam=os.path.join(config["output_dir"],"{patient}",  "{sample}", "11_consensus_unfiltered", "{sample}.grouped.consensus.bam")
    params:
        min_reads=lambda wc: config.get("call_molecular_consensus_reads", {}).get("min_reads", 1),
        min_input_base_quality=lambda wc: config.get("call_molecular_consensus_reads", {}).get("min_input_base_quality", 30)
    threads: 4
    shell:
        r"""
        mkdir -p $(dirname {output.bam})
        java -Xmx4g -jar /usr/share/java/fgbio-2.2.1.jar CallMolecularConsensusReads \
          --input {input.bam} \
          --output {output.bam} \
          --min-reads {params.min_reads} \
          --min-input-base-quality {params.min_input_base_quality} \
          --threads {threads}
        """


# 12 - Filtrage consensus
rule filter_consensus_reads_12:
    input:
        bam=rules.call_molecular_consensus_reads_11.output.bam
    output:
        bam=os.path.join(config["output_dir"],"{patient}",  "{sample}", "12_consensus_filtered","{sample}.filtered_consensus.unmapped.bam")
    params:
        reference=lambda wc: ref(),
        min_reads=lambda wc: config.get("filter_consensus_reads", {}).get("min_reads", 3),
        min_base_quality=lambda wc: config.get("filter_consensus_reads", {}).get("min_base_quality", 30)
    threads: 4
    shell:
        r"""
        mkdir -p $(dirname {output.bam})
        java -Xmx4g -jar /usr/share/java/fgbio-2.2.1.jar FilterConsensusReads \
          --input {input.bam} \
          --output {output.bam} \
          --ref {params.reference} \
          --min-reads {params.min_reads} \
          --min-base-quality {params.min_base_quality}
        """


# 13 - Consensus BAM -> FASTQ
rule consensus_bam_to_fastq_13:
    input:
        bam=rules.filter_consensus_reads_12.output.bam
    output:
        fastq1=os.path.join(config["output_dir"],"{patient}",  "{sample}", "13_consensus_fastq", "{sample}.R1.consensus.fastq.gz"),
        fastq2=os.path.join(config["output_dir"],"{patient}",  "{sample}", "13_consensus_fastq", "{sample}.R2.consensus.fastq.gz")
    threads: 4
    shell:
        r"""
        mkdir -p $(dirname {output.fastq1})
        java -XX:ParallelGCThreads={threads} -Xmx8g -Xms8g -jar /usr/share/java/picard-2.8.1.jar SamToFastq \
          I={input.bam} \
          F={output.fastq1} \
          F2={output.fastq2}
        """


# 14 - Alignement consensus
rule align_consensus_fastq_14:
    input:
        fastq1=rules.consensus_bam_to_fastq_13.output.fastq1,
        fastq2=rules.consensus_bam_to_fastq_13.output.fastq2
    output:
        sam=os.path.join(config["output_dir"],"{patient}",  "{sample}","14_align_consensus","{sample}.consensus.aligned.sam")
    threads: 8
    params:
        rg=lambda wc: '@RG\\tID:{0}\\tSM:{1}\\tPL:{2}'.format(
            sample_cfg(wc).get("read_group_id", wc.sample),
            sample_cfg(wc).get("sample_name", wc.sample),
            sample_cfg(wc).get("platform", "ILLUMINA").upper(),
        ),
        reference=lambda wc: ref()
    shell:
        r"""
        mkdir -p $(dirname {output.sam})
        bwa-mem2 mem -t {threads} -Y \
          -R "{params.rg}" \
          {params.reference} \
          {input.fastq1} \
          {input.fastq2} \
          > {output.sam}
        """


# 15 - SAM -> BAM consensus
rule sam_to_bam_consensus_15:
    input:
        sam=rules.align_consensus_fastq_14.output.sam
    output:
        bam=os.path.join(config["output_dir"],"{patient}",  "{sample}", "15_consensus_sam_to_bam","{sample}.consensus.aligned.bam")
    threads: 8
    shell:
        r"""
        mkdir -p $(dirname {output.bam})
        samtools view -@ {threads} -bS {input.sam} > {output.bam}
        """


# 15b - Supprime secondary/supplementary du consensus aligné
rule remove_secondary_supplementary_consensus_15b:
    input:
        bam=rules.sam_to_bam_consensus_15.output.bam
    output:
        bam=os.path.join(config["output_dir"],"{patient}",  "{sample}","15b_consensus_primary_only", "{sample}.consensus.aligned.primary.bam")
    threads: 4
    shell:
        r"""
        mkdir -p $(dirname {output.bam})
        samtools view -@ {threads} -h -F 0x900 {input.bam} | samtools view -@ {threads} -b -o {output.bam}
        """


# 16 - Tri queryname consensus
rule sort_unmapped_queryname_consensus_16:
    input:
        bam=rules.filter_consensus_reads_12.output.bam
    output:
        bam=os.path.join(config["output_dir"],"{patient}",  "{sample}", "16_unmapped_queryname_sort_consensus","{sample}.unmapped.queryname.bam")
    threads: 4
    shell:
        r"""
        mkdir -p $(dirname {output.bam})
        java -Xmx4g -jar /usr/share/java/picard.jar SortSam \
          I={input.bam} \
          O={output.bam} \
          SORT_ORDER=queryname
        """


rule sort_aligned_queryname_consensus_16b:
    input:
        bam=rules.remove_secondary_supplementary_consensus_15b.output.bam
    output:
        bam=os.path.join(config["output_dir"],"{patient}",  "{sample}","16b_aligned_queryname_sort_consensus","{sample}.aligned.queryname.bam")
    shell:
        r"""
        mkdir -p $(dirname {output.bam})
        java -Xmx4g -jar /usr/share/java/picard.jar SortSam \
          I={input.bam} \
          O={output.bam} \
          SORT_ORDER=queryname
        """


# 17 - Retire @RG avant zipper consensus
rule remove_rg_header_consensus_17:
    input:
        bam=rules.sort_aligned_queryname_consensus_16b.output.bam
    output:
        header=os.path.join(config["output_dir"],"{patient}",  "{sample}","17_reheader_consensus", "{sample}.aligned.noRG.header.sam")
    threads: 1
    shell:
        r"""
        mkdir -p $(dirname {output.header})
        samtools view -H {input.bam} | grep -v '^@RG' > {output.header}
        """


rule reheader_aligned_queryname_consensus_17b:
    input:
        header=rules.remove_rg_header_consensus_17.output.header,
        bam=rules.sort_aligned_queryname_consensus_16b.output.bam
    output:
        bam=os.path.join(config["output_dir"],"{patient}",  "{sample}", "17b_reheader_consensus", "{sample}.aligned.queryname.noRG.bam")
    threads: 1
    shell:
        r"""
        mkdir -p $(dirname {output.bam})
        samtools reheader {input.header} {input.bam} > {output.bam}
        """


# 18- Zipper consensus
rule zipper_bams_consensus_18:
    input:
        unmapped=rules.sort_unmapped_queryname_consensus_16.output.bam,
        aligned=rules.reheader_aligned_queryname_consensus_17b.output.bam
    output:
        bam=os.path.join(config["output_dir"],"{patient}",  "{sample}","18_consensus_zipper","{sample}.consensus.zipped.bam")
    params:
        reference=lambda wc: ref(),
        tags_to_reverse=lambda wc: config.get("Zipper_Bams", {}).get("tags-to-reverse", "Consensus"),
        tags_to_revcomp=lambda wc: config.get("Zipper_Bams", {}).get("tags-to-revcomp", "Consensus")
    threads: 4
    shell:
        r"""
        mkdir -p $(dirname {output.bam})
        java -Xmx4g -jar /usr/share/java/fgbio-2.2.1.jar ZipperBams \
          --unmapped {input.unmapped} \
          --input {input.aligned} \
          --output {output.bam} \
          --ref {params.reference} \
          --tags-to-reverse {params.tags_to_reverse} \
          --tags-to-revcomp {params.tags_to_revcomp}
        """


# 19 - Tri final
rule sort_post_zipper_consensus_19:
    input:
        bam=rules.zipper_bams_consensus_18.output.bam
    output:
        bam=os.path.join(config["output_dir"],"{patient}",  "{sample}","19_consensus_final","{sample}.consensus.zipped.coord.bam")
    threads: 4
    shell:
        r"""
        mkdir -p $(dirname {output.bam})
        samtools sort --threads {threads} -o {output.bam} {input.bam}
        """



rule index_post_zipper_consensus_19:
    input:
        bam=rules.sort_post_zipper_consensus_19.output.bam
    output:
        bai=os.path.join(config["output_dir"],"{patient}",  "{sample}","19_consensus_final","{sample}.consensus.zipped.coord.bam.bai")
    threads: 4
    shell:
        r"""
        samtools index -@ {threads} {input.bam}
        """


#Variant_calling_Mutect2_01
rule Variant_calling_Mutect2_01:
    input:
        BAM_T=lambda wc: consensus_bam(wc.patient, wc.tumor),
        BAM_N=lambda wc: consensus_bam(wc.patient, wc.normal),
        BAI_T=lambda wc: consensus_bai(wc.patient, wc.tumor),
        BAI_N=lambda wc: consensus_bai(wc.patient, wc.normal)
    output:
        vcf=os.path.join(config["output_dir"], "{patient}","{normal}_vs_{tumor}", "Variant_calling_Mutect2_01","{normal}_vs_{tumor}.somatic.vcf.gz"),
        tar=os.path.join(config["output_dir"], "{patient}","{normal}_vs_{tumor}", "Variant_calling_Mutect2_01","{normal}_vs_{tumor}.f1r2.tar.gz"),
        log=os.path.join(config["output_dir"], "{patient}","{normal}_vs_{tumor}","log_mutect2.log")
    threads: 4
    params:
        thread=lambda wc: config.get("mutect2_calling", {}).get("native_pair_hmm_threads", 4),
        min_base_quality=lambda wc: config.get("mutect2_calling", {}).get("min_base_quality_score", 10),
        f1r2_min_bq =lambda wc: config.get("mutect2_calling", {}).get("f1r2_min_bq", 10),
        normal=lambda wc: normal_cfg(wc)["sample_name"],
        tumor=lambda wc: tumor_cfg(wc)["sample_name"],
        ref=lambda wc: ref()

    shell:
        r"""
        mkdir -p $(dirname {output.vcf})
        java -jar /usr/share/java/GenomeAnalysisTK.jar  Mutect2 \
          -R {params.ref} \
          -I:tumor {input.BAM_T} \
          -I:normal {input.BAM_N} \
          -normal {params.normal} \
          -tumor {params.tumor} \
          -O {output.vcf} \
          --native-pair-hmm-threads {params.thread} \
          --f1r2-tar-gz  {output.tar} \
          --min-base-quality-score {params.min_base_quality} \
          --f1r2-min-bq {params.f1r2_min_bq} \
          2>&1 | tee -a {output.log}        
        """
        


#Learn_ReadOrientation_Model_02
rule Learn_ReadOrientation_Model_02:
    input:
        tar=rules.Variant_calling_Mutect2_01.output.tar
    output:
        tar=os.path.join(config["output_dir"], "{patient}","{normal}_vs_{tumor}", "Learn_ReadOrientation_Model_02","{normal}_vs_{tumor}.read-orientation-model.tar.gz")
    threads: 1
    shell:
        r"""
        mkdir -p $(dirname {output.tar})
        java -jar /usr/share/java/GenomeAnalysisTK.jar LearnReadOrientationModel \
        -I {input.tar} \
        -O {output.tar}
        """


#GetPileupSummaries_Tumor_03
rule Get_Pileup_Summaries_Tumor_03:
    input:
        BAM_T=lambda wc: consensus_bam(wc.patient, wc.tumor),
        BAI_T=lambda wc: consensus_bai(wc.patient, wc.tumor)

    output:
        pileup=os.path.join(config["output_dir"], "{patient}","{normal}_vs_{tumor}", "Get_Pileup_Summaries_03","{normal}_vs_{tumor}.tumor.pileups.table")
    
    params:
        gnome=lambda wc: ref_gnomad(),
        bed=lambda wc: ref_bed(),
        ref=lambda wc: ref()
    threads: 2
    shell:
        r"""
        mkdir -p $(dirname {output.pileup})
        java -jar /usr/share/java/GenomeAnalysisTK.jar GetPileupSummaries \
          -R {params.ref} \
          -I {input.BAM_T} \
          -V {params.gnome} \
          -L {params.bed} \
          -O {output.pileup}
        """


#GetPileupSummaries_Normal_04
rule Get_Pileup_Summaries_Normal_04:
    input:
        BAM_N=lambda wc: consensus_bam(wc.patient, wc.normal),
        BAI_N=lambda wc: consensus_bai(wc.patient, wc.normal)

    output:
        pileup=os.path.join(config["output_dir"], "{patient}","{normal}_vs_{tumor}", "Get_Pileup_Summaries_04","{normal}_vs_{tumor}.normal.pileups.table")
    
    params:
        gnome=lambda wc: ref_gnomad(),
        bed=lambda wc: ref_bed(),
        ref=lambda wc: ref()
    threads: 2
    shell:
        r"""
        mkdir -p $(dirname {output.pileup})
        java -jar /usr/share/java/GenomeAnalysisTK.jar GetPileupSummaries \
          -R {params.ref} \
          -I {input.BAM_N} \
          -V {params.gnome} \
          -L {params.bed} \
          -O {output.pileup}
        """

#Calculate_Contamination_05
rule Calculate_Contamination_05:
    input:
        tumor=rules.Get_Pileup_Summaries_Tumor_03.output.pileup,
        matched_normal=rules.Get_Pileup_Summaries_Normal_04.output.pileup
    output:
        segments=os.path.join(config["output_dir"], "{patient}","{normal}_vs_{tumor}", "Calculate_Contamination_05","{normal}_vs_{tumor}.tumor.segments.table"),
        contamination=os.path.join(config["output_dir"], "{patient}","{normal}_vs_{tumor}","Calculate_Contamination_05","{normal}_vs_{tumor}.tumor.contamination.table")
    threads: 1
    shell:
        r"""
        mkdir -p $(dirname {output.contamination})
        java -jar /usr/share/java/GenomeAnalysisTK.jar CalculateContamination \
         -I {input.tumor} \
         -matched {input.matched_normal} \
         --tumor-segmentation {output.segments} \
         -O {output.contamination} 
        """

#Filter_Mutect_Calls_06
rule Filter_Mutect_Calls_06:
    input:
        vcf=rules.Variant_calling_Mutect2_01.output.vcf,
        prior=rules.Learn_ReadOrientation_Model_02.output.tar,
        contamination=rules.Calculate_Contamination_05.output.contamination,
        segments=rules.Calculate_Contamination_05.output.segments
    output:
        vcf=os.path.join(config["output_dir"], "{patient}","{normal}_vs_{tumor}","Filter_Mutect_Calls_06","{normal}_vs_{tumor}.filtered.vcf.gz")
    params:
        ref=lambda wc: ref()
    threads: 1
    shell:
        r"""
        mkdir -p $(dirname {output.vcf})
        java -jar /usr/share/java/GenomeAnalysisTK.jar FilterMutectCalls -R {params.ref} -V {input.vcf} --ob-priors {input.prior} --contamination-table {input.contamination} --tumor-segmentation {input.segments} -O {output.vcf}

        """

#Select_variants_07
rule Select_variants_07:
    input:
        vcf=rules.Filter_Mutect_Calls_06.output.vcf

    output:
        vcf=os.path.join(config["output_dir"], "{patient}","{normal}_vs_{tumor}","Select_variants_07","{normal}_vs_{tumor}.PASS.vcf.gz"),
        done=touch(os.path.join(config["output_dir"], "{patient}","{normal}_vs_{tumor}","Select_variants_07","done.py"))

    params:
        ex=lambda wc: config.get("Select_variants", {}).get("exclude-filtered","true"),
        ref=lambda wc: ref()

    threads: 1

    shell:
        r"""
        mkdir -p $(dirname {output.vcf})
        java -jar /usr/share/java/GenomeAnalysisTK.jar SelectVariants \
          -R {params.ref} \
          -V {input.vcf} \
          --exclude-filtered {params.ex} \
          -O {output.vcf}
        
        
        """

#Msisensor_score_08
rule Msisensor_score_08:
    input:
        BAM_T=lambda wc: consensus_bam(wc.patient, wc.tumor),
        BAM_N=lambda wc: consensus_bam(wc.patient, wc.normal),
        done=rules.Select_variants_07.output.done

    output:
        prefix=os.path.join(config["output_dir"], "{patient}","{normal}_vs_{tumor}", "Msisensor_score_08", "{normal}_vs_{tumor}.somatic.prefix"),
        done=touch(os.path.join(config["output_dir"], "{patient}","{normal}_vs_{tumor}","Msisensor_score_08","done.py"))

    params:
        microsatellite=lambda wc: ref_satellite_list(),
        bed=lambda wc: ref_bed()

    threads: 2

    shell:
        r"""
        mkdir -p $(dirname {output.prefix})
        msisensor2 msi -d {params.microsatellite}  -n {input.BAM_N} -t {input.BAM_T} -b {threads} -o {output.prefix}  -e {params.bed}

       
        """



#CNV_batch_calling_09
rule CNV_batch_calling_09:
    input:
        BAM_T=lambda wc: consensus_bam(wc.patient, wc.tumor),
        BAM_N=lambda wc: consensus_bam(wc.patient, wc.normal),
        done=rules.Msisensor_score_08.output.done

    output:
        DIR=os.path.join(config["output_dir"], "{patient}","{normal}_vs_{tumor}","CNV_batch_calling_09","CNVkit_outputs"),
        done=touch(os.path.join(config["output_dir"], "{patient}","{normal}_vs_{tumor}","CNV_batch_calling_09","done.py"))

    params:
        ref=lambda wc: ref(),
        bed=lambda wc: ref_bed(),
        accessible_bed=lambda wc: ref_accessible_bed(),
        gene_labels=lambda wc: ref_gene_labels()

    threads:2
    shell:
        r"""
        mkdir -p $(dirname {output.done})
        cnvkit batch {input.BAM_T} \
         -n {input.BAM_N} \
         -t {params.bed} \
         -f {params.ref} \
         --annotate {params.gene_labels} \
         --access {params.accessible_bed} \
         -p {threads} \
         -d {output.DIR} \
         --scatter \
         --diagram 

        
        """


#CNV_depth_coverage__10
rule CNV_depth_coverage_10:
    input:
        BAM_T=lambda wc: consensus_bam(wc.patient, wc.tumor),
        bed=lambda wc: ref_bed(),
        

    params:
        ref=lambda wc: ref(),
        tumor=lambda wc: tumor_cfg(wc)["sample_name"]

    threads: 8
    output:
        dir_coverage=os.path.join(config["output_dir"], "{patient}","{tumor}","CNV_depth_coverage_10"),
        done=os.path.join(config["output_dir"], "{patient}","{tumor}","CNV_depth_coverage_10","done.py")

    shell:
        r"""
        mkdir -p $(dirname {output.done})
        cd {output.dir_coverage} 
        java -jar /usr/share/java/GenomeAnalysisTK.jar DepthOfCoverage \
         -R {params.ref} \
         -O {params.tumor} \
         -I {input.BAM_T} \
         -L {input.bed} 
        touch {output.done}

        """



#viscap_cnv_calling_11 
rule viscap_cnv_calling_11:
    input:
        coverage=lambda wc : [os.path.join(config["output_dir"],v["patient"],v['tumor'],"CNV_depth_coverage_10") for v in VALID_COMPARISONS if v["patient"] == wc.patient],
        bed=lambda wc : ref_annotated_bed()
        

    output:
        done=os.path.join(config["output_dir"],"{patient}","viscap_cnv_calling_11","done.py"),
        output_dir=os.path.join(config["output_dir"],"{patient}","viscap_cnv_calling_11","results")
        

    threads:4
    log:
        viscap=os.path.join(config["output_dir"],"{patient}","logs","viscap.log")

    params:
        viscap_script=config["viscap"]["script"],
        run_dir=os.path.join(config["output_dir"],"{patient}","viscap_cnv_calling_11","run_inputs"),
        coverage_dir=os.path.join(config["output_dir"],"{patient}","viscap_cnv_calling_11","run_inputs","coverage"),
        interval_dir=os.path.join(config["output_dir"],"{patient}","viscap_cnv_calling_11","run_inputs","intervals"),
        cfg=config["viscap"]["cfg"]

    shell:
        r"""
        set -euo pipefail

        
        mkdir -p {params.coverage_dir}
        mkdir -p {params.interval_dir}
        mkdir -p {output.output_dir}
        

for covdir in {input.coverage}; do
    for cov in "$covdir"/*.sample_interval_summary; do
        [ -e "$cov" ] || continue

        awk -F "," 'BEGIN{{OFS="\t"}} {{print $1,$4}}' "$cov" \
        | sed 's/^chr//' \
        > "{params.coverage_dir}/$(basename "$cov")"
    done
done
       
            cp "$(realpath {input.bed})" "{params.interval_dir}/targets.bed"

sed 's|^[[:space:]]*interval_list_dir[[:space:]]*<-[[:space:]]*.*|interval_list_dir <- "{params.interval_dir}"|' \
    "{params.cfg}" > "{params.run_dir}/VisCap.cfg"

        cp {params.viscap_script} {params.run_dir}

        cd {params.run_dir}

        Rscript   $(basename {params.viscap_script}) {params.coverage_dir} {output.output_dir} {params.interval_dir} > {log.viscap} 2>&1 
        touch {output.done}
        """

#lookup_and_export_viscap_xls_to_bed_12
rule lookup_and_export_viscap_xls_to_bed_12:
    input:
        done=rules.viscap_cnv_calling_11.output.done,
        cnv_dir=rules.viscap_cnv_calling_11.output.output_dir
        
        
    output:
        cnv_bed=os.path.join(config["output_dir"],"{patient}","{tumor}","lookup_and_export_viscap_xls_to_bed_12","{tumor}.tab.cnvs.bed"),
        done=os.path.join(config["output_dir"],"{patient}","{tumor}","lookup_and_export_viscap_xls_to_bed_12","done.py")
        

    params:
        script=config["viscap"]["tobed_script"],
        sample_name=lambda wc: tumor_cfg(wc)["sample_name"]

    threads:1
    shell:
      r"""
      python3 {params.script} \
      --cnv_dir {input.cnv_dir} \
      --cnv_bed {output.cnv_bed}  \
      --sample_name {params.sample_name}
      touch {output.done}
      """


#Annot_SV_CNV_viscap_annotation_13
rule Annot_SV_CNV_viscap_annotation_13:
    input:
        cnv_bed=rules.lookup_and_export_viscap_xls_to_bed_12.output.cnv_bed,
        done=rules.lookup_and_export_viscap_xls_to_bed_12.output.done
    output:
        tsv_annotated=os.path.join(config["output_dir"], "{patient}","{tumor}","Annot_SV_CNV_viscap_annotation_13","{tumor}.viscapCNV.annotated.tsv"),
        done=os.path.join(config["output_dir"],"{patient}","{tumor}","Annot_SV_CNV_viscap_annotation_13","done.py")

    params:
        prefix=os.path.join(config["output_dir"], "{patient}","{tumor}","Annot_SV_CNV_viscap_annotation_13","{tumor}.viscapCNV.annotated"),
        genome_build=lambda wc: config.get("annotSV", {}).get("genomeBuild"),
        SVminSize=lambda wc: config.get("annotSV", {}).get("SVminSize"),
        annotation_mode=lambda wc: config.get("annotSV", {}).get("annotationMode"),
        annotation_dir=lambda wc: config.get("annotSV", {}).get("annotationDir"),
        samplesidBEDcol=lambda wc: config.get("annotSV", {}).get("samplesidBEDcol"),
        svtBEDcol=lambda wc: config.get("annotSV", {}).get("svtBEDcol")

    threads:1
    run:
        cnv_path = input.cnv_bed
        os.makedirs(os.path.dirname(params.prefix), exist_ok=True)
        if cnv_path:
            shell(
                "AnnotSV "
                "-SVinputFile {input.cnv_bed} "
                "-samplesidBEDcol {params.samplesidBEDcol} "
                "-svtBEDcol {params.svtBEDcol} "
                "-genomeBuild {params.genome_build} "
                "-outputFile {params.prefix} "
                "-SVminSize {params.SVminSize} "
                "-annotationsDir {params.annotation_dir} "
                "-annotationMode {params.annotation_mode} "
            
            )
            shell ("sed -i '1s/user#1/gene_input/; 1s/user#2/log2/' {output.tsv_annotated}")
        else:
          
            shell("touch {output.tsv_annotated}")

        shell("touch {output.done}")


#SV_workflow_14
rule SV_workflow_14:
    input:
        BAM_T=lambda wc: consensus_bam(wc.patient, wc.tumor),
        BAM_N=lambda wc: consensus_bam(wc.patient, wc.normal),
        done=rules.CNV_batch_calling_09.output.done
    output:
        done=touch(os.path.join(config["output_dir"], "{patient}","{normal}_vs_{tumor}","SV_workflow_14","done.py")),
        DIR=os.path.join(config["output_dir"], "{patient}","{normal}_vs_{tumor}","SV_workflow_14")
    params:
        ref=lambda wc: ref()
    threads:2
    shell:
        r"""
        
        mkdir -p $(dirname {output.done})
        configManta.py --normalBam {input.BAM_N} \
         --tumorBam {input.BAM_T}  \
         --referenceFasta {params.ref} \
         --runDir {output.DIR} 

        """

#SV_fusion_run_15
rule SV_run_workflow_15:
    input:
        done=rules.SV_workflow_14.output.done,
        dir=rules.SV_workflow_14.output.DIR
    output:
        done=os.path.join(config["output_dir"],'{patient}',"{normal}_vs_{tumor}","SV_run_workflow_15","done.py")

    threads:1
    
    shell:
        r"""
        mkdir -p $(dirname {output.done})
        cd {input.dir} 
        python runWorkflow.py
        touch {output.done}

        """




#Annovar_Snp_and_indel_annotation_16
rule Annovar_Snp_and_indel_annotation_16:
    input:
        vcf=rules.Select_variants_07.output.vcf,
        done=rules.SV_run_workflow_15.output.done
    output:
        done=os.path.join(config["output_dir"],"{patient}","{normal}_vs_{tumor}","Annovar_Snp_and_indel_annotation_16","done.py"),
        uncompressed_vcf=os.path.join(config["output_dir"], "{patient}","{normal}_vs_{tumor}","Annovar_Snp_and_indel_annotation_16","{normal}_vs_{tumor}.PASS.vcf"),
        multianno=os.path.join(config["output_dir"], "{patient}","{normal}_vs_{tumor}","Annovar_Snp_and_indel_annotation_16","{normal}_vs_{tumor}.hg38_multianno.txt"),
        vcf_multianno=os.path.join(config["output_dir"], "{patient}","{normal}_vs_{tumor}","Annovar_Snp_and_indel_annotation_16","{normal}_vs_{tumor}.hg38_multianno.vcf")
    params:
        prefix=os.path.join(config["output_dir"], "{patient}","{normal}_vs_{tumor}","Annovar_Snp_and_indel_annotation_16","{normal}_vs_{tumor}"),
        annovar=lambda wc: config.get("annovar", {}).get("annovar_database"),
        buildver=lambda wc: config.get("annovar", {}).get("buildver"),
        protocol=lambda wc: config.get("annovar", {}).get("protocol"),
        operation=lambda wc: config.get("annovar", {}).get("operation"),
        nastring=lambda wc: config.get("annovar", {}).get("nastring")

    shell:
        r"""
        mkdir -p $(dirname {params.prefix})
        gzip -kdc {input.vcf} > {output.uncompressed_vcf} 
        table_annovar.pl {output.uncompressed_vcf} {params.annovar} \
       -buildver {params.buildver} \
       -out {params.prefix} \
       -remove \
       -protocol {params.protocol} \
       -operation {params.operation} \
       -nastring {params.nastring} \
       -vcfinput 
        touch {output.done}

        """

#snp_Eff_annotation_17
rule snp_Eff_annotation_17:
    input:
        vcf=rules.Select_variants_07.output.vcf,
        done=rules.Annovar_Snp_and_indel_annotation_16.output.done
    output:
        done=os.path.join(config["output_dir"],"{patient}","{normal}_vs_{tumor}","snp_Eff_annotation_17","done.py"),
        vcf=os.path.join(config["output_dir"],"{patient}","{normal}_vs_{tumor}","snp_Eff_annotation_17","{normal}_vs_{tumor}.ann_snpEff.vcf"),
        stats=os.path.join(config["output_dir"],"{patient}","{normal}_vs_{tumor}","snp_Eff_annotation_17","{normal}_vs_{tumor}.stats.html")
        

    params:
        genome_version=lambda wc: config.get("snp_Eff", {}).get("genome_version"),
        snp_Eff=lambda wc: config.get("snp_Eff", {}).get("config_database")

    shell:
        r"""
        mkdir -p $(dirname {output.vcf})
        snpEff  \
          -lof \
          -cancer \
          -c {params.snp_Eff} \
          -stats {output.stats} \
          {params.genome_version} \
          {input.vcf} | sed -n '/^##fileformat=/,$p' > {output.vcf}
        touch {output.done}

        """

#Snp_sift_annotation_18
rule snp_sift_annotation_18:
    input:
        vcf=rules.snp_Eff_annotation_17.output.vcf,
        done=rules.snp_Eff_annotation_17.output.done
    output:
        vcf=os.path.join(config["output_dir"], "{patient}","{normal}_vs_{tumor}","snp_sift_annotation_18","{normal}_vs_{tumor}.ann_snpEff_snpSift.vcf"),
        done=os.path.join(config["output_dir"],"{patient}","{normal}_vs_{tumor}","snp_sift_annotation_18","done.py")

    params:
        clinvar_name=recherche_databse("snp_Sift","database",config,"clinvar","name"),
        clinvar_path=recherche_databse("snp_Sift","database",config,"clinvar","path"),
        dbsnp_name=recherche_databse("snp_Sift","database",config,"dbsnp","name"),
        dbsnp_path=recherche_databse("snp_Sift","database",config,"dbsnp","path"),
        cosmic_name=recherche_databse("snp_Sift","database",config,"cosmic","name"),
        cosmic_path=recherche_databse("snp_Sift","database",config,"cosmic","path")


    shell:
       r"""
       mkdir -p $(dirname {output.vcf})
       snpSift Annotate -name {params.clinvar_name} {params.clinvar_path} {input.vcf} \
       | snpSift Annotate -name {params.dbsnp_name} {params.dbsnp_path}  \
         | snpSift Annotate -name {params.cosmic_name} {params.cosmic_path} > {output.vcf}
    
       touch {output.done}

       """

#AnnotSV_SV_annotation_19
rule AnnotSV_SV_annotation_19:
    input:
        vcf_dir=rules.SV_workflow_14.output.DIR,
        done=rules.snp_sift_annotation_18.output.done


    output:
        tsv_annotated=os.path.join(config["output_dir"], "{patient}","{normal}_vs_{tumor}","AnnotSV_SV_annotation_19","{normal}_vs_{tumor}.somaticSV.annotated.tsv"),
        done=os.path.join(config["output_dir"],"{patient}","{normal}_vs_{tumor}","AnnotSV_SV_annotation_19","done.py")

    params:
        prefix=os.path.join(config["output_dir"], "{patient}","{normal}_vs_{tumor}","AnnotSV_SV_annotation_19","{normal}_vs_{tumor}.somaticSV.annotated"),
        genome_build=lambda wc: config.get("annotSV", {}).get("genomeBuild"),
        SVminSize=lambda wc: config.get("annotSV", {}).get("SVminSize"),
        annotation_mode=lambda wc: config.get("annotSV", {}).get("annotationMode"),
        annotation_dir=lambda wc: config.get("annotSV", {}).get("annotationDir")
        

    run:
        vcf_path = os.path.join(input.vcf_dir, "results", "variants", "somaticSV.vcf.gz")
        os.makedirs(os.path.dirname(params.prefix), exist_ok=True)
        if vcf_has_variants(vcf_path):
            shell(
                "AnnotSV "
                "-SVinputFile {vcf_path} "
                "-genomeBuild {params.genome_build} "
                "-outputFile {params.prefix} "
                "-SVminSize {params.SVminSize} "
                "-annotationsDir {params.annotation_dir} "
                "-annotationMode {params.annotation_mode} "
            )
        else:
  
            shell("touch {output.tsv_annotated}")

        shell("touch {output.done}")


#CNVkit_export_cns_to_vcf_20
rule CNVkit_export_cns_to_vcf_20:
    input:
        cns_dir=rules.CNV_batch_calling_09.output.DIR,
        done=rules.AnnotSV_SV_annotation_19.output.done

    output:
        vcf=os.path.join(config["output_dir"], "{patient}","{normal}_vs_{tumor}","CNVkit_export_cns_to_vcf_20","{normal}_vs_{tumor}.cns.vcf"),
        done=os.path.join(config["output_dir"],"{patient}","{normal}_vs_{tumor}","CNVkit_export_cns_to_vcf_20","done.py")

    params:
        tumor=lambda wc: tumor_cfg(wc)["sample_name"]

    shell:
       r"""
       mkdir -p $(dirname {output.vcf})
       cnvkit export vcf {input.cns_dir}/{params.tumor}.consensus.zipped.coord.call.cns -i {params.tumor} -o  {output.vcf}
       touch {output.done}

        """
#Annot_SV_CNV_annotation_21
rule Annot_SV_CNV_annotation_21:
    input:
        vcf=rules.CNVkit_export_cns_to_vcf_20.output.vcf,
        done=rules.CNVkit_export_cns_to_vcf_20.output.done
    output:
        done=os.path.join(config["output_dir"],"{patient}","{normal}_vs_{tumor}","Annot_SV_CNV_annotation_21","done.py"),
        tsv_annotated=os.path.join(config["output_dir"], "{patient}","{normal}_vs_{tumor}","Annot_SV_CNV_annotation_21","{normal}_vs_{tumor}.cns.vcf.annotated.tsv")
    params:
        prefix=os.path.join(config["output_dir"], "{patient}","{normal}_vs_{tumor}","Annot_SV_CNV_annotation_21","{normal}_vs_{tumor}.cns.vcf.annotated"),
        genome_build=lambda wc: config.get("annotSV", {}).get("genomeBuild"),
        SVminSize=lambda wc: config.get("annotSV", {}).get("SVminSize"),
        annotation_mode=lambda wc: config.get("annotSV", {}).get("annotationMode"),
        annotation_dir=lambda wc: config.get("annotSV", {}).get("annotationDir")


    run:
        vcf_path = input.vcf
        os.makedirs(os.path.dirname(params.prefix), exist_ok=True)
        if vcf_has_variants(vcf_path):
            shell(
                "AnnotSV "
                "-SVinputFile {input.vcf} "
                "-genomeBuild {params.genome_build} "
                "-outputFile {params.prefix} "
                "-SVminSize {params.SVminSize} "
                "-annotationsDir {params.annotation_dir} "
                "-annotationMode {params.annotation_mode}"
            )
        else:

            shell("touch {output.tsv_annotated}")

        shell("touch {output.done}")




