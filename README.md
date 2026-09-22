
#  UMI aware consensus pipeline

Le pipeline UMI-aware consensus permet la détéction et la quantification des altérations somatiques au cours du temps.
Dans la permière partie du pipeline, les fastqs des échantillons tumeur et normal vont être analysés en série pour génerer les consensus BAMS à partir des UMIs.

Dans un deuxième temps, les consensus bams des échantillons tumoraux sont intégrés simultanément avec le consensus bam de l'échantillon normal pour faire tourner les outils de variant calling en mode Tumeur versus Normal dans le but d'analyser 5 types de variants candidats :

  - les SNP et les indels
  - les altération du nombre de copies de variants
  - les variants structuraux
  - les instabilités microsatellite

En dernier lieu, les sorties des outils de variant calling vont être annotées avec les outils d'annotation adaptés au format de sortie .


## Description des principales étapes du pipeline

### Traitement des UMIs


| Étape | Rule | Dossier de sortie| Fichier produit| Rôle|
|-------|--------|-------|-------|-------|
| 0     | input initial                            | libre, défini dans le JSON                    | FASTQ R1/R2                                                    | données brutes                                        |
| 1     | fastq_to_ubam                            | 01_fastq_to_ubam/                             | {sample}.unmapped.bam                                          | créer un uBAM contenant UMI + métadonnées             |
| 2     | ubam_to_fastq                            | 02_ubam_fastq/R1/ et 02_ubam_fastq/R2/        | {sample}.R1.fastq.gz, {sample}.R2.fastq.gz                     | reconvertir le uBAM en FASTQ pour l’alignement        |
| 3     | align_raw_fastq                          | 03_align_raw/                                 | {sample}.aligned.sam                                           | aligner les reads bruts                               |
| 4     | sam_to_bam_raw                           | 04_sam_to_bam_raw/                            | {sample}.aligned.bam                                           | compresser le SAM en BAM                              |
| 5     | remove_secondary_supplementary           | 05_primary_only/                              | {sample}.aligned.primary.bam                                   | ne garder que les alignements principaux              |
| 6     | sort_unmapped_queryname                  | 06_queryname_sort/                            | {sample}.unmapped.queryname.bam                                | trier le BAM non aligné par nom de lecture            |
| 6b    | sort_primary_queryname                   | 06_queryname_sort/                            | {sample}.aligned.primary.queryname.bam                         | trier le BAM aligné primaire par nom de lecture       |
| 7     | remove_rg_header                         | 07_reheader/                                  | {sample}.aligned.noRG.header.sam                               | produire un header sans @RG                           |
| 7b    | reheader_primary_queryname               | 07_reheader/                                  | {sample}.aligned.primary.queryname.noRG.bam                    | réinjecter ce header dans le BAM aligné               |
| 8     | zipper_bams_raw                          | 08_zipper_raw/                                | {sample}.zipped.bam                                            | fusionner métadonnées uBAM + coordonnées d’alignement |
| 9     | sort_zipped_raw                          | 09_sort_index_raw/                            | {sample}.zipped.coord.bam                                      | trier le BAM zippé par coordonnées                    |
|       | index_zipped_raw                         | 09_sort_index_raw/                            | {sample}.zipped.coord.bam.bai                                  | indexer le BAM trié                                   |
| 10    | group_reads_by_umi                       | 10_grouped_by_umi/                            | {sample}.zipped.grouped.bam + {sample}.tag-family-sizes.txt    | regrouper les reads par famille UMI                   |
| 11    | call_molecular_consensus_reads           | 11_consensus_unfiltered/                      | {sample}.grouped.consensus.bam                                 | construire les reads consensus                        |
| 12    | filter_consensus_reads                   | 12_consensus_filtered/                        | {sample}.filtered_consensus.unmapped.bam                       | filtrer les consensus insuffisants                    |
| 13    | consensus_bam_to_fastq                   | 13_consensus_fastq/                           | {sample}.R1.consensus.fastq.gz, {sample}.R2.consensus.fastq.gz | reconvertir les consensus non alignés en FASTQ        |
| 14    | align_consensus_fastq                    | 14_align_consensus/                           | {sample}.consensus.aligned.sam                                 | réaligner les reads consensus                         |
| 15    | sam_to_bam_consensus                     | 15_consensus_sam_to_bam/                      | {sample}.consensus.aligned.bam                                 | compresser le SAM consensus en BAM                    |
| 15b   | remove_secondary_supplementary_consensus | 15b_remove_secondary_supplementary_consensus/ | {sample}.consensus.aligned.primary.bam                         | ne garder que les alignements principaux              |
| 16    | sort_unmapped_queryname_consensus_16     | 16_unmapped_queryname_sort_consensus /        | {sample}.unmapped,queryname.bam                                | trier le BAM non aligné par nom de lecture            |
| 16b   | sort_aligned_queryname_consensus_16b     | 16b_sort_aligned_queryname_consensus/         | {sample}.aligned.queryname.bam                                 | trier le BAM aligné primaire par nom de lecture       |
| 17    | remove_rg_header_consensus               | 17_remove_rg_header_consensus/                | {sample},aligned.noRG.header.sam                               | produire un header sans @RG                           |
| 17b   | reheader_aligned_queryname_consensus     | 17b_reheader_aligned_queryname_consensus/     | {sample},aligned.noRG.header.bam                               | réinjecter ce header dans le BAM aligné               |
| 18    | zipper_bams_consensus                    | 18_consensus_zipper/                          | {sample}.consensus.zipped.bam                                  | fusionner consensus non aligné + consensus aligné     |
| 19    | sort_consensus_final                     | 19_consensus_final                            | {sample}.consensus.zipped.coord.bam                            | trier le BAM final                                    |
| 19    | index_consensus_final                    | 19_consensus_final                            | {sample}.consensus.zipped.coord.bam.bai                        | indexer le BAM final                                  |



### Variant calling des SNP,Indel,CNV,SV,Instabilités microsatellites





### Annotation 



## Dossiers de srotie

- Par échantillon unique (normal ou tumeur) :

```

patient1
    ├── Tumeur/Normal
    │   ├── 01_fastq_to_ubam 
    │   ├── 02_ubam_fastq 
    │   ├── 03_align_raw 
    │   ├── 04_sam_to_bam_raw 
    │   ├── 05_primary_only 
    │   ├── 06_queryname_sort 
    │   ├── 07_reheader 
    │   ├── 08_zipper_raw 
    │   ├── 09_sort_index_raw 
    │   ├── 10_grouped_by_umi 
    │   ├── 11_consensus_unfiltered 
    │   ├── 12_consensus_filtered 
    │   ├── 13_consensus_fastq 
    │   ├── 14_align_consensus 
    │   ├── 15_consensus_sam_to_bam 
    │   ├── 15b_consensus_primary_only 
    │   ├── 16_unmapped_queryname_sort_consensus 
    │   ├── 16b_aligned_queryname_sort_consensus 
    │   ├── 17_reheader_consensus 
    │   ├── 17b_reheader_consensus 
    │   ├── 18_consensus_zipper 
    │   ├── 19_consensus_final 

```

- Par comparaison binaire entre échantillon tumeur vs échantillon Normal :

```

patient1
    ├── Normal_vs_Tumeur
    │   ├── Annot_SV_annotation_19 
    │   ├── AnnotSV_CNV_annotation_21 
    │   ├── Annovar_Snp_and_indel_annotation_16 
    │   ├── CNV_batch_calling_09 
    │   ├── CNVkit_export_cns_to_vcf_20 
    │   ├── Calculate_Contamination_05 
    │   ├── Filter_Mutect_Calls_06 
    │   ├── Get_Pileup_Summaries_03 
    │   ├── Get_Pileup_Summaries_04 
    │   ├── Learn_ReadOrientation_Model_02 
    │   ├── Msisensor_score_08 
    │   ├── SV_run_workflow_15 
    │   ├── SV_workflow_14 
    │   ├── Select_variants_07 
    │   ├── Variant_calling_Mutect2_01 
    │   ├── log_mutect2.log
    │   ├── snp_Eff_annotation_17
    │   └── snp_sift_annotation_18

```


- Par échantillon tumeur uniquement :

```

patient1
    ├── Tumeur
          ├── Annot_SV_CNV_viscap_annotation_13
          ├── CNV_depth_coverage_10
          └── lookup_and_export_viscap_xls_to_bed_12

```

- Par batch d'échantillon tumeur :

```

patient1
    ├── Tumeur 
    └── viscap_cnv_calling_11

```

## Outils de traitement des UMIs

Fgbio

## Outils d'alignement

Bwa-mem 2

## Outils d'annotation

ANNOVAR

SnpEFF

SnpSift

AnnotSV

## Outils de variant calling

Mutect2

CNVkit

VisCap

MSISensor

Manta

## Bases de données d'annotation

Cosmic

Clinvar

dbSNP

RefGene

Base de données Amazon GRCh38 de snpEFF



## Génerer le samplesheet à partir d'un formulaire

Utiliser le formulaire samplesheet.html afin de génerer le config file

## Génerer le config file

json generator permet de compléter un config file avec les paramètres adéquats  en lui spécifiant comme arguments :

- une sample sheet génerée à partir du formulaire_samplesheet  (--sample_sheet )

- un fichier json en entrée (--json)

- un fichier json de sortie (--output)

```{python}

python3 scripts/build_patient_config.py --samplesheet files/sample_sheet_Patient_001.tsv --json files/test_patients.json --output outdir/config_patients_001.json

```

## Lancer le pipeline

Le pipeline opère suivant 2 modes distincts selon le nombre de fastqs d'ADN tumoral circulant spécifiés au niveau du config file  :

- [ Mode Tumeur versus Normal sans Viscap ](DAG/full_pipeline_consensus_2samples.pdf) (1 gDNA normal, 1 gDNA tumeur, au plus 1 cfDNA tumeur circulant):
  
    Dans ce cas le workflow suivant va etre ignoré :

      CNV_depth_coverage_10 -> viscap_cnv_calling_11 -> lookup_and_export_viscap_xls_to_bed_12 -> Annot_SV_CNV_viscap_annotation_13 
 
- [ Mode Tumeur versus Normal avec Viscap ](DAG/full_pipeline_consensus_5samples.pdf) (1 gDNA normal, 1gDNA tumeur, au moins 2 cfDNA tumeur circulant)




```

#!/bin/bash



cd $DIR
/svc/snakemake/miniconda/miniconda3/envs/snakemake3/bin/snakemake \-s $SNAKEFILE \
--configfile $CONFIG \
--cluster-config $CLUSTER_CONFIG \
--cluster 'godjob.py create -n {cluster.name} -t {cluster.tags} --external_image -v {cluster.volume_snakemake} -v {cluster.volume_home} -v {cluster.volume_workdir} -v {cluster.volume_annotations} -c {cluster.cpu} -r {cluster.mem} -i {cluster.image} -s' \-j 5 -w 400 -k -p \2>&1 | tee -a logs/snake_full_pipeline_5p.log

```



## Script pour fusionner les datasets des variants annotés 

Le script merge_annotation_piepline.py fusionne les annotations SNP/indel, CNV/SV et microsatellites d'un pipeline
somatique (ANNOVAR, snpEff/SnpSift, AnnotSV, VisCap, MSIsensor...) en un seul tableau final, à partir de fichiers passés en argument.

Chaque type de fichier est fourni via une option répétable, avec l'outil de
variant calling et l'outil d'annotation donnés explicitement en argument :
    --vcf SAMPLE PATTERN CALLER        VCF annoté (ANNOVAR, snpEff/SnpSift...)
    --manta-vcf SAMPLE PATTERN CALLER  VCF brut de sortie Manta (SV non annotées)
    --cnv SAMPLE PATTERN CALLER        Tableau CNV/SV annoté (AnnotSV, VisCap, Manta...)
    --microsat SAMPLE PATTERN CALLER   Sortie MSIsensor (microsatellites)


PATTERN peut être un chemin exact ou un motif glob (ex: "*.vcf") : tous les fichiers correspondants sont lus et empilés pour le sample donné.

Exemple  pour l'échantillon IO890P :

```

python3 merge_annotation_pipeline.py \
--cnv  IO890P Annot_SV_CNV_annotation_21/UISTZ998_vs_IO890P.cns.vcf.annotated.tsv cnvkit \
--cnv  IO890P Annot_SV_CNV_viscap_annotation_13/IO890P.viscapCNV.annotated.tsv viscap  \
--vcf IO890P snp_sift_annotation_18/UISTZ998_vs_IO890P.ann_snpEff_snpSift.vcf mutect2  \
--vcf IO890P Annovar_Snp_and_indel_annotation_16/UISTZ998_vs_IO890P.hg38_multianno.vcf mutect2  \
--manta_vcf IO890P SV_workflow_14/results/variants/somaticSV.vcf manta \
--microsat IO890P  Msisensor_score_08/T_489UIIIOP_VS_N_IO890P.somatic.prefix_somatic msisensor \
--outdir files/final_dataset_merge.tsv

```




## License
For open source projects, say how it is licensed.
