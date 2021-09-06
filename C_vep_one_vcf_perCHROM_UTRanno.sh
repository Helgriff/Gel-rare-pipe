#!/bin/bash

## DEFINE PARAMETERS + PATHS + FILES ##
CHRNO=$1
CHROM="chr${CHRNO}"

REG_ID="PCDx92" #ID or name of gene group
NUMBERF="123" #Number of vcf files (participants in cohort)
SAMPLE_ID="NCFB_probands"  #ID or name of (vcf file) cohort

OUTPATH2="/file/path/to/analysis/directory/${SAMPLE_ID}"

## Annotate vcf with VEP
module load bio/VEP/99.1-foss-2019a-Perl-5.28.1
module load bio/tabix/0.2.6-GCCcore-7.3.0

	vfile="${OUTPATH2}/perCHROM_norm_${REG_ID}/InHouse_MAFs_${SAMPLE_ID}_norm_${REG_ID}_x${NUMBERF}_${CHROM}.vcf" ###per chromosome
	OUTFILE3="${OUTPATH2}/perCHROM_norm_${REG_ID}/InHouse_MAFs_${SAMPLE_ID}_norm-vep_${REG_ID}_x${NUMBERF}_UTRanno_${CHROM}.vcf.gz"
        
	vep --input_file "${vfile}.gz" \
	--output_file ${OUTFILE3} \
	--vcf \
	--compress_output bgzip \
	--offline \
	--cache \
	--dir_cache /path/to/vep/cache \
	--cache_version 99 \
	--species homo_sapiens \
	--assembly GRCh38 \
	--force_overwrite \
	--no_stats --fasta /path/to/ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
	--everything \
	--plugin UTRannotator,/path/to/uORF_starts_ends_GRCh38_PUBLIC.txt \
	--use_given_ref
	####################################
	tabix -f -p vcf "${OUTFILE3}"
        ####################################
