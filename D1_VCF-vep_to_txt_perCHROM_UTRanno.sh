#!/bin/bash

## DEFINE PARAMETERS + PATHS + FILES ##
CHRNO=$1
CHROM="chr${CHRNO}"

REG_ID="PCDx92" #ID or name of gene group
GROUP="NCFB_probands" #ID or name of (vcf file) cohort
NUMBERF="123" #Number of vcf files (participants in cohort)

SAMPLE_PATH="/file/path/to/analysis/directory/${GROUP}/perCHROM_norm_${REG_ID}"
SCRIPT_PATH="/file/path/to/scripts"

SYMBOLS="/file/path/to/analysis/directory/info_files/${REG_ID}_GeneSymbols.csv" #file listing gene symbols

OUTFILE3="InHouse_MAFs_${GROUP}_norm-vep_${REG_ID}_x${NUMBERF}_UTRanno_${CHROM}.vcf.gz"

##Run perl script
module load lang/Perl/5.30.0-GCCcore-8.3.0
perl "${SCRIPT_PATH}/D2_VCF-vep_to_txt_UTRanno.pl" ${SAMPLE_PATH} ${OUTFILE3} ${SYMBOLS} ${NUMBERF}
