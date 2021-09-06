#!/bin/bash

## DEFINE PARAMETERS + PATHS + FILES ##
CHRNO=$1
CHROM="chr${CHRNO}"

REG_ID="PCDx92" #ID or name of gene group
GROUP="NCFB_probands" #ID or name of (vcf file) cohort
NUMBERF="123" #Number of vcf files (participants in cohort)
FREQ="MAF-noFilt"

SAMPLE_PATH="/file/path/to/analysis/directory/${GROUP}/perCHROM_norm_${REG_ID}"

OUTFILE3="${SAMPLE_PATH}/InHouse_MAFs_${GROUP}_norm-vep_${REG_ID}_x${NUMBERF}_UTRanno_${CHROM}.vcf.gz_Variants.txt"
PREpVARIANTS="${SAMPLE_PATH}/preVariants_${GROUP}_${FREQ}_${REG_ID}_${CHROM}.txt_allCILgeneVars-extended"

## extract intial set of variant info and counts to combine case / control genotype counts using R-script
module load lang/Perl/5.30.0-GCCcore-8.3.0
perl -lane 'print "$F[0]\t$F[1]\t$F[3]\t$F[4]\t$F[13]\t$F[7]\t$F[8]\t$F[9]\t$F[11]\t$F[20]\t$F[21]\t$F[22]\t$F[30]\t$F[31]\t$F[33]\t$F[34]\t$F[35]\t$F[41]"' "${OUTFILE3}" > "${PREpVARIANTS}"

##