#!/bin/bash

## DEFINE PARAMETERS + PATHS + FILES ##
REG_ID="PCDx92" #ID or name of gene group
GROUP="NCFB_probands" #ID or name of (vcf file) cohort
NUMBERF="123" #Number of vcf files (participants in cohort)
FREQ="MAF-noFilt"
SAMPLE_ID="${GROUP}_norm_${REG_ID}"

SAMPLE_PATH="/file/path/to/analysis/directory/${GROUP}/perCHROM_norm_${REG_ID}"
SCRIPT_PATH="/file/path/to/analysis/scripts"

VCFs="${SAMPLE_PATH}/../${GROUP}_${REG_ID}_norm_vcfs.list"
utrVARIANTSc="ALLunfiltered_Counts_perVariant_noMAF-filter_NonBronchRels_NCFB_BRR_combined_PCDx92_GRCh38_allCHROMS.txt"

##Print text output of variants and participant genotypes in genes of interest
module load lang/Perl/5.30.0-GCCcore-8.3.0
perl "${SCRIPT_PATH}/G2_Get_proband_relative_PASS_genotypes_from_norm_VCFs.pl" ${SAMPLE_PATH} ${VCFs} ${SAMPLE_ID} ${utrVARIANTSc} ${FREQ} ${REG_ID}

##
