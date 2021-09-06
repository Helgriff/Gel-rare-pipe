#!/bin/bash

## DEFINE PARAMETERS + PATHS + FILES ##
REG_ID="PCDx92" #ID or name of gene group
NUMBERF="123" #Number of vcf files (participants in cohort)
SAMPLE_ID="NCFB_probands" #ID or name of (vcf file) cohort

OUTPATH2="/file/path/to/analysis/directory/${SAMPLE_ID}"
OUTPATH3="${OUTPATH2}/perCHROM_norm_${REG_ID}"

if [ ! -d "${OUTPATH3}" ]; then
	mkdir "${OUTPATH3}"
fi

## Normalised (bcftools) vcf  file with no multi-allelic sites
vfile="${OUTPATH2}/InHouse_MAFs_${SAMPLE_ID}_norm_${REG_ID}_x${NUMBERF}.vcf"

####################################
module load bio/BCFtools/1.9-foss-2018b
module load bio/SAMtools/1.9-foss-2018b

tabix --list-chroms "${vfile}.gz" > "${OUTPATH2}/chromosomes.txt"  # save all the chromosome names into a file

while IFS= read -r line; do
	
chrvcf="${OUTPATH3}/InHouse_MAFs_${SAMPLE_ID}_norm_${REG_ID}_x${NUMBERF}_${line}.vcf"
	tabix "${vfile}.gz" $line > ${chrvcf};
	
	bgzip -f "${chrvcf}"
	tabix -f -p vcf "${chrvcf}.gz"
	
	done < "${OUTPATH2}/chromosomes.txt"  # make an individual vcf for each chromosome
