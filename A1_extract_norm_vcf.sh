# !/bin/bash

## DEFINE FILES + PARAMETERS + PATHS ##
REGIONS="/file/path/Genes_coordinates.bed" #list of gene coordinates in bed format
VCFS="/file/path/vcf_files.list" #list of vcf file names including full path to file

REG_ID="PCDx92" #ID or name of gene group
NUMBERF="123" #Number of vcf files (participants in cohort)
GROUP="NCFB_probands" #ID or name of (vcf file) cohort
SAMPLE_ID="${GROUP}_norm_${REG_ID}"

SCRIPT_PATH="/file/path/to/scripts"

SAMPLE_PATH="/file/path/to/analysis/directory/${GROUP}"
if [ ! -d "${SAMPLE_PATH}" ]; then
	mkdir "${SAMPLE_PATH}"
fi

OUTPATH="${SAMPLE_PATH}/filtered_${REG_ID}"
if [ ! -d "${OUTPATH}" ]; then
	mkdir "${OUTPATH}"
fi

OUTPATH2="${SAMPLE_PATH}/norm_${REG_ID}"
if [ ! -d "${OUTPATH2}" ]; then
	mkdir "${OUTPATH2}"
fi

#######################################
module load bio/BEDTools/2.27.1-foss-2018b
module load bio/BCFtools/1.9-foss-2018b
module load bio/SAMtools/1.9-foss-2018b

## 1. EXTRACT  variants from vcf with bedtools and 2. NORMALISE vcf with bcftools (ie. splits multi-allelic sites) ##
while read -r vcf; do
	vfile="${vcf##*/}"
	OUTFILE="${OUTPATH}/${vfile}.${REG_ID}-region.vcf"
	bedtools intersect -a ${vcf} \
	-b ${REGIONS} \
	-u \
	-header > ${OUTFILE}

	OUTFILE2="${OUTPATH2}/${vfile}.${REG_ID}-norm.vcf.gz"	
	bcftools norm -o "${OUTFILE2}" -Oz -m-any "${OUTFILE}"
	tabix -f -p vcf "${OUTFILE2}"
done < ${VCFS}

## remove intermediate filtered vcfs
rm ${OUTPATH}/*.vcf

#####################################
module load lang/Perl/5.30.0-GCCcore-8.3.0

## 3. CREATE SINGLE 'vcf' of per cohort variants (includes HET/HOM genotype counts from total cohort)
VCFs="${SAMPLE_PATH}/${GROUP}_${REG_ID}_norm_vcfs.list"
if [ ! -f "${VCFs}" ]; then
	ls ${OUTPATH2}/*gz >${VCFs}
fi

perl "${SCRIPT_PATH}/A2_make_single_perCohort_VCF.pl" ${SAMPLE_PATH} ${VCFs} ${SAMPLE_ID}

#####################################
module load bio/tabix/0.2.6-GCCcore-7.3.0

## 4. GZIP and Index vcf file  ##############
vfile="${SAMPLE_PATH}/InHouse_MAFs_${GROUP}_norm_${REG_ID}_x${NUMBERF}.vcf"
bgzip -f "${vfile}"
tabix -f -p vcf "${vfile}.gz"

####################################
