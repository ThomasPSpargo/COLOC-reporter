#!/bin/bash
#SBATCH --job-name=CleanCOLOCsumstats
#SBATCH --mem-per-cpu=20G
#SBATCH --time=0-2:00
#SBATCH --ntasks=8
#SBATCH --output=/users/k1802739/SUMSTATS/%x.log
## This header is written for a slurm scheduler

#Argument to supply indicating to rerun completed steps if TRUE
override=${1:-FALSE}

#Root directory to store all summary statistic files
rootdir="/users/k1802739/SUMSTATS"
cd $rootdir

###########

### NOTE: Presence of the following directory, containing required scripts, is assumed:
#${rootdir}/scripts
#The '00_summary_statistic_processing' directory is used in the GitHub repository for organisational purposes; this can be renamed as 'scripts'.

############

#Paths to reference plink and frequency files [not in the main directory]
plinkref=/scratch/users/k1802739/COLOC-reporter/ld_reference/EUR_phase3_chr
freqref=/scratch/users/k1802739/COLOC-reporter/ld_reference/freq/EUR_phase3_chr

#Run the preprocessing scripts for each of the sumstats, if their output file isnt already detected or override is TRUE
mkdir -p ${rootdir}/SZ/
if [[ ! -f ${rootdir}/SZ/SZ_sumstats_PGC.tsv || $override == TRUE ]]; then
	echo "---- Preprocessing SZ ----"
	bash ${rootdir}/scripts/prep_SZ2022_sumstats.sh ${rootdir}/SZ
fi

mkdir -p ${rootdir}/ALS/
if [[ ! -f ${rootdir}/ALS/ALS_sumstats_vanRheenen.tsv || $override == TRUE ]]; then
	echo "--- nPreprocessing ALS ----"
	bash ${rootdir}/scripts/prep_ALS_sumstats.sh ${rootdir}/ALS
fi

mkdir -p ${rootdir}/FTD/
if [[ ! -f ${rootdir}/FTD/FTD_GWAS_META.txt || $override == TRUE ]]; then
	echo "---- Preprocessing FTD ----"
	bash ${rootdir}/scripts/prep_FTD_sumstats.sh ${rootdir}/FTD
fi

mkdir -p ${rootdir}/PD/
if [[ ! -f ${rootdir}/PD/PD_sumstats_Nalls.tsv || $override == TRUE ]]; then
	echo "---- Preprocessing PD ----"
	bash ${rootdir}/scripts/prep_PD_sumstats.sh ${rootdir}/PD
fi

mkdir -p ${rootdir}/AD/
if [[ ! -f ${rootdir}/AD/AD_sumstats_Kunkle.txt || $override == TRUE ]]; then
	echo "---- Preprocessing AD ----"
	bash ${rootdir}/scripts/prep_AD_sumstats.sh ${rootdir}/AD
fi

#Space delimited list of GWAS sumstat files to clean in standardised protcol after initial preprocessing
declare -a sumstats=($rootdir/SZ/SZ_sumstats_PGC.tsv $rootdir/ALS/ALS_sumstats_vanRheenen.tsv $rootdir/AD/AD_sumstats_Kunkle.txt $rootdir/FTD/FTD_GWAS_META.txt $rootdir/PD/PD_sumstats_Nalls.tsv)

#Note that:
#The ID matching protocol recognises the columns:
#'SNP','CHR','BP','A1','A2','BETA''FREQ'

#The main sumstat_cleaner script recognises the columns:
#'SNP','A1','A2','P','OR','BETA','SE','N','FREQ','INFO'

######
### CHECK ID FORMATTING FOR whole-genome version of REFERENCE PANEL
### [NOTE: id formatting is checked per GWAS in their respective bash scripts]
######

#head /scratch/users/k1802739/COLOC-reporter/ld_reference/EUR_autosomes.bim

#awk 'BEGIN{FS = OFS = "\t"} /^[0-9]+\trs.*$/{print}' /scratch/users/k1802739/COLOC-reporter/ld_reference/EUR_autosomes.bim | wc -l
#14722395 of 14752705 records begin with rs

#awk 'BEGIN{FS = OFS = "\t"} /^[0-9]+\tesv.*$/{print}' /scratch/users/k1802739/COLOC-reporter/ld_reference/EUR_autosomes.bim | wc -l
#12139 records begin with esv

#awk 'BEGIN{FS = OFS = "\t"} /^[0-9]+\tss.*$/{print}' /scratch/users/k1802739/COLOC-reporter/ld_reference/EUR_autosomes.bim | wc -l
#18142 records begin with ss

#awk 'BEGIN{FS = OFS = "\t"} /^[0-9]+\t12_.*$/{print}' /scratch/users/k1802739/COLOC-reporter/ld_reference/EUR_autosomes.bim | wc -l
#28 records begin with 12_

#One record remains, beginning with 17_:
#awk 'BEGIN{FS = OFS = "\t"} !/^[0-9]+\t(rs|esv|ss|12_).*$/{print}' /scratch/users/k1802739/COLOC-reporter/ld_reference/EUR_autosomes.bim | head

#Summary:
#14722394 of 14752705 records begin with rs
#12139 records begin with esv
#18142 records begin with ss
#28 records begin with 12_
#One record remains, beginning with 17_

#THEREFORE, any CHR:POS ids in the sumstats must be converted to an rsID if possible. This is done using the sumstat_IDmatcher script. See each individual prep_*_sumstats file for a breakdown of how SNP IDs have been preprocessed.

#Base sumstat_cleaner script programmed by Oliver Pain for GenoPredPipe. I have made some tweaks to the script but the workflow is effectively the same.

#Rscript /users/k1802739/SUMSTATS/sumstat_cleaner_TS.R
#  --sumstats #Path to summary statistics file [required]  
#  --ref_plink_chr #Path to per chromosome PLINK files [required]
#  --ref_freq_chr #Path to per chromosome PLINK frequency (.frq) files [optional]
#  --info #INFO threshold [optional; default 0.9]
#  --maf #MAF threshold [optional; default 0.01]
#  --maf_diff #Difference between reference and reported MAF threshold [optional; default 0.2]
#  --insert_ref_maf #Set to T to insert reference allele frequency [optional; default T]
#  --gz #Set to T to gzip summary statistics [optional; default T]
#  --output #Path for output files [optional; default ./Output]
#  --SDcheck" #Specify 'quant' to perform check of whether SD is reasonable for all snps in a quantitative trait, or 'binary' to check for a binary trait. see doi: 10.1093/bioinformatics/btaa1029, section 3.4"),
 
#loop through the preprocessed summary statistic files
for sum in ${sumstats[@]}; do
		
	#Determine output directory, based on input name, dropping the file extension
	sout=$(echo $sum | sed 's/\..*//')

	echo "------------------------------------------"
	echo "Summary file is: $sum"
	
	if [[ ! -f ${sout}.idmatch.gz || $override == TRUE ]]; then 

		#Assign rsIDs that match with the reference; all datasets are GRCh37 aligned
		Rscript ${rootdir}/scripts/sumstat_IDmatcher.R \
			--sumstats $sum \
			--ref_plink_chr $plinkref \
			--writeConflicts TRUE \
			--output ${sout}.idmatch \
			--gz T

		echo "---- ID matching complete for $(basename $sum) ----"

	else
		echo "ID-matched file already exists for $(basename $sum). Jumping to main summary statistic cleaning"
	fi
	
	#Run main sumstat cleaning pipeline - retaining variants at MAF >0.005
	
	if [[ ! -f $sout.clean.pt05.gz || $override == TRUE ]]; then 
	
		#--SDcheck option has been implemented in the scripts but not performed in final run, reflecting unsuitability for GWAS meta-analysis
		Rscript ${rootdir}/scripts/sumstat_cleaner.R \
		  --sumstats ${sout}.idmatch.gz \
		  --ref_plink_chr $plinkref \
		  --ref_freq_chr $freqref \
		  --output ${sout}.clean.pt05 \
		  --info 0.9 \
		  --maf 0.005 \
		  --maf_diff 0.2 \
		  --insert_ref_maf T \
		  --gz T
		  
		  echo "---- Sumstat QC complete for $(basename $sum) ----"
	  
 	else
	  	echo "Outfile already exists for $(basename $sum), no action taken"
	fi
	
done