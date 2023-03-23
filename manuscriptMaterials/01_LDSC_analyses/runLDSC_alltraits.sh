#!/bin/bash
#SBATCH --job-name=ldsc_neuro_rg
#SBATCH --mem-per-cpu=5G
#SBATCH --time=0-1:30
#SBATCH --ntasks=1
#SBATCH --output=/scratch/users/k1802739/LAVA/%x.log

### Pre-setup NOTES:
# This script should run within the ldsc conda environment
# Running as a batch job expects to return log files within an existing 'logs' directory: mkdir -p /scratch/users/k1802739/LAVA/ldsc/logs

######
### SETUP
######

#Declare:
#directory for LDSC software
LDSCdir=/users/k1802739/ldsc

#path to ldsc reference info, precomputed for the 1KG EUR sample, and hm3 snplist
#Note that removal of the trailing '/' in ldref will cause ldsc.py to fail
ldref=/scratch/users/k1802739/ldsc/reference/eur_w_ld_chr/
mungesnplist=/scratch/users/k1802739/ldsc/reference/w_hm3.snplist

#Declare the main LAVA directory which contains the majority of files
rootdir=/scratch/users/k1802739/LAVA

#Declare path to (symbolic links for) sumstats
statdir=${rootdir}/inputs

#Create directory for LDSC analysis and set as WD
mkdir -p ${rootdir}/ldsc
cd ${rootdir}/ldsc

#Create directories in which to return univ and bivar results
mkdir -p ./ldsc_h2
mkdir -p ./ldsc_rg

######
### Prepare SUMSTATS, run univariate LDSC, and then bivariate LDSC for traits pairwise
######

##==========##==========##==========##==========##
#Munge all sumstats

#Declare summary statistic identifying prefixes
#In the same order as the sumstats, declare population prevalence per trait
declare -a sumfiles=(ALS AD PD SZ FTD)
declare -a traitPREV=(0.0029 0.1 0.027 0.004 0.00134)

#ALS 1/350=0.0029; Midpoint of the 1/300-1/400 estimates for ALS. doi(s): 10.1007/s00415-006-0195-y, 10.1111/j.1468-1331.2009.02586.x
#AD 1/10=0.1; Fig 1B assuming age ~85, doi: 10.1016/j.jalz.2013.10.005
#PD 1/37=0.027; https://www.parkinsons.org.uk/sites/default/files/2018-01/Prevalence%20%20Incidence%20Report%20Latest_Public_2.pdf
#SZ 4/1000=0.004; doi: 10.1371/journal.pmed.0020141
#FTD 1/742=0.00134; doi: 10.1212/WNL.0000000000002638

for trait in ${!sumfiles[@]}; do

	#Preprocess each set of sumstats
	python ${LDSCdir}/munge_sumstats.py \
	--sumstats  ${statdir}/${sumfiles[$trait]}.sumstats.gz \
	--out ./${sumfiles[$trait]}_munge \
	--chunksize 500000 \
	--merge-alleles $mungesnplist \
	--snp SNP \
	--a1 A1 \
	--a2 A2 \
	--p P \
	--N-col N \
	--signed-sumstats BETA,0 
	
	#Run univariate LDSC for h2 estimates
	python ${LDSCdir}/ldsc.py \
	--h2 ./${sumfiles[$trait]}_munge.sumstats.gz \
	--ref-ld-chr ${ldref} \
	--w-ld-chr ${ldref} \
	--samp-prev 0.5 \
	--pop-prev ${traitPREV[$trait]} \
	--out ./ldsc_h2/${sumfiles[$trait]}_h2
	
done

##==========##==========##==========##==========##
#Genetic correlation for ALS dataset and other traits:

declare -a nextTrait=(AD PD SZ FTD)
for next in ${nextTrait[@]}; do
	python ${LDSCdir}/ldsc.py \
	--rg ./ALS_munge.sumstats.gz,./${next}_munge.sumstats.gz \
	--ref-ld-chr ${ldref} \
	--w-ld-chr ${ldref} \
	--out ./ldsc_rg/ALS_${next}_rg
done

##==========##==========##==========##==========##
#Genetic correlation for AD:
declare -a nextTrait=(PD SZ FTD)
for next in ${nextTrait[@]}; do
	python ${LDSCdir}/ldsc.py \
	--rg ./AD_munge.sumstats.gz,./${next}_munge.sumstats.gz \
	--ref-ld-chr ${ldref} \
	--w-ld-chr ${ldref} \
	--out ./ldsc_rg/AD_${next}_rg
done

##==========##==========##==========##==========##
#Genetic correlation for PD dataset
declare -a nextTrait=(SZ FTD)
for next in ${nextTrait[@]}; do
	#vs SZ data
	python ${LDSCdir}/ldsc.py \
	--rg ./PD_munge.sumstats.gz,./${next}_munge.sumstats.gz \
	--ref-ld-chr ${ldref} \
	--w-ld-chr ${ldref} \
	--out ./ldsc_rg/PD_${next}_rg
done

##==========##==========##==========##==========##
#Final genetic correlation for SZ dataset vs FTD
python ${LDSCdir}/ldsc.py \
--rg ./SZ_munge.sumstats.gz,./FTD_munge.sumstats.gz \
--ref-ld-chr ${ldref} \
--w-ld-chr ${ldref} \
--out ./ldsc_rg/SZ_FTD_rg

##==========##==========##==========##==========##
