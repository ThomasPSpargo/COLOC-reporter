#!/bin/bash

#### OBTAIN FILES
#wget calls commented out to prevent accidental rerunning

#${1} Declares the directory to work within

SUMSTAT=${1}/ALS_sumstats_vanRheenen.tsv
cd ${1}

wget -q http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90027001-GCST90028000/GCST90027164/GCST90027164_buildGRCh37.tsv.gz

zcat GCST90027164_buildGRCh37.tsv.gz > $SUMSTAT
#cp /scratch/users/k1802739/ldsc/ALS_2022/full_sumstat/european/GCST90027164_buildGRCh37.tsv $SUMSTAT

#Extract original header from file
cat $SUMSTAT | awk 'NR==1{print $0}' > ${SUMSTAT}.header

#Replace header with LAVA-friendly colnames

#Original headers
#rsid	chromosome	base_pair_location	effect_allele	other_allele	effect_allele_frequency	beta	standard_error	p_value	N_effective

#New headers
cols="SNP	CHR	BP	A1	A2	FREQ	BETA	SE	P	N"

#Substitute entirety of first row with new colnames
sed -i "1s/.*/$cols/" $SUMSTAT


######
### PROCESS THE SNP ID column to retain rsIDs and remove IDs which would not match to the reference
######

#0 records remain after filtering on these params
#awk 'BEGIN{FS = OFS = "\t"} NR>1 && !/^(rs|chr).*$/{print}' ~/SUMSTATS/ALS/ALS_sumstats_vanRheenen.tsv | wc -l

#awk 'BEGIN{FS = OFS = "\t"} NR>1 && /^rs.*$/{print}' ~/SUMSTATS/ALS/ALS_sumstats_vanRheenen.tsv | wc -l
#10439566 of 10461755 variants are rsIDs

#awk 'BEGIN{FS = OFS = "\t"} NR>1 && /^chr.*$/{print}' ~/SUMSTATS/ALS/ALS_sumstats_vanRheenen.tsv | wc -l
#22189 are chr:pos...

#Change the chrpos variants to 'NA'
### gsub the first column where the marker is not the first row and rsID is not matched
awk 'BEGIN{FS = OFS = "\t"} NR>1 && !/^rs.*$/{gsub(".*","NA",$1)} {print}' $SUMSTAT > temp.tmp && mv temp.tmp $SUMSTAT





