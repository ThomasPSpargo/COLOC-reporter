#!/bin/bash

#### OBTAIN FILE
#gdown call commented out to prevent accidental rerunning

# ${1} Declares the directory to work within
SUMSTAT=${1}/PD_sumstats_Nalls.tsv
cd ${1}

#https://stackoverflow.com/questions/25010369/wget-curl-large-file-from-google-drive
gdown -O . 1FZ9UL99LAqyWnyNBxxlx6qOUlfAnublN

unzip -d . nallsEtAl2019_excluding23andMe_allVariants.tab.zip

mv nallsEtAl2019_excluding23andMe_allVariants.tab $SUMSTAT
#######

#Extract original header from file
cat $SUMSTAT | awk 'NR==1{print $0}' > ${SUMSTAT}.header

#Replace header with LAVA-friendly colnames

#Original headers
#SNP	A1	A2	freq	b	se	p	N_cases	N_controls

##Not run
##Identify effect allele based on labelling from manuscript table
##cat ~/SUMSTATS/PD/PD_sumstats_Nalls.tsv | grep -e '161469054' 

#New headers
cols="CHR:BP\tA1\tA2\tFREQ\tBETA\tSE\tP\tNCAS\tNCON"

#Substitute entirety of first row with new colnames
sed -i "1s/.*/$cols/" $SUMSTAT

#Compute effective N for each SNP based on the case control sample sizes
awk 'BEGIN{ FS = OFS = "\t" } {s=(NR==1)?"N":4*($8 + $9) * ($8 /($8 + $9)) * ($9 / ($8 + $9));$0=$0 OFS  s}1' ${SUMSTAT} > temp.tmp && mv temp.tmp $SUMSTAT

#Split the chr:BP column based on ':' separator; sed the 1st column to drop the 'chr' character prefix 
awk '{n=split($0,a,/[:]/); print a[1]"\t"a[2]}' ${SUMSTAT} | sed 's/chr//g' > ${SUMSTAT}.chrpos

#Overwrite original file to avoid redundancy 
mv ${SUMSTAT}.chrpos ${SUMSTAT}

######
### There is no SNP id column for these sumstats
######
