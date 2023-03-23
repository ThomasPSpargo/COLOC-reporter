#!/bin/bash
#https://rdr.ucl.ac.uk/articles/dataset/IFGC_Summary-statistics_Data-sharing/13042166/1?file=24952892

#Declare the name of the file to work with, and directory to work within
SUMSTAT=${1}/FTD_GWAS_META.txt
cd $(dirname $SUMSTAT)


#### OBTAIN FILES
#wget calls commented out to prevent accidental rerunning

#Stage 1 sumstats
wget -q https://rdr.ucl.ac.uk/ndownloader/articles/13042166/versions/1

unzip 1
 #extracting: FTD_GWAS_bvFTD.txt      
 #extracting: FTD_GWAS_META.txt       
 #extracting: FTD_GWAS_MND.txt        
 #extracting: FTD_GWAS_PNFA.txt       
 #extracting: FTD_GWAS_SD.txt         
 #extracting: IFGC_Members-Affiliations-Acknowledgements.pdf  
 
#Display all invisible characters
#cat -A FTD_GWAS_META.txt 
 
#Remove carriage return from all unzipped sumstats [this is caused issues with AWK: (^M character)]
perl -p -i -e 's/\r\n$/\n/g' FTD_GWAS_*.txt

#Extract original header from file
cat ${SUMSTAT} | awk 'NR==1{print $0}' > ${SUMSTAT}.header

#Original headers
#marker Allele1	Allele2	beta1	SE	pValue	chr	Bp

#New headers, 'SNP' column contains a mix of rsIDs and chrpos records
cols="SNP	A1	A2	BETA	SE	P	CHR	BP"

#Substitute entirety of first row with new colnames
sed -i "1s/.*/$cols/" $SUMSTAT

#Formula: EffN=4*n_k*v_k*(1-v_k) v_k sampling proportion, n_k total samples

#Effective N calculated per-stratum and summed, based in per-phenotype Ns in paper/supplement
###Ncases=2154
###Ncontrols=4308
###totn=$((Ncases+Ncontrols))

### bvFTD
Ncases=1377
Ncontrols=2754
totn=$((Ncases+Ncontrols))
bvftd_N=$(python -c "print( 4*$totn*($Ncases / float($totn))*$Ncontrols / float($totn) )")

### Semantic dementia
Ncases=308
Ncontrols=616
totn=$((Ncases+Ncontrols))
SD_N=$(python -c "print( 4*$totn*($Ncases / float($totn))*$Ncontrols / float($totn) )")

### PNFA
Ncases=269
Ncontrols=538
totn=$((Ncases+Ncontrols))
PNFA_N=$(python -c "print( 4*$totn*($Ncases / float($totn))*$Ncontrols / float($totn) )")

### FTD - MND  
Ncases=200
Ncontrols=400
totn=$((Ncases+Ncontrols))
FTDMND_N=$(python -c "print( 4*$totn*($Ncases / float($totn))*$Ncontrols / float($totn) )")

#Sum for the meta-analysis effective N
#effN=$(($bvftd_N+$SD_N+$PNFA_N+$FTDMND_N))
effN=$(python -c "print( float($bvftd_N)+float($SD_N)+float($PNFA_N)+float($FTDMND_N) )")

echo "EFFN is $effN"

### Add a column detailing total N
awk -v num=$effN 'BEGIN{ FS = OFS = "\t" } {s=(NR==1)?"N":num;print $0 OFS  s}' $SUMSTAT > temp.tmp && mv temp.tmp $SUMSTAT

#####
### PROCESS THE SNP ID column to retain rsIDs and remove IDs which would not match to the reference
#####
#0 records remain after filtering on these params
#awk 'BEGIN{FS = OFS = "\t"} NR>1 && /^rs.*$/{print}'  ~/SUMSTATS/FTD/FTD_GWAS_META.txt | wc -l
#4812662 of 6026384 are rsid 

#awk 'BEGIN{FS = OFS = "\t"} NR>1 && /^chr.*$/{print}' ~/SUMSTATS/FTD/FTD_GWAS_META.txt | wc -l
#1213722 variants are chr

### gsub the first column where the marker is not the first row and rsID is not matched
awk 'BEGIN{FS = OFS = "\t"} NR>1 && !/^rs.*$/{gsub(".*","NA",$1)} {print}' $SUMSTAT > temp.tmp && mv temp.tmp $SUMSTAT

