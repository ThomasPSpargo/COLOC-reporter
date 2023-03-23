#!/bin/bash

#### OBTAIN FILES
#wget calls commented out to prevent accidental rerunning

#${1} Declares the directory to work within

SUMSTAT=${1}/AD_sumstats_Kunkle.txt
cd ${1}

#Stage 1 sumstats
wget -q -O Kunkle_etal_Stage1_results.txt https://www.niagads.org/system/tdf/public_docs/Kunkle_etal_Stage1_results.txt?file=1
  
#Step 2 sumstats - not run
#wget -O ~/SUMSTATS/AD/AD_stg2_sumstats_Kunkle.txt https://www.niagads.org/system/tdf/public_docs/Kunkle_etal_Stage2_results.txt?file=1

#wget -O SUMSTATS/AD/Kunkle_readme.docx .https://www.niagads.org/system/tdf/public_docs/Kunkle_etal_2019_IGAP_summary_statistics_README_0.docx?file=1

#Create new version of the file 
cp Kunkle_etal_Stage1_results.txt $SUMSTAT

#Extract original header from file
cat $SUMSTAT | awk 'NR==1{print $0}' > ${SUMSTAT}.header

#Replace header with LAVA-friendly colnames

#Original headers
#Chromosome Position MarkerName Effect_allele Non_Effect_allele Beta SE Pvalue
#New headers
cols="CHR BP SNP A1 A2 BETA SE P"

#Substitute entirety of first row with new colnames
sed -i "1s/.*/$cols/" $SUMSTAT

#Effective N


#Effective N calculated per-stratum and summed, based in per-phenotype Ns in supplement
#Ncases=21982
#Ncontrols=41944
#totn=$((Ncases+Ncontrols))
#N=$(python -c "print( 4*$totn*($Ncases / float($totn))*$Ncontrols / float($totn) )")

### ADGC
Ncases=14428
Ncontrols=14562
totn=$((Ncases+Ncontrols))
ADGC_N=$(python -c "print( 4*$totn*($Ncases / float($totn))*$Ncontrols / float($totn) )")

### CHARGE
Ncases=2137
Ncontrols=13474
totn=$((Ncases+Ncontrols))
CHARGE_N=$(python -c "print( 4*$totn*($Ncases / float($totn))*$Ncontrols / float($totn) )")

### EADI
Ncases=2240
Ncontrols=6631
totn=$((Ncases+Ncontrols))
EADI_N=$(python -c "print( 4*$totn*($Ncases / float($totn))*$Ncontrols / float($totn) )")

### GERAD 
Ncases=3177
Ncontrols=7277
totn=$((Ncases+Ncontrols))
GERAD_N=$(python -c "print( 4*$totn*($Ncases / float($totn))*$Ncontrols / float($totn) )")

#Sum for the meta-analysis effective N
effN=$(python -c "print( float($ADGC_N)+float($CHARGE_N)+float($EADI_N)+float($GERAD_N) )")

echo "EFFN is $effN"

#Add a column detailing total N, rounded to 2dp
awk -v num=$effN '{s=(NR==1)?"N":sprintf("%.2f",num);$0=$0 OFS  s}1' $SUMSTAT > temp.tmp && mv temp.tmp $SUMSTAT

######
### PROCESS THE SNP ID column to retain rsIDs and remove IDs which would not match to the reference
######

#0 records remain after filtering on these params
#awk 'BEGIN{FS = OFS = " "} NR>1 && !/^[0-9]+ [0-9]+ (rs|chr).*$/{print $3}' ~/SUMSTATS/AD/AD_sumstats_Kunkle.txt | wc -l

#awk 'BEGIN{FS = OFS = " "} NR>1 && /^[0-9]+ [0-9]+ rs.*$/{print $3}' ~/SUMSTATS/AD/AD_sumstats_Kunkle.txt | wc -l
#10534426 of 11480632 variants are rsIDs

#awk 'BEGIN{FS = OFS = " "} NR>1 && /^[0-9]+ [0-9]+ chr.*$/{print $3}' ~/SUMSTATS/AD/AD_sumstats_Kunkle.txt | wc -l
#940980 are chr:pos...

#awk 'BEGIN{FS = OFS = " "} NR>1 && /^[0-9]+ [0-9]+ merged_del.*$/{print $3}' ~/SUMSTATS/AD/AD_sumstats_Kunkle.txt | wc -l
#2939 are merged_del

#awk 'BEGIN{FS = OFS = " "} NR>1 && /^[0-9]+ [0-9]+ NA .*$/{print $3}' ~/SUMSTATS/AD/AD_sumstats_Kunkle.txt | wc -l
#2287 are NA

### REMOVE ALL NON RSID identifiers
### gsub the first column where the marker is not the first row and rsID is not matched
awk 'BEGIN{FS = OFS = " "} NR>1 && !/^[0-9]+ [0-9]+ rs.*$/{gsub(".*","NA",$3)} {print}' $SUMSTAT > temp.tmp && mv temp.tmp $SUMSTAT

