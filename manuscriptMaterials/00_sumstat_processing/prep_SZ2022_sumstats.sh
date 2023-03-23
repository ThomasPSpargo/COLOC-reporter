#!/bin/bash

#${1} declares the directory in which to work

#Declare name of the sumstats file to be returned
SUMSTAT=${1}/SZ_sumstats_PGC.tsv
cd ${1}

#### OBTAIN FILES
#Source: https://figshare.com/articles/dataset/scz2022/19426775?file=34865091
	#PGC3_SCZ_wave3.european.autosome.public.v3.tsv.gz
wget -q -O PGC3_SCZ_wave3.european.autosome.public.v3.tsv.gz https://figshare.com/ndownloader/files/34517828

#extract the VCF header
zcat PGC3_SCZ_wave3.european.autosome.public.v3.tsv.gz | awk '{ if ($0 ~ /^#/) { print } }' > ${1}/extract.sumstat.header 

#extract the body (i.e. summary statistics to the main file)
zcat PGC3_SCZ_wave3.european.autosome.public.v3.tsv.gz | awk '{ if ($0 !~ /^#/) { print } }' > ${SUMSTAT}
 
#Extract summary statistic header from file
cat ${SUMSTAT} | awk 'NR==1{print $0}' > ${SUMSTAT}.header

#Replace headers
#Original: CHROM	ID	POS	A1	A2	FCAS	FCON	IMPINFO	BETA	SE	PVAL	NCAS	NCON	NEFF
#New headers
cols="CHR	SNP	BP	A1	A2	FCAS	FCON	INFO	BETA	SE	P	NCAS	NCON	N"

#Substitute entirety of first row with new colnames
sed -i "1s/.*/$cols/" $SUMSTAT

#Convert NEFF (which appears to be NEFF/2) into NEFF*2, rounding to 2dp with sprintf. Replace column 14, which is NEff
awk 'BEGIN{ OFS = "\t" } {s=(NR==1)? "N":sprintf("%.2f", $14*2);$14=s}1' $SUMSTAT > temp.tmp && mv temp.tmp $SUMSTAT


#Determine total allele frequency across cases and controls by weighted sum, store in new column FREQ,
cat $SUMSTAT | awk 'BEGIN{ OFS = "\t" } {s=(NR==1)? "FREQ":sprintf("%.5f",$6*($12/($12+$13))+$7*($13/($12+$13)));$0=$0 OFS  s}1' > temp.tmp && mv temp.tmp $SUMSTAT


######
### PROCESS THE SNP ID column to retain rsIDs and remove IDs which would not match to the reference
######
#0 records remain after filtering on these params
#awk 'BEGIN{FS = OFS = "\t"} NR>1 && !/^(rs|chr).*$/{print}' ~/SUMSTATS/SZ_2022/SZ_sumstats_PGC.tsv | wc -l

#awk 'BEGIN{FS = OFS = "\t"} NR>1 && /^[0-9]+\trs.*$/{print $2}' ~/SUMSTATS/SZ_2022/SZ_sumstats_PGC.tsv | wc -l
#7638748 of 7659767 variants are rsIDs

#awk 'BEGIN{FS = OFS = "\t"} NR>1 && /^[0-9]+\t[0-9]+:[0-9+].*$/{print $2}' ~/SUMSTATS/SZ_2022/SZ_sumstats_PGC.tsv | wc -l
#21019 are chr:pos not denoted with 'chr'  prefix

### REMOVE ALL NON RSID variants variants
### gsub the first column where the marker is not the first row and rsID is not matched
awk 'BEGIN{FS = OFS = "\t"} NR>1 && !/^[0-9]+\trs.*$/{gsub(".*","NA",$2)} {print}' $SUMSTAT > temp.tmp && mv temp.tmp $SUMSTAT

