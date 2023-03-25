#!/bin/bash

######
### inputs
######

## Required setup

# Set directory containing all COLOC-reporter scripts
scriptpath=/scratch/users/k1802739/COLOC-reporter/scripts

# Set path to per-chromosome genotype reference files
LDREFERENCE=/scratch/users/k1802739/COLOC-reporter/ld_reference/EUR_phase3_chr

# #If plink executable cannot be automatically identified by system, specify the relevant file path
plinkpath=plink

#Directory containing all analyses. Assumes the working directory unless passed as a trailing argument to the script.
dir=${1:-$(pwd)}

#Declare the Rds file in which to return the summary Rds file. '.Rds' will be automatically appended.
#Defaults to the working directory with the prefix 'allSummaryFigures'
out=${2:-$(pwd)/allCSLDFigures}

######
### Run
######

#Normalise input and output filepaths
inputDir=$(readlink -f ${dir})
output=$(readlink -f ${out})

#Run in temporary directory
tmp_dir=$(mktemp -d -t ci-XXXXXXXXXX)
cd  $tmp_dir

#Loop across all analyses
for i in $(basename -a ${inputDir}/*/); do

	printf "Generating plots for the analysis: ${i}\n"
	
	#Generate cross-trait Credible-set LD heatmaps per analysis
	Rscript ${scriptpath}/crosstrait_CS_LDHeatmaps.R \
		--susieFitsDir ${inputDir}/${i}/data/finemapping \
		--LDreference $LDREFERENCE \
		--plink $plinkpath \
		--harmonisedSumstats ${inputDir}/${i}/data/datasets/harmonised_sumstats.csv \
		--helperFunsDir "${scriptpath}/helper_functions" \
		--rdsOut "./plots_$i" \
		--rdsOnly TRUE
	
done

#If Rds files have been output, read-in and concatenate them into a single list which will be returned in the output directory
matchRds=(./*.Rds)
if [[ -f ${matchRds[0]} ]]; then

R -s -e 'cat("Concatenating across per-analysis Rds files... ")
separateplots <- list.files(".",full.names=TRUE,pattern=".Rds")
collectedplots <- lapply(separateplots,readRDS)

#Assign names based on credible sets extracted and genomic region spanned
chr_bp<- gsub(".*\\\((.*)\\\)","\\\\1",sapply(collectedplots,function(x)x$heatmap_allSNPs$labels$x))
names(collectedplots) <- paste0("csFor_",lapply(collectedplots,function(x) paste0(unique(gsub(".*\\\((.*)\\\:.*","\\\\1",x$heatmap_allSNPs$data$leadSNP)),collapse="_")),"_",chr_bp)

saveRDS(collectedplots,file="allPlots.Rds")
cat("Done!\n")'
	
	#Ensure the output directory exists and copy across the summary file
	mkdir -p $(dirname $output)
	cp allPlots.Rds ${output}.Rds

else 
	printf	"No output Rds files found to return.\n"
fi

#Remove all files and then the temp directory
rm $tmp_dir/*
rmdir $tmp_dir
