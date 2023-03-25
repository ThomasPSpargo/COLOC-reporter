#!/bin/bash

######
### inputs
######

## Required setup

# Set directory containing all COLOC-reporter scripts
scriptpath=/scratch/users/k1802739/COLOC-reporter/scripts




#Directory within which to search. Assumes the working directory unless passed as a trailing argument to the script.
dir=${1:-$(pwd)}

#Declare the Rds file in which to return the summary Rds file. '.Rds' will be automatically appended.
#Defaults to the working directory with the prefix 'allSummaryFigures'
out=${2:-$(pwd)/allSummaryFigures}

#Options:

#Metrics to extract. For details, see: Rscript ./scripts/gwasSummaryPlotter.R --help
measures="p,PIP,z,beta"

#Indicate TRUE to mark in colour any finemapping credible sets or FALSE otherwise
includeFinemapping=TRUE

#Indicate TRUE to include 'full' trait names on the plots, otherwise, IDs will be used
plotTraitLabels=FALSE




######
### Run
######

#Normalise input and output filepaths
inputDir=$(readlink -f ${dir})
outputDir=$(readlink -f ${out})

#Run in temporary directory
tmp_dir=$(mktemp -d -t ci-XXXXXXXXXX)
cd  $tmp_dir

#Identify all harmonised sumstats files to concatenate across
ls -1 ${inputDir}/**/data/datasets/harmonised_sumstats.csv > collectfiles.txt

#Optionally set parameter which controls whether traits have full labels or their assigned IDs. If TRUE, extract trait labels
[[ $plotTraitLabels == TRUE ]] && params+=(--GWASconfig "${scriptpath}/GWAS_samples.txt")

#Loop across all files identified
loops=$(cat collectfiles.txt | wc -l)
for i in $(seq 1 $loops); do

	#Extract the summary statistics file for the loop
	file=$(awk  -F "\t" 'NR=='$i' {print $1}' collectfiles.txt)

	#Identify the column containing traits information
	traitColumn=$(awk 'NR==1{print $0}' $file | tr ',' '\n' | cat -n | grep trait | cut -f1)
	
	#Identify unique traits in the data and pivot into comma separated list
	traits=$(awk -F ',' 'NR>1{print $0}' $file | cut -d"," -f$traitColumn | sort | uniq | xargs | tr ' ' ',')

	printf "Generating plots for: ${file}\nTraits detected in dataset are: ${traits}\n"
	
	#This script takes harmonised summary statistics as input and plots comparisons across the measures indicated in --GWASsumplots across pairs of GWAS summary statistics
	Rscript ${scriptpath}/gwasSummaryPlotter.R \
			--traits $traits \
			--harmonisedSumstats ${file} \
			--GWASsumplots $measures \
			--GWASsumplots_incfinemapping $includeFinemapping \
			--helperFunsDir "${scriptpath}/helper_functions" \
			--rdsOut "./plots_$i" \
			--rdsOnly TRUE "${params[@]}"
done

#If Rds files have been output, read-in and concatenate them into a single list which will be returned in the output directory
matchRds=(./*.Rds)
if [[ -f ${matchRds[0]} ]]; then

R -s -e 'cat("Concatenating across summary figure lists... ")
separateplots <- list.files(".",full.names=TRUE,pattern=".Rds")
collectedplots <- lapply(separateplots,readRDS)

#Extract genomic region spanned by each plot and assign trait comparison names
chr_bp<- gsub(".*\n\\\((.*)\\\)","\\\\1",sapply(collectedplots, function(x)x[[1]]$bpfigure$labels$x))
#Apply meaningful names to the plots to indicate their analysis
names(collectedplots) <- paste0(lapply(collectedplots,function(x)paste0(levels(x$p$bpfigure$data$trait),collapse="_")),"_",chr_bp)

#Save the file
saveRDS(collectedplots,file="allSummaryPlots.Rds")
cat("Done!\n")'
	
	#Ensure the output directory exists and copy across the summary file
	mkdir -p $(dirname $outputDir)
	cp allSummaryPlots.Rds ${outputDir}.Rds

else 
	printf	"No output Rds files found to return.\n"
fi

#Remove all files and then the temp directory
rm $tmp_dir/*
rmdir $tmp_dir
