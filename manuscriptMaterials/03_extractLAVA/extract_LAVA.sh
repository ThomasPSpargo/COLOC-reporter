#!/bin/bash
# This script extracts LAVA regions for colocalisation analysis and produces some summary figures of the regions in which local genetic correlation analysis was performed


####
## Set parameters for script
####

MAIN=/scratch/users/k1802739/LAVA 					#Main directory to work within
extractLAVAscripts=${MAIN}/scripts/03_extractLAVA	#Directory containing the scripts for extracting LAVA results
LAVAresult_dir=${MAIN}/results						#Directory containing LAVA output files
PlotOutput_dir=${MAIN}/plots						#Directory in which to return the output plots


####
## Determine number of bivariate comparisons performed and concatenate all bivariate outputs
####

#Initial number of comparisons
ncompars=0

echo "" #Echo blank line here and below for readability of command line output


declare -a LAVARes=($(echo ${LAVAresult_dir}/*.*.bivar)) #Identify all the lava bivariate results files
cat ${LAVARes[0]} | head -n 1 > ${MAIN}/all.bivar #Extract the header from the first file

#Ignoring file header, sum the total number of rows across the results file, indicating number of comparisons
for file in ${LAVARes[@]}; do
	
	#Add up across loops for total number of comparisons
	newcompars=$(awk 'NR>1{print $0}' $file | wc -l)
	ncompars=$((ncompars+newcompars)) 
	
	#Concatenate across files dropping header
	awk 'NR>1{print $0}' $file >> ${MAIN}/all.bivar
done

#Calculate the bonferroni P-value threshold and write to console
pthresh=$(echo "scale=10 ; 0.05 / $ncompars" | bc)
echo "N of bivariate comparisons: ${ncompars}, pvalue threshold set at 0.05/${ncompars}: ${pthresh}"

####
## Obtain COLOC-reporter set.regions.txt format from LAVA results
####			  

#For R script documentation call: Rscript extractLAVA.R --help

## Extract regions significant after strict bonferroni correction (this way of calling loops across individual scripts)
Rscript ${extractLAVAscripts}/getRGregions.R \
	--indir $LAVAresult_dir \
	--outfile ${MAIN}/set.regions_bonferroni.txt \
	--pthresh $pthresh \
	--outFormat genomicRegion \
	--returnPvalue TRUE
	
## Extract regions after correction for FDR (this points specifically to the all.bivar summary file and adjusts for all comparisons within that file)
#Write out a LAVA results summary file for the significant rows only with --fullRowsToFile option. Write out all rows with --everythingToFile option (this is useful for obtaining the FDR adjusted p-value)
Rscript ${extractLAVAscripts}/getRGregions.R \
	--indir ${MAIN}/ \
	--extension "all.bivar" \
	--outfile ${MAIN}/set.regions_globalFDR.txt \
	--useFDR TRUE \
	--outFormat genomicRegion \
	--returnPvalue TRUE \
	--fullRowsToFile ${MAIN}/LAVAsigrows.tsv \
	--everythingToFile ${MAIN}/LAVAall.tsv
	
####
## Plot LAVA results, looping across files
####
#These figures aim to give an overview of the overall LAVA results but are not necessarily ideal visualisations.

#For R script documentation call:
#Rscript plotLAVAoutputs.R --help

mkdir -p $PlotOutput_dir

for file in ${LAVAresult_dir}/*; do
	#Loop across files ending in .univ and .bivar from the directory
	if [[ $file == *.univ ]] || [[ $file == *.bivar ]]; then
		
		#Drop the filepath to file
		fname=$(basename $file) 	
		echo ${fname}
		
		#Drop the suffix
		P1P2=$(echo $fname | sed 's/\(\.bivar\|\.univ\)//g') 
				
		#Run plotting script only when the files don't already exist
		if [[ $file == *.univ ]] && [[ -f ${PlotOutput_dir}/${fname}.manhat.pdf ]]; then
			echo "Univariate output figure for ${P1P2} already exists"
		
		elif [[ $file == *.bivar ]] && [[ -f ${PlotOutput_dir}/${fname}.manhat.pdf ]] && [[ -f ${PlotOutput_dir}/${P1P2}.rgcor.pdf ]]; then
			echo "Bivariate output figures for ${P1P2} already exist"
		else
		
			Rscript ${extractLAVAscripts}/plotLAVAoutputs.R \
			--infile $file \
			--outdir ${PlotOutput_dir}/
		fi
	fi
done