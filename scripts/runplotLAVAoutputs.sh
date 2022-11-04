#!/bin/bash
# This script plots univariate heritability and bivariate genetic correlations between pairs of phenotypes as returned from LAVA software, iterating across files ending in '.univ' and '.bivar' within the defined "$LAVAresult_dir" variable


MAIN=/scratch/users/k1802739/COLOC-reporter 			#Main directory to work within
LAVAresult_dir=${MAIN}/LAVA/results						#Directory containing LAVA output files
PlotOutput_dir=${MAIN}/LAVA/plots						#Directory in which to return the output plots


#For R script documentation call:
#Rscript ./scripts/plotLAVAoutputs.R --help



####
## Run
####

mkdir -p $PlotOutput_dir

for file in ${LAVAresult_dir}/*; do
	#Loop across files ending in .univ and .bivar from the directory
	if [[ $file == *.univ ]] || [[ $file == *.bivar ]]; then
		
		fname=$(basename $file) 							 #Drop the filepath to file
		echo ${fname} 										 #Echo filename
		P1P2=$(echo $fname | sed 's/\(\.bivar\|\.univ\)//g') #Drop the suffix
				
		#Run plotting script only when the files don't already exist
		if [[ $file == *.univ ]] && [[ -f ${PlotOutput_dir}/${fname}.manhat.pdf ]]; then
			echo "Univariate output figure for ${P1P2} already exists"
		
		elif [[ $file == *.bivar ]] && [[ -f ${PlotOutput_dir}/${fname}.manhat.pdf ]] && [[ -f ${PlotOutput_dir}/${P1P2}.rgcor.pdf ]]; then
			echo "Bivariate output figures for ${P1P2} already exist"
		else
		
			Rscript ${MAIN}/scripts/plotLAVAoutputs.R \
			--infile $file \
			--outdir ${PlotOutput_dir}/
		fi
	fi
done