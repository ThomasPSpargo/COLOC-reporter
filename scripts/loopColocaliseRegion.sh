#!/bin/bash
# This script submits batch jobs based on input configurations and those of "${MAIN}/scripts/runColocaliseRegion.job". Batch jobs are submitted for each analysis indicated in the file directed to within the '$phenos' variable.

#####
# Input settings
#####
MAIN=/scratch/users/k1802739/COLOC-reporter 		#Main directory to work within
LDREFERENCE=${MAIN}/ld_reference/EUR_phase3_chr		#Location and prefix for per-chromosome reference files in plink-binary format
plinkpath=plink										#If plink executable cannot be automatically identified by system, specify the relevant file path.

scriptpath=${MAIN}/scripts		#Path identifying the the location of scripts for running the job and the GWAS_samples.txt input file

comparisons=${MAIN}/scripts/set.regions.txt			#Input file dictating phenotype pairings and loci to analyse


sumplots=PIP,p,beta,z 								#A comma-separated string passed to --GWASsumplots, see option help documentation for the gwasSummaryPlotter.R script

runMode=doBoth		#Which analysis configuration to run (see documentation).


#####
# Run jobs
#####
mkdir -p ${MAIN}/logs 			#Make a directory in which to store logfiles
loops=$(sed 1d $comparisons | wc -l) 	#Determine how many times to loop to submit job for all unique phenotype pairs

for i in $(seq 1 $loops); do	#Loop for each row of $phenos
	
	# Extract phenotype IDs and locus to analyse.
	# The field separator is set explicitly to tab to ensure that the input file is correctly delimited.
	# i+1 is extracted since the is expected to have a header
	traits=$(awk  -F "\t" 'NR=='$(($i+1))' {print $1}' $comparisons)
	locus=$(awk  -F "\t" 'NR=='$(($i+1))' {print $2}' $comparisons)
	
	#Convert commas to a better delimiter for filepaths and the log
	analysis="$(echo $traits | sed -e 's/,/_/g').$(echo $locus | sed -e 's/,/_/g')"
	
	#Set filepath and prefix for the analysis
	outpath=${MAIN}/coloc/results/${analysis}

	#Print configuration of current analysis
	printf "########################\nTraits: ${traits} at locus: ${locus}\n"
	
	
	if [[ -f "${outpath}_coloc/colocalisation.log" ]]; then
		printf "The colocalisation.log file exists in the output directory, analysis is not run\n"
	else
		#If no output already, submit as a slurm job
		sbatch -J ${analysis} -o ${MAIN}/logs/%x.log ${MAIN}/scripts/runColocaliseRegion.job $scriptpath $traits $locus $LDREFERENCE $outpath $sumplots $plinkpath $runMode

		printf "\nOutput will be saved in directory with path and prefix:\n${outpath}\n"

	fi
	
	printf "____________\n"
	
		
done