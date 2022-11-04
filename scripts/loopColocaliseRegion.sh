#!/bin/bash
# This script submits batch jobs based on input configurations and those of "${MAIN}/scripts/runColocaliseRegion.job". Batch jobs are submitted for each analysis indicated in the file directed to within the '$phenos' variable.



#For KCL CREATE HPC Platform: Load required modules - adapt as required
module load r
module load plink

#####
# Input settings
#####
MAIN=/scratch/users/k1802739/COLOC-reporter 		#Main directory to work within
LDREFERENCE=${MAIN}/ld_reference/EUR_phase3_chr		#Location and prefix for per-chromosome reference files in plink-binary format
GWASinfo=${MAIN}/scripts/GWAS_samples.txt			#Input file detailing configuration of traits and GWAS summary statistic inputs
plinkpath=plink										#If plink executable cannot be automatically identified by system, specify the relevant file path.

scriptpath=${MAIN}/scripts/colocaliseRegion.R		#Path identifying the rscript to run when submitting job
phenos=${MAIN}/scripts/coloc.phenopairs.txt			#Input file dictating phenotype pairings and loci to analyse


sumplots=PIP,p,beta 								#A comma-separated string passed to --GWASsumplots, see option help documentation
sumplots_onefile=FALSE								#Logical setting, set TRUE to provide all plots requested in $sumplots as a single figure


#####
# Run jobs
#####
mkdir -p ${MAIN}/logs 			#Make a directory in which to store logfiles
loops=$(cat $phenos | wc -l) 	#Determine how many times to loop to submit job for all unique phenotype pairs

for i in $(seq 1 $loops); do	#Loop for each row of $phenos
	
	# extract phenotype IDs and locus to analyse
	p1=$(awk 'NR=='$i' {print $1}' $phenos)
	p2=$(awk 'NR=='$i' {print $2}' $phenos)
	locus=$(awk 'NR=='$i' {print $3}' $phenos)
			
	#Set filepath and prefix for the analysis
	outpath=${MAIN}/coloc/results/${p1}.${p2}.${locus}

	#Print configuration of current analysis
	printf "########################\nTraits $p1 and $p2 at locus ${locus}\n"
	
	
	if [[ -f "${outpath}_coloc/colocalisation.log" ]]; then
		printf "The colocalisation.log file exists in the output directory, analysis is not run\n"
	else
		#If no output already, submit as a slurm job
		sbatch -J ${p1}.${p2}.${locus} -o ${MAIN}/logs/%x.log ${MAIN}/scripts/runColocaliseRegion.job $scriptpath $p1 $p2 $locus $GWASinfo $LDREFERENCE $outpath $sumplots $sumplots_onefile $plinkpath

		printf "\nOutput will be saved in directory with path and prefix:\n${outpath}\n"

	fi
	
	printf "____________\n"
	
		
done