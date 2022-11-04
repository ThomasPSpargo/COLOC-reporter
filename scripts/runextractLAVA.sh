#!/bin/bash

#Input options [default values]:
#  --indir ["."]
#  --extension [".bivar"]
#  --outfile ["./coloc.phenopairs.txt"]
#  --pthresh [NULL] 
#  --outFormat ["genomicRegion"]
  
#For R script documentation call:
#Rscript ./scripts/extractLAVA.R --help
              
			  

MAIN=/scratch/users/k1802739/COLOC-reporter #Set main directory to work within
			  
#Obtain significant RG pairs
Rscript ${MAIN}/scripts/extractLAVA.R \
	--indir ${MAIN}/LAVA/results/ \
	--outfile ${MAIN}/scripts/coloc.phenopairs.txt \
	--outFormat genomicRegion



