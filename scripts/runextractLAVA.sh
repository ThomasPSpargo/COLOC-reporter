#!/bin/bash

#Input options [default values]:
#	--indir ["."]
#	--extension [".bivar"]
#	--outfile ["./set.regions.txt"]
#	--pthresh [NULL] 
#	--useFDR [FALSE]
#	--outFormat ["genomicRegion"]
#	--returnPvalue [FALSE]
#	--regionWindow [0]
  
#For R script documentation call:
#Rscript ./scripts/extractLAVA.R --help
              
			  

MAIN=/scratch/users/k1802739/COLOC-reporter #Set main directory to work within
			  
#Obtain significant RG pairs
Rscript ${MAIN}/scripts/extractLAVA.R \
	--indir ${MAIN}/LAVA/results/ \
	--outfile ${MAIN}/scripts/set.regions.txt \
	--outFormat genomicRegion



