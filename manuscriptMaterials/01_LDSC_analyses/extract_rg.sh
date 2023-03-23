#!/bin/bash

#Set the WD as the first trailing argument passed to this script
cd ${1}
pwd

#Final location for the bivar intercept summary
destination=${2}

#Path to R script to generate heatmap, saved in the directory containing the RG summaries
heatmapscript=${3}

#Script from LAVA vignettes
	#Directory to files has been added

FILES=($(ls ./*_rg.log))		  # assuming the output format is [phenotypes]_rg.log

echo $FILES

#N=$(echo ${#FILES[@]})  		# and that all combinations of penotypes have been analysed

for I in ${FILES[@]}; do
        PHEN=$(echo $I | sed 's/_rg\.log//')
        
        # subset log files to relevant output
		tail -n 5 $I | head -n 2 > $PHEN.rg
        
        # add to single data set
        if [[ $I == ${FILES[0]} ]]; then
        	cat $PHEN.rg > all.rg		# only including the header for the first phenotypes
        else
        	cat $PHEN.rg | sed '1d' >> all.rg
        fi
done

#Note the reshape2 dependency

###THEN RUN R 
R -e '
scor = read.table("all.rg",header=T)              # read in
scor = scor[,c("p1","p2","gcov_int")]             # retain key headers
scor$p1 = gsub("_munge.sumstats.gz","",scor$p1)   # assuming the munged files have format [phenotype]_munge.sumstats.gz
scor$p2 = gsub("_munge.sumstats.gz","",scor$p2)   # (adapt as necessary)

scor$p1 = gsub("\\./","",scor$p1)
scor$p2 = gsub("\\./","",scor$p2)


traits <- unique(c(scor$p1,scor$p2)) #Identify unique traits analysed

newscor<- data.frame(p1=traits,p2=traits,gcov_int=1) #Simulate 100% correlation within comparison of matching traits

scor <- rbind(scor,newscor) #Append to scor

scor<- reshape2::acast(scor,  p1 ~p2,value.var = "gcov_int") #use reshape2 cast function to move to wide format

#Some pairings may be inserted into the lower triangle rather than upper.
reversed <- which(!is.na(scor[lower.tri(scor)])) #Identify which of the lower matrix has been filled

if(length(reversed)>0){
 unfilled <- which(is.na(scor[upper.tri(scor)])) #Identify which of the upper matrix has not been filled
 to_move <- scor[lower.tri(scor)][reversed] #Extract values from lower matrix 
 scor[upper.tri(scor)][unfilled] <- to_move #Reinsert into values upper matrix 
}

#Then use custom function to mirror upper covariance matrix 
makeSymm <- function(m) {
 m[lower.tri(m)] <- t(m)[lower.tri(m)]
 return(m)
}

mat <- makeSymm(scor)

mat = round(cov2cor(mat),5)                       # standardise
write.table(mat, "sample.overlap.txt", quote=F)   # save
'

##Move the output file to the LAVA output directory
cp sample.overlap.txt $destination


#Produce a heatmap for the total genetic correlation
Rscript $heatmapscript "./all.rg"