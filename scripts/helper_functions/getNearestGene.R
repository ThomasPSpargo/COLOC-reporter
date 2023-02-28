#####
# Author: Thomas Spargo (thomas.spargo@kcl.ac.uk)
# Obtained from GitHub repository: https://github.com/ThomasPSpargo/COLOC-reporter
#
# This script defines a function to identify the nearest genes upstream and downstream from a position and their basepair distance from that position
#####

## Inputs
#ID is an identifier for a single SNP
#positions is either a single numeric or an index of chromosomal positions and corresponding IDs, across the columns 'snp' and 'pos'
#genesInChr is a list of all genes to compare against in the chromosome

#Returns a list with the elements "upstream" and "downstream" which respectively are a single character string of the format: "ID (<nearest gene>; ±<sdistance>bp)"
getNearestGene<- function(ID,positions,genesInChr){
  if(!is.numeric(positions)){ varPos<- positions$pos[positions$snp==ID]
  } else { varPos  <- positions[1] }
  
  #Nested pasting to allow for possibility of multiple tags
  
  isWithin<- varPos>genesInChr$start_position & varPos<genesInChr$end_position
  if(any(isWithin)){ #If the snp occurs WITHIN a gene, return this gene as both upstream and downstream result, with 0 bp distance
    
    downstreamGene <- upstreamGene <- paste0(ID, " (",paste(genesInChr$external_gene_name[isWithin],collapse=", "),"; 0bp)")
    
  } else { #Otherwise, identify the nearest upstream and downstream gene relative to the position of the SNP
    genesUpstream <- genesInChr[varPos<genesInChr$start_position,]
    distancesUpstream <- genesUpstream$start_position-varPos
    upstreamGene <- paste0(ID, " (",paste(genesUpstream$external_gene_name[which(distancesUpstream==min(distancesUpstream))],collapse=", "),"; +", min(distancesUpstream),"bp)")
    
    genesDownstream <- genesInChr[varPos>genesInChr$end_position,]
    distancesDownstream <- varPos-genesDownstream$end_position
    downstreamGene <- paste0(ID, " (",paste(genesDownstream$external_gene_name[which(distancesDownstream==min(distancesDownstream))],collapse=", "),"; -", min(distancesDownstream),"bp)")
    
  }
  return(list(upstream=upstreamGene,downstream=downstreamGene))
}
