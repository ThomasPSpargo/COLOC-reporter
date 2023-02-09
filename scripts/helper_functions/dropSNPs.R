#####
# Author: Thomas Spargo (thomas.spargo@kcl.ac.uk)
# Obtained from GitHub repository: https://github.com/ThomasPSpargo/COLOC-reporter
#
# This script defines a function for filtering down a coloc-format dataset including an LD matrix
# based on a logical index which identifies SNPs to remove based on the snp element
#####

#Input:
# x should be a coloc-dataset represented in a list and which includes the elements 'snp' and 'LD'. The LD matrix dims should correspond to IDs for SNPs in the 'snp' element.
# dropIndex is a logical vector containing TRUE and FALSE values where TRUE indicates positions to be removed from the dataset.

#Output:
# The coloc-dataset x after filtering out snps marked by dropIndex

dropSNPs<- function(x,dropIndex){
  dropIDs<- x$snp[dropIndex] #Extract IDs to match Dropped SNPs to an LD matrix
  
  #Identify which elements of the dset list are snp-wise vectors and retain only the elements of each of these which are not marked for exclusion
  SNPwise<- sapply(x,length)==length(x$snp)
  x[SNPwise] <- lapply(x[SNPwise],function(x)x[!dropIndex])
  
  #Using on the SNP IDs, drop corresponding elements from the LD matrix
  #dimnames(dset$LD) 1 and 2 should be identical, but matching performed independently by dim for certainty
  x$LD <- x$LD[!dimnames(x$LD)[[1]] %in% dropIDs,!dimnames(x$LD)[[2]] %in% dropIDs]        
  
  return(x)
}