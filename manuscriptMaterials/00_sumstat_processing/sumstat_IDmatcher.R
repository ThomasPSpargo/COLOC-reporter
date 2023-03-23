#!/usr/bin/Rscript
# This script was written by Thomas Spargo (thomas.spargo@kcl.ac.uk),
# The intent for this script is a GWAS summary statistic preprocessing step before passing to the sumstat_cleaner script from GenoPredPipe (https://github.com/opain/GenoPred).
# Please see GenoPredPipe for further details regarding summary statistic formatting.

#GWAS INPUT FORMAT
### Expect a gwas summary statistic dataset with a minimum of the columns:
# CHR: chromosome 
# BP: position
# A1: effect allele
# A2: reference allele
# BETA: snp effect estimate

#Optionally, include the columns:
#SNP, wherein the script will attempt to match any missing SNP ids with the reference panel
#FREQ, where, if --isMAF is set to FALSE, 1-FREQ will be returned for any alleles with reversed allele order


start.time <- Sys.time()
library(optparse)
library(data.table)
library(dplyr)

option_list = list(
  make_option("--sumstats", action="store", default=NA, type='character',
              help="Path to summary statistics file [required]"),
  make_option("--ref_plink_chr", action="store", default=NA, type='character',
              help="Path to  PLINK binary reference files. If provided per-chromosome, prefix before chromosome number must end in 'chr' [required]"),
  make_option("--isMAF", action="store", default=F, type='logical',
              help="If FREQ column is provided and refers to minor allele frequency, set true to suppress FREQ column reversal if allele order is reversed."),
  make_option("--writeConflicts", action="store", default=T, type='logical',
              help="If SNP column is present and id matching is performed to fill in variants without IDs, set TRUE to return a file detailing ID conflicts between the GWAS and reference."),
  make_option("--gz", action="store", default=T, type='logical',
              help="Set to T to gzip summary statistics [optional]"),
  make_option("--output", action="store", default='./Output', type='character',
              help="Path and prefix for output files [optional]")
)

opt = parse_args(OptionParser(option_list=option_list))

opt$output_dir<-paste0(dirname(opt$output),'/')
system(paste0('mkdir -p ',opt$output_dir))

sink(file = paste(opt$output,'.log',sep=''), append = F)
cat(
  '#################################################################
# sumstat_IDmatcher.R
# This script was written by Thomas Spargo (thomas.spargo@kcl.ac.uk) please reach out with any queries.
#################################################################\n')
cat('Options are:\n')
print(opt)
cat('Analysis started at',as.character(start.time),'\n')
sink()

#####
# Read in sumstats
#####

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Reading in GWAS sumstats.\n')
sink()

if(substr(opt$sumstats,(nchar(opt$sumstats)+1)-3,nchar(opt$sumstats)) == '.gz'){
  GWAS<-data.frame(fread(cmd=paste0('zcat ',opt$sumstats)))
} else {	
  GWAS<-data.frame(fread(opt$sumstats))
}

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('GWAS contains',dim(GWAS)[1],'variants.\n')
sink()

#####
# Read in ref_bim
#####

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Reading in ref_plink_chr\n')
sink()

#If the string 'chr' is detected, then expect per-chromosome reference files, otherwise, read from a single file
if(grepl("chr$",opt$ref_plink_chr)){
  ref_bim<-NULL
  for(i in 1:22){
    bim<-fread(paste0(opt$ref_plink_chr,i,'.bim'))
    ref_bim<-rbind(ref_bim, bim)
  }
} else {
  ref_bim<-fread(opt$ref_plink_chr)
}

names(ref_bim)<-c('CHR','SNP','POS','BP','A1','A2')
ref_bim$POS<-NULL

#Ensure snps are all upper case
ref_bim$A1 <- toupper(ref_bim$A1)
ref_bim$A2 <- toupper(ref_bim$A2)

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('ref_plink_chr contains',dim(ref_bim)[1],'variants.\n')
sink()

#####
# Match snps to RSIDs from the reference panel based on chromosomal position 
#####

if(all(c("CHR", "BP") %in% colnames(GWAS)) #&& 
   #(!"SNP" %in% colnames(GWAS) || any(is.na(GWAS$SNP)))
){
  sink(file = paste(opt$output,'.log',sep=''), append = T)
  if(!"SNP" %in% colnames(GWAS)){ cat('SNP column not detected in file, ')
  } else if(any(is.na(GWAS$SNP))){cat('Missing values were detected in the SNP column, ')
  }
  cat('passing GWAS sumstats and reference to bigsnpr snp_match function to assign SNP ids based on CHR and POS columns.\n')
  cat('Please note that that this step is necessarily predicated on the assumption of corresponding alignment between the summary statistics and reference.\n')
  cat("Aligning snp column...")
  sink()
  
  ### TS NOTES: if the stat is OR, then bigsnpr may not reverse.... consider transforming OR into beta before passing through bigsnpr
  library(bigsnpr)
  library(dplyr)
  
  #Store the old colnames
  oldcols <- colnames(GWAS)
  
  #Assign bigsnpr-friendly column names to the reference. a1 = alt allele (allele1), a0 = ref allele (allele2)
  names(ref_bim)<-c("chr","id","pos","a1","a0")
  
  #Rename SUMSTAT cols to correspond with bigsnpr, match IDs, revert colnames, and return GWAS as before 
  #return_flip_and_rev = TRUE is set in order to ascertain when to reverse allele frequencies, if they have been reversed
  GWAS<- GWAS %>%
    rename(chr=CHR,
           pos=BP,
           a1=A1,
           a0=A2,
           beta=BETA,
           p=P) %>%
    mutate(a1=toupper(a1), #Ensure characters are upper case
           a0=toupper(a0)) %>%
    bigsnpr::snp_match(.,ref_bim,remove_dups = FALSE,return_flip_and_rev = TRUE)#,match.min.prop = 0) 
  
  if("SNP" %in% colnames(GWAS)){
    
    #Compare ids which are provided as a sanity check
    bothmatch<- which(!is.na(GWAS$SNP) & !is.na(GWAS$id))
    mismatch<- GWAS$SNP[bothmatch]!=GWAS$id[bothmatch]
    
    sink(file = paste(opt$output,'.log',sep=''), append = T)
    cat("Done!\n")
    cat('IDs were found for',sum(is.na(GWAS$SNP) & !is.na(GWAS$id)),'SNPs which were with previously missing \n')
    cat('SNP ID conflicts occurred for',sum(mismatch),"of", dim(GWAS[!is.na(GWAS$SNP),])[1],'SNPs for which rsIDs were originally provided the GWAS position matching with reference. The reference panel ID has been assigned in these circumstances.\n')
    sink()
    
    if(sum(mismatch)>0 && opt$writeConflicts){
      #Create a file for the SNP id conflicts, based on the output name, removing any file extension provided
      conflictfile <- file.path(dirname(opt$output),paste0(gsub('\\..*','',basename(opt$output)),"_IDconflicts.tsv"))
      
      sink(file = paste(opt$output,'.log',sep=''), append = T)
      cat("Details for variants with conflicting IDs have been written to:", conflictfile,"\n")
      sink()
      
      #Add file path as option?
      GWAS[bothmatch[mismatch],] %>%
        select(SNP,id,chr,pos,a1,a0) %>%
        rename(GWAS_ID=SNP,
               referenceID=id,
               CHR=chr,
               pos=pos,
               A1=a1,
               A2=a0) %>%
        fwrite(.,conflictfile, sep='\t')
    } else if(sum(mismatch)>0 && !opt$writeConflicts){
      sink(file = paste(opt$output,'.log',sep=''), append = T)
      cat("Set the --writeConflicts option to TRUE to return a summary of conflicting SNPs.\n")
      sink()
    }
    
    #Harmonise the SNP ID columns, prioritising the original IDs
    #GWAS <- GWAS %>%
      #mutate(SNP=if_else(!is.na(SNP),SNP,id))
    GWAS <- GWAS %>%
      mutate(SNP=id)
  } else {
    sink(file = paste(opt$output,'.log',sep=''), append = T)
    cat("Done!\n")
    cat('IDs were found for,',sum(!is.na(GWAS$id)),'SNPs\n')
    sink()
    
    # #Rename the column
    GWAS <- GWAS %>%
      rename(SNP=id)
  }

  #Optionally reverse the allele frequency for variants where allele order has been flipped
  if(any(GWAS[["_REV_"]])){
    sink(file = paste(opt$output,'.log',sep=''), append = T)
    cat("Allele order was flipped for", sum(GWAS[["_REV_"]]),"variants. ")
    sink()
    
    if("FREQ" %in% colnames(GWAS)){
      GWAS[GWAS[["_REV_"]],] <- GWAS[GWAS[["_REV_"]],] %>%
        mutate(FREQ=1-FREQ)
      
      sink(file = paste(opt$output,'.log',sep=''), append = T)
      if(!opt$isMAF){
        cat("The FREQ column has been flipped to match, ")
      } else { 
        cat("The FREQ column has not been flipped to match, ")
      }
      cat("in accordance with the --isMAF option.\n")
      sink()
    }
  }
  
  #Revert the column names and restrict to cols present in the original sumstats file
  GWAS <- GWAS %>%
    rename(CHR=chr,
           BP=pos,
           A1=a1,
           A2=a0,
           BETA=beta,
           P=p) %>%
    dplyr::select(c("SNP",all_of(oldcols)))
  
  sink(file = paste(opt$output,'.log',sep=''), append = T)
  cat('After bigsnpr alignment of SNP IDs,',dim(GWAS)[1],'variants remain.\n')
  sink()
  
  #Revert ref_bim names to avoid any downstream issues
  names(ref_bim) <- c('CHR','SNP','BP','A1','A2')

  # Depreciated setting to skip matching when a complete record of SNP IDs are already present provided
  #} #else if("SNP" %in% colnames(GWAS) && all(!is.na(GWAS$SNP))){
  # #Print number dropped
  # sink(file = paste(opt$output,'.log',sep=''), append = T)
  # cat('SNP column matching was not performed since IDs were provided for for all SNPs.\n')
  # cat('Returning sumstats file unaltered.\n')
  # sink()
  
} else if(any(!c("CHR", "BP") %in% colnames(GWAS))){
  #Print number dropped
  sink(file = paste(opt$output,'.log',sep=''), append = T)
  cat('SNP positional matching was not performed since the columns "CHR" and "BP" were not identified in the GWAS summary statistics.\n')
  cat('Returning sumstats file unaltered.\n')
  sink()
}

#####
# Write out results
#####

if(file.exists(paste0(opt$output,'.gz'))){
  system(paste0('rm ',opt$output,'.gz'))
}
if(file.exists(paste0(opt$output))){
  system(paste0('rm ',opt$output))
}

fwrite(GWAS, opt$output, sep='\t')

if(opt$gz == T){
  # Compress the GWAS sumstats
  system(paste0('gzip ', opt$output))
}

end.time <- Sys.time()
time.taken <- end.time - start.time
sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Analysis finished at',as.character(end.time),'\n')
cat('Analysis duration was',as.character(round(time.taken,2)),attr(time.taken, 'units'),'\n')
sink()