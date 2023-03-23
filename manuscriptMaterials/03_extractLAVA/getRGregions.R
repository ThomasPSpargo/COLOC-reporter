#####
# Author: Thomas Spargo (thomas.spargo@kcl.ac.uk)
# Obtained from GitHub repository: https://github.com/ThomasPSpargo/COLOC-reporter
#
# The purpose of this script is to extract genomic regions highly correlated between two traits as produced by the LAVA software
# The main output is a plain-text file formatted for use with COLOC-reporter, labelling regions for fine-mapping and colocalisation analysis
#####


suppressMessages(library(optparse))
library(data.table) #Use data table for automatic detection of field separators 

option_list = list(
  make_option("--indir", action="store", default=".", type='character',
              help="The file path containing local genetic correlation files from lava. Files are extracted by selecting all of those within the directory with the extension set in --extension, which is by default '.bivar'"),
  make_option("--extension", action="store", default=".bivar", type='character',
              help="File extension used across the files to read in. Expects '.bivar' by default."),
  make_option("--outfile", action="store", default="./set.regions.txt", type='character',
              help="The name of the file to store the pheno pairs"),
  make_option("--pthresh", action="store", default=NULL, type='numeric',
              help="Manually assign the p-value threshold for loci to extract. If not set, the p-value threshold is 0.05/nrow(<file>). Option is ignored if useFDR is true."),
  make_option("--useFDR", action="store", default=FALSE, type='logical',
              help="Override the p-value threshold option and determine loci to extract with False-Discovery Rate p-value adjustment, adjusting for nrow(<file>)."),
  make_option("--outFormat", action="store", default="genomicRegion", type='character',
              help="Character string which may be any of: 'genomicRegion', 'LAVAlocus', 'traitsOnly'.\n'genomicRegion' is the default, and extracts regions below the p-value threshold in comma separated list of the format 'chromosome,start_position,end_position' (e.g. 17,43460501,44865832).\n'LAVAlocus' functions like genomicRegion, but instead extracts the locus number associated with the region.\n'traitsOnly' ignores any p-value thresholding, and instead extracts every pair of traits for which LAVA has an output - the output will only have 2 columns."),
  make_option("--returnPvalue", action="store", default=FALSE, type='logical',
              help="False by default. Set TRUE to include p-values for the regions extracted. Note that if --useFDR is true, then two columns will be returned. first p and the p following fdr adjustment."),
  make_option("--fullRowsToFile", action="store", default=NULL, type='character',
              help="Specify a file path to return an additional tab-separated output file including all rows below the threshold set."),
  make_option("--everythingToFile", action="store", default=NULL, type='character',
              help="Specify a file path to return an additional tab-separated output file including all rows. This is primarily for writing out the file with all accompanying FDR adjusted p-values."),
  make_option("--writeUnivariate", action="store", default=FALSE, type='logical',
              help="Logical, defaulting to FALSE. If TRUE, return the region list for traits in a univariate, rather than pairwise, manner. Automatically deduplicates repeated univariate analyses, and is not compatible with the --returnPvalue TRUE option; p-values are not returned."),
  make_option("--regionWindow", action="store", default=0, type='numeric',
              help="0 by default. Specify an additional number of base pairs to include around the region identified by LAVA")
)
opt = parse_args(OptionParser(option_list=option_list))

######
### Perform sanity checks that output setting is valid
######
if(!opt$outFormat %in% c('genomicRegion', 'LAVAlocus', 'traitsOnly')){
  stop("The --outFormat option has not been set correctly; please specify one of 'genomicRegion', 'LAVAlocus', or 'traitsOnly'. Refer to the help documentation for further details.")  
}

#Make the input directory the working directory
setwd(opt$indir)

######
### Identify files, if any are identified set output file and loop across them to identify regions to extract
######
files<-list.files(path=opt$indir,pattern=paste0(opt$extension,"$"),full.names = TRUE)

if(!length(files)>0){
  stop("No files with the extension '",opt$extension,"' were identified in the  directory:\n",getwd(),"\nPlease check documentation for the --indir and --extension options.")
}

for(i in 1:length(files)){
  lavabivar <- data.frame(fread(files[i]))
  
  #Set cols to extract in all settings of outFormat
  cols <- c("phen1","phen2")
  
  if(opt$outFormat=="traitsOnly"){
    rows <- 1 #Extract only row 1, to get each trait pairwise
  } else {
    
    ## Significance threshold options
    #If no p<sig (or p.fdr<0.05 when --useFDR is TRUE) then no rows will be identified, and no output will be written
    if(opt$useFDR){
      #Determine significant rows with FDR adjustment 
      lavabivar$p.fdr<- p.adjust(lavabivar$p, method="fdr")
      cat("Using FDR correction across",nrow(lavabivar), "comparisons.\n")
      
      rows <- which(lavabivar$p.fdr<0.05)
    } else { 
      #Go by p-thresholds
      if(is.null(opt$pthresh)){
        sig <- 0.05/nrow(lavabivar) 
      } else {
        sig <- opt$pthresh
      }
      rows <- which(lavabivar$p<sig)
    }
    
    #Extract additional columns according to the outFormat options
    if(opt$outFormat=="LAVAlocus"){
      cols <- c(cols,"locus")
    } else if (opt$outFormat=="genomicRegion") {
      cols <- c(cols,"chr", "start", "stop")
    }
    
    #Return p-values per the chosen setting, but override this if writing out univariate results (mainly because p-values would interfere with any deduplication)
    if(opt$returnPvalue && !opt$writeUnivariate){
      pcols<- colnames(lavabivar)[colnames(lavabivar) %in% c("p","p.fdr")]
      cols <- c(cols,pcols)
      cat("Output file will contain:", paste(pcols,collapse=" & "),"\n")
    }
    
  }
  
  if(!is.null(opt$everythingToFile)){
    write.table(lavabivar,
                file=opt$everythingToFile,
                quote=FALSE,
                append=FALSE,
                sep="\t",
                col.names = TRUE,
                row.names = FALSE)
    
  }
  
  if(length(rows)>0){ #Write to file if there are any rows to include
    
    if(!is.null(opt$fullRowsToFile)){
      write.table(lavabivar[rows,],
                  file=opt$fullRowsToFile,
                  quote=FALSE,
                  append=FALSE,
                  sep="\t",
                  col.names = TRUE,
                  row.names = FALSE)
      
    }
    
    
    toWrite <- lavabivar[rows,cols]
    
    #Widen the snp region extracted by the number indicated in opt$regionWindow
    if(all(c("start","stop") %in% cols) && opt$regionWindow!=0){
      toWrite$start <- toWrite$start-opt$regionWindow
      toWrite$stop  <- toWrite$stop+opt$regionWindow
    }
    
    #Define traits into comma-separated pairs
    toWrite[,1] <- apply(toWrite[,1:2],1,paste,collapse = ',')
    toWrite<- toWrite[,-2]
    names(toWrite)[1] <- "traits"
    
    #Optionally convert outfile into a univariate list
    if(opt$writeUnivariate){
      #Pivot out traits
      toWrite<- tidyr::separate_rows(toWrite,"traits",sep=",")
      
      #Check for and remove any duplicate rows if they occur
      deduplicate<- unique.data.frame(toWrite)
      if(nrow(deduplicate)!=nrow(toWrite)){
        cat(paste0("Deduplication prior to writing out univariate results removed ",nrow(toWrite)-nrow(deduplicate)," row(s).\n"))
        toWrite <- deduplicate
      }
    }
    
    if(opt$outFormat=="genomicRegion") {
      #Concatenate the genomic region into a comma separated list, then drop the old columns
      toWrite[,2] <- apply(toWrite[,2:4],1,paste,collapse = ',')
      toWrite<- toWrite[,-(3:4)]
      names(toWrite)[2] <- "region"
    }
    
    write.table(toWrite,
                file=opt$outfile,
                quote=FALSE,
                append=file.exists(opt$outfile),
                sep="\t",
                col.names = !file.exists(opt$outfile),
                row.names = FALSE)
  } else {
    message("No output returned for file ", basename(files[i]) ," since no loci were found with local genetic correlation p-value below the current threshold.",
        "If this is not expected, please check the --outFormat, --pthresh, and --useFDR options.")
  }
}

