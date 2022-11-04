#####
# Author: Thomas Spargo (thomas.spargo@kcl.ac.uk)
# Obtained from GitHub repository: https://github.com/ThomasPSpargo/COLOC-reporter
#
# The purpose of this script is to extract genomic regions highly correlated between two traits as produced by the LAVA software.
#####


suppressMessages(library(optparse))

option_list = list(
  make_option("--indir", action="store", default=".", type='character',
              help="The file path containing local genetic correlation files from lava. Files are extracted by selecting all of those within the directory with the extension set in --extension, which is by default '.bivar'"),
  make_option("--extension", action="store", default=".bivar", type='character',
              help="File extension used across the files to read in"),
  make_option("--outfile", action="store", default="./coloc.phenopairs.txt", type='character',
              help="The name of the file to store the pheno pairs"),
  make_option("--pthresh", action="store", default=NULL, type='numeric',
              help="Manually assign the p-value threshold for loci to extract. If not set, the p-value threshold is 0.05/nrow(<file>)"),
  make_option("--outFormat", action="store", default="genomicRegion", type='character',
              help="Character string which may be any of: 'genomicRegion', 'LAVAlocus', 'traitsOnly'.\n'genomicRegion' is the default, and extracts regions below the p-value threshold in comma separated list of the format 'chromosome,start_position,end_position' (e.g. 17,43460501,44865832).\n'LAVAlocus' functions like genomicRegion, but instead extracts the locus number associated with the region.\n'traitsOnly' ignores any p-value thresholding, and instead extracts every pair of traits for which LAVA has an output - the output will only have 2 columns.")
)
opt = parse_args(OptionParser(option_list=option_list))

suppressMessages(library(dplyr))

######
### Do some sanity checks that initial input checks
######
if(!opt$outFormat %in% c('genomicRegion', 'LAVAlocus', 'traitsOnly')){
  stop("The --outFormat option has not been set correctly; please specify one of 'genomicRegion', 'LAVAlocus', or 'traitsOnly'. Refer to the help documentation for further details.")  
}

#Make the input directory the working directory
setwd(opt$indir)

# #Check that the options are correctly called - maybe not necessary, people might want to use input directoy
# if(opt$indir=="."){
#   warning(paste0("--indir is not specified, looking for LAVA output files with the extension '",opt$extension,"' the current directory:\n",getwd()))  
# }

######
### Identify files, if any are identified set output file and loop across them to identify regions to extract
######
files<-list.files(path=opt$indir,pattern=opt$extension)

if(length(files)>0){
  files <- paste0(opt$indir,files)
  file.create(opt$outfile)   #Create the output file
} else {
  stop("No files with the extension '",opt$extension,"' were identified in the  directory:\n",getwd(),"\nPlease check documentation for the --indir and --extension options.")
}

for(i in 1:length(files)){
  lavabivar <- read.table(files[i],header=T)
  
  #Set cols to extract in all settings of outFormat
  cols <- c("phen1","phen2")
  
  if(opt$outFormat=="traitsOnly"){
    rows <- 1 #Extract only row 1, to get each trait pairwise
  } else {
    
    #Identify rows according to the p-value threshold. If no p<sig then no rows will be identified, and no output will be written
    if(is.null(opt$pthresh)){
      sig <- 0.05/nrow(lavabivar) 
    } else {
      sig <- opt$pthresh
    }
    rows <- which(lavabivar$p<sig)
    
    #Extract additional columns according to the outFormat options
    if(opt$outFormat=="LAVAlocus"){
      cols <- c(cols,"locus")
    } else if (opt$outFormat=="genomicRegion") {
      cols <- c(cols,"chr", "start", "stop")
    }
    
  }
  
  if(length(rows)>0){ #Write to file if there are any rows to include
    
    toWrite<- lavabivar[rows,cols]
    
    if(opt$outFormat=="genomicRegion") {
      #Concatenate the genomic region into a comma separated list, in the column 'locus', then drop the chr, start, stop columns
      toWrite$locus <- apply(toWrite[,3:5],1,paste,collapse = ',')
      toWrite<- toWrite[,-(3:5)]
    }
    
    write.table(toWrite,
                file=opt$outfile,
                quote=FALSE,
                append=TRUE,
                sep=" ",
                col.names = FALSE,
                row.names = FALSE)
  } else {
    message("No output returned for file ", basename(files[i]) ," since no loci were found with local genetic correlation p-value below the current threshold.",
        "If this is not expected, please check the --outFormat and --pthresh options.")
  }
}

