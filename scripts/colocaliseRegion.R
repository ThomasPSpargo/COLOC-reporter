#####
# Author: Thomas Spargo (thomas.spargo@kcl.ac.uk)
# Obtained from GitHub repository: https://github.com/ThomasPSpargo/COLOC-reporter
#
# This script performs colocalisation analysis between two traits at a specified genomic region. 
# Preprocessing and finemapping are handled per-trait and colocalisation analysis is performed with the COLOC software, either using coloc.abf or coloc.susie.
#####
suppressPackageStartupMessages(library(optparse))

option_list = list(
  
  #General options
  make_option("--plink", action="store", default='plink', type='character',
              help="Path to PLINK executable (syntax written for PLINK 1.9). By default, has the value 'plink'"),
  make_option("--traits", action="store", default=NULL, type='character',
              help="Comma separated list of IDs for traits to analyse (e.g. 'P1,P2'), as relevant to the summary statistics defined in the  --GWASconfig input file."),
  make_option("--set_locus", action="store", default=NA, type='character',
              help='Specify target locus as a comma separated list of the format "chromosome,start_position,end_position" (e.g. 17,43460501,44865832).'),
  make_option("--LDreference", action="store", default=NA, type='character',
              help="Path to, and prefix for, per-chromosome PLINK binary files used to compute LD matrix for SNPs in region. This input is used by SuSiE. Expected format of <prefix>i.bim, where i is the chromosome number"),
  make_option("--GWASconfig", action="store", default=NULL, type='character',
              help="Path to file specifying configuration of GWAS sumstats. Configuration will be determined by identifying rows whose IDs match the --traits option."),
  make_option("--runMode", action="store", default="trySusie", type='character',
              help="Character string, any of 'trySusie', 'skipSusie', 'doBoth','finemapOnly'.\nIf 'doBoth', both coloc.abf and coloc.susie will attempt to run.\nIf 'trySusie' coloc.susie will be used if SuSiE finemapping identifies at least 1 credible set in each trait and coloc.abf is returned otherwise.\nIf 'skipSusie', only coloc.abf will be applied, and processes necessary for coloc.susie are skipped (e.g. no need to call to plink and generate LD matrix); the LD reference will however still still used for SNP alignment.\nIf 'finemapOnly' then univariate finemapping is performed which can be applied across any number of traits. Note however that only snps in common across all traits and the reference panel will be retained. Therefore, for colocalisation analysis it may be preferable to harmonise summary statistics used in colocalisation analysis pairwise."),
  make_option("--force_matrix", action="store", default=FALSE, type='logical',
              help="If TRUE, the LD matrix will always be recomputed using PLINK. If FALSE, the default, the LD matrix will only be computed if the expected ld_matrix.snplist and ld_matrix.ld files are absent from the <output>/LDmatrix directory."),
  make_option("--finemapQC_handleBetaFlips", action="store", default="drop", type='character',
              help="The SuSiE kriging_rss function is used to check whether observed Z-scores (calculated from beta and SE) correspond with expected Z-scores based on information from the LD matrix and other SNPs. This check indicates if any SNP effect estimates appear to be reversed. Set this option to 'flip' to automatically reverse the direction of betas with test statistics identified as inverted on a trait-by-trait basis, or to 'drop' to remove them from all traits across analyses; the option defaults to 'drop'. The check for potentially flipped alleles will always be performed at least once but the data used for subsequent analysis will be unchanged if any other strings are passed to this option. Note that the setting applied here will affect both the finemapping and colocalisation analysis steps, but will only be applied in circumstances when SuSiE finemapping is performed. Note also that the recursion behavior of this check is adjusted using the --finemapQC_limitCheckBetaFlips option."),
  make_option("--finemapQC_limitCheckBetaFlips", action="store", default=10, type='numeric',
              help="Related to the --finemapQC_handleBetaFlips option. Define the max number of recursive checks for reversed statistic encoding before breaking the loop. This arbitrarily defaults to 10, a limit which seems unlikely to be hit in a dataset with largely correct allele encoding. Setting 0 means that the check will be performed once only."),
  make_option("--finemap_CScoverage", action="store", default=0.95, type='numeric',
              help="Numeric between 0 and 1 to indicate the credible set threshold to use for finemapping; defaults to 0.95, which returns 95% credible sets."),
  make_option("--finemap_refine", action="store", default=FALSE, type='logical',
              help="Logical, defaulting to FALSE. Specify TRUE to add a refinement step to SuSiE finemapping call to avoid identification of local maxima."),
  make_option("--finemap_initialL", action="store", default=10, type='numeric',
              help="Numeric, defaulting to 10, the susie default, to indicate the maximum number of non-zero effects to allow in the initial susie regression model."),
  make_option("--finemap_increaseL", action="store", default=TRUE, type='logical',
              help="Logical, defaulting to TRUE which indicates that higher values of L should be attempted in susie finemapping (with L+10 for each loop) when the number of credible sets found match the current L value. Setting this to FALSE will run susie once only, using the L setting specified in the --finemap_initialL option."),
  make_option("--restrict_source_for_nearest_finemap_gene", action="store", default=NULL, type='character',
              help="Comma delimited character string indicating sources to consider when identifying the nearest gene to the top SNP identified within any finemapping credible sets. Included as an option to prioritise tagging of nearby genes which are e.g. HGNC symbols. If NULL, the default, all sources will be evaluated."),
  make_option("--priors_susie", action="store", default=NULL, type='numeric',
              help="Set a prior probability that a snp is causal for susie finemapping (equivalent to coloc p1 or p2). Set to NULL by default, the default for the coloc::runsusie wrapper around susie_rss; this is not usual default for susie_rss - see coloc::runsusie documentation for details."),
  make_option("--priors_coloc.abf", action="store", default="1e-4,1e-4,1e-5", type='character',
              help="Comma separated list of 3 numerics indicating priors to set respectively for p1,p2,p12 in coloc.abf function; default values are '1e-4,1e-4,1e-5' the coloc.abf default settings."),
  make_option("--priors_coloc.susie", action="store", default="1e-04,1e-04,5e-06", type='character',
              help="Comma separated list of 3 numerics indicating priors to set respectively for p1,p2,p12 in coloc.susie function; default values are '1e-04,1e-04,5e-06' the coloc.susie default settings (passed through to coloc.bf_bf function)."),
  make_option("--out", action="store", default="./COLOC-reporter", type='character',
              help="Path and prefix for directory in which to return all outputs. Defaults to ./COLOC-reporter. When running multiple analyses, unique output directories are essential for tidy file organisation."),
  make_option("--rdsOut",action="store", default=NULL, type='character',
              help="Specify a directory in which to return Rds files containing r markdown report parameters for use in generating a global analysis report [optional]"),
  make_option("--genomeAlignment", action="store", default=37, type='numeric',
              help="Indicate the reference genome to which summary statistics and LD reference are aligned. Defaults to 37."),
  make_option("--gene_tracks", action="store", default=NULL, type='numeric',
              help="Numeric. Specify the threshold at which plotting of nearby genes changes from one gene per row to plotting multiple genes across individual rows, labelling genes with geom_text_repel. Defaults to NULL, where per-row plotting has no upper limit. Specify 1 to always plot multiple genes per row, provided their genomic positions do not overlap [work in progress]."),
  make_option("--restrict_nearby_gene_plotting_source", action="store", default=NULL, type='character',
              help="Comma delimited character string indicating sources from which nearby genes included in plots are retrieved. Check output tables 'external_gene_source' column for options. Included as an option to reduce overplotting of obscure gene symbols. If NULL, the default, all sources will plot"),
  make_option("--scriptsDir", action="store", default=NULL, type='character',
              help="Filepath to directory script and subdirectory dependencies (note: expects the subdirectories 'helper_functions' (containing custom functions) and 'rmd' (for rmarkdown report templates)")
)

opt = parse_args(OptionParser(option_list=option_list))

#######
### Begin script
#######

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table) #To read in the datasets
  library(R.utils)    #Required reading gz files directly with fread
  library(coloc)
  library(susieR)
  library(biomaRt)    #For ensembl library
  library(ggrepel)
  library(patchwork)     #arranging summary plot
})

test <- FALSE # test <- TRUE
if(test){
  #####
  setwd("/Users/tom/OneDrive - King's College London/PhD/PhD project/COLOC/git.local.COLOC-reporter/testing")
  
  # # #TESTING V3
  opt$set_locus <- "17,43460501,44865832"
  opt$GWASconfig <- "./GWAS_samples_testing.txt"
  opt$out <- "./tidy_processingTEST_PD.SZ.chr17"
  opt$traits <- "PD,SZ"#,ALS"
  
  opt$LDreference <- './EUR_phase3_chr'
  opt$runMode <- 'doBoth'
  opt$force_matrix <- FALSE
  
  opt$gene_tracks <- 40
  opt$restrict_source_for_nearest_finemap_gene <- opt$restrict_nearby_gene_plotting_source <- "HGNC Symbol"
  opt$genomeAlignment <- 37
  opt$rdsOut <-"./tidy_processingTEST_PD.SZ.chr17_coloc/data/globalReportInputs"
  
  opt$scriptsDir <- "/Users/tom/OneDrive - King's College London/PhD/PhD project/COLOC/git.local.COLOC-reporter/scripts"
  
}
P01<-list(saveopts=opt)

#Extract names of traits compared, first dropping file path, and then any prefixes indicated by a preceding underscore
P01$traits <- traits <- strsplit(opt$traits,",")[[1]]

coveragePcent<- paste0(opt$finemap_CScoverage*100,"%")

#Read in the configuration options
GWASconfig<- fread(opt$GWASconfig)
Config_msg <- paste0("Summary statistic configuration options set according to the specification of --GWASconfig for trait IDs: ", paste0(1:length(traits)," = ",traits,collapse = ", "),".\n",sep='')

#Run some simple config checks
checkdups<- duplicated(GWASconfig$ID)
if(any(checkdups)){ stop("Terminating script because the file supplied to --GWASconfig contained duplicates in the ID column. Please assign unique IDs for all traits. Problem found with: ", paste0(GWASconfig$ID[checkdups],collapse=", "),"\n"); q("no") }
checkunknownID<- !traits %in% GWASconfig$ID
if(any(checkunknownID)){ stop("Terminating script because the IDs indicated in --traits were not found in the ID column of the file supplied to --GWASconfig. Problem found with: ",paste0(traits[checkunknownID], collapse=", "),"\n"); q("no") }
rm(checkunknownID,checkdups)

#Assign names to traits
names(traits) <- GWASconfig[traits,on="ID"]$traitLabel

#Set up directories across which outputs are returned
opt$out <- paste0(opt$out,"_coloc")
if(!dir.exists(opt$out)){dir.create(opt$out,recursive = TRUE)}

figdir <- file.path(opt$out,"plots")
if(!dir.exists(figdir)){dir.create(figdir)}
tabdir <- file.path(opt$out,"tables")
if(!dir.exists(tabdir)){dir.create(tabdir)}
datadir <- file.path(opt$out,"data")
if(!dir.exists(datadir)){dir.create(datadir)}
reportdir <- file.path(datadir,"reports")
if(!dir.exists(reportdir)){dir.create(reportdir)}
if(opt$runMode %in% c("trySusie", "doBoth","finemapOnly")){
  finemapQCdir <- file.path(opt$out,"finemapQC")
  if(!dir.exists(finemapQCdir)){dir.create(finemapQCdir)}
}

if(!is.null(opt$rdsOut) && !dir.exists(opt$rdsOut)) dir.create(opt$rdsOut,recursive = TRUE)

logfile <- file.path(opt$out,'colocalisation.log')

#Setup a colour palette to be used across ggplots (colours taken from plot SuSiE)
ggpalette = c("dodgerblue2", "green4", "#6A3D9A", "#FF7F00", 
              "gold1", "skyblue2", "#FB9A99", "palegreen2", "#CAB2D6", 
              "#FDBF6F", "gray70", "khaki2", "maroon", "orchid1", "deeppink1", 
              "blue1", "steelblue4", "darkturquoise", "green1", "yellow4", 
              "yellow3", "darkorange4", "brown")

sink(file = logfile, append = F)
cat(
  '#################################################################
# colocaliseRegion.R
# Perform colocalisation analysis between GWAS summary statistics in a predetermined genomic region.
# Colocalisation leverages the COLOC and SuSiE software packages.
# Genetic regions must either identified manually or can be selected based on output from LAVA software.
# Documentation for options for running the script can be called with: Rscript ./colocaliseRegion.R --help
# 
# This script was written by Thomas Spargo (thomas.spargo@kcl.ac.uk), please get in touch with any queries.
#
#################################################################
Analysis started at',as.character(Sys.time()),'\n',Config_msg,'Options are:\n')
print(opt)

cat("\n######\n### Setup\n######\n\n")

cat("All outputs will be returned in the directory: ", opt$out,"\n")
cat("Figures from the main analysis steps are returned within the subdirectory: ",basename(figdir),"\n")
cat("Tabular summaries are in the subdirectory: ",basename(tabdir),"\n")
if(opt$runMode %in% c("trySusie", "doBoth","finemapOnly")){
  cat("Finemap quality control checks are in the subdirectory: ",basename(finemapQCdir),"\n")
}
cat("Resources from each main analysis step are in the subdirectory: ",basename(datadir),"\n\n")

sink()

#Read-in custom helper functions stored in directory specified by opt$scriptsDir
list.files(file.path(opt$scriptsDir,"helper_functions"),full.names = TRUE,pattern=".R") %>%
  lapply(.,source) %>%
  invisible(.)

  
######
### Import datasets
######
#Import datasets and convert to tibbles
sums <- lapply(GWASconfig[traits,on="ID"]$FILEPATH,function(x){tibble(fread(x))})
names(sums) <- traits

#Loop across traits to test whether column names match those indicated in the config file
colcheck<- mapply(function(x,trait){

  #Define mandatory columns
  setCols<- c(GWASconfig[trait,on="ID"]$p_col,
                GWASconfig[trait,on="ID"]$stat_col,
                GWASconfig[trait,on="ID"]$N_col,
                GWASconfig[trait,on="ID"]$se_col,
                GWASconfig[trait,on="ID"]$snp_col,
                GWASconfig[trait,on="ID"]$A1_col,
                GWASconfig[trait,on="ID"]$A2_col,
                GWASconfig[trait,on="ID"]$freq_col)
  
  #If not NA, check for chr and pos
  if(!is.na(GWASconfig[trait,on="ID"]$pos_col)) setCols <- c(setCols,GWASconfig[trait,on="ID"]$pos_col)
  if(!is.na(GWASconfig[trait,on="ID"]$chr_col)) setCols <- c(setCols,GWASconfig[trait,on="ID"]$chr_col)
  
  colcheck <-  which(!(setCols%in% colnames(x)))
  return(colcheck)
},x=sums,trait=traits,SIMPLIFY = FALSE)

colProblem<- sapply(colcheck,length)>0
if(any(colProblem)){
  sink(file = logfile, append = T)
  warning("The column names expected based on the --GWASconfig options set for trait(s): ", paste0(names(colProblem)[colProblem],collapse=", "), " do not match columns detected in the file. Please check the option(s) specified.")
  sink()
}

#Adjust column names into COLOC format.
#mapply is set to SIMPLIFY=FALSE to ensure a list result
sums<- mapply(function(x,trait){
  #Unless specified to be an odds ratio, rename as beta
  if(GWASconfig[trait,on="ID"]$stat_col=="OR"){
    x$beta <- log(x$OR)
    x$OR <- NULL
  } else {
    names(x)[names(x)==GWASconfig[trait,on="ID"]$stat_col] <- "beta"
  }
  
  #Rename the other columns  
  names(x)[names(x)==GWASconfig[trait,on="ID"]$p_col] <- "pvalues"
  names(x)[names(x)==GWASconfig[trait,on="ID"]$N_col] <- "N"
  names(x)[names(x)==GWASconfig[trait,on="ID"]$chr_col] <- "chr"
  names(x)[names(x)==GWASconfig[trait,on="ID"]$pos_col] <- "pos"
  names(x)[names(x)==GWASconfig[trait,on="ID"]$se_col] <- "SE"
  names(x)[names(x)==GWASconfig[trait,on="ID"]$snp_col] <- "snp"
  names(x)[names(x)==GWASconfig[trait,on="ID"]$A1_col] <- "A1"
  names(x)[names(x)==GWASconfig[trait,on="ID"]$A2_col] <- "A2"
  names(x)[names(x)==GWASconfig[trait,on="ID"]$freq_col] <- "MAF"
  
  #Filter down to recognised cols only
  x <- x[colnames(x) %in% c("pvalues","N","chr","pos","beta","SE","snp","A1","A2","MAF")]
    
  return(x)
},x=sums,trait=traits,SIMPLIFY=FALSE)


########
#### Identify region for analysis
########
if(!is.na(opt$set_locus)){
  opt$set_locus <- as.numeric(strsplit(opt$set_locus,",")[[1]])
  
  if(length(opt$set_locus)==3){
    #Set genomic positions manually
    reg_range <-c(chr=opt$set_locus[1],start=opt$set_locus[2],stop=opt$set_locus[3])
  } else {
    stop('The option --set_locus is defined but cannot be identified as specifying a genomic region. The expected format is a comma separated list of "chromosome,start_position,end_position"')
  }
} else {
  stop('Please define a genomic region with the --set_locus option. This expects a comma separated list of numerics in the format "chromosome,start_position,end_position"')
}

P01$target_region <- paste0("Chr",reg_range["chr"],":",reg_range["start"],"-",reg_range["stop"])

sink(file = logfile, append = T)
cat(paste0("Colocalisation analysis will be performed for the region: ", P01$target_region," \n"))
sink()

######
# Harmonise summary statistics with reference
######

# Define function to flip sumstat direction to match alleles across gwas and reference,
# Allele frequency is not reversed since coloc utilises MAF.
alignSS <- function(ss,bim){
  ss_bim_match<-merge(ss, bim, by.x=c('snp','A1','A2'), by.y=c('V2','V5','V6'))
  ss_bim_swap<-merge(ss, bim, by.x=c('snp','A1','A2'), by.y=c('V2','V6','V5'))
  
  ss<-ss[ss$snp %in% ss_bim_match$snp | ss$snp %in% ss_bim_swap$snp,] 
  ss$beta[ss$snp %in% ss_bim_swap$snp] <- -ss$beta[ss$snp %in% ss_bim_swap$snp]
  
  return(ss)
}

# Read in reference SNP data to match alleles 
bim<-fread(paste0(opt$LDreference,reg_range["chr"],'.bim'))
names(bim)[c(1,4)] <- c("chr","pos") #Assign names for chromosomal positions in case of matching

### Check for missing chromosomal positions column info from the GWAS config file
sums<- mapply(function(ss,trait){
  if(any(is.na(GWASconfig[trait,on="ID",c(pos_col,chr_col)]))){
    ss<- inner_join(ss[,!names(ss) %in% c("chr","pos")],
               bim[,c("V2","chr","pos")],
               by=c('snp'= "V2"))
  }
  return(ss)
},ss=sums,trait=traits,SIMPLIFY=FALSE)




#Return total snps available
P01$initial_data_overview<- lapply(sums,function(x){
  y <- x[x$chr==reg_range["chr"] &
         x$pos>=reg_range["start"] &
         x$pos<=reg_range["stop"],]
  datasetSummary(y)}) %>%
  do.call(rbind,.) %>%
  cbind(traits,.)

sink(file = logfile, append = T)
cat("Total number of SNPs for each sumstats in tested region:\n", paste0(P01$initial_data_overview$traits," = ",P01$initial_data_overview$NSNP,collapse = "\n"),"\n\n",sep='')
sink()


#Filter to snps in the define chromosome:position range and present in the LD reference data
#Mutate to calculate minor allele frequency and varbeta.
sums.region <- lapply(sums, function(x){
  x %>%
    filter(chr==reg_range["chr"],
           pos>=reg_range["start"],
           pos<=reg_range["stop"],
           snp %in% bim[["V2"]]
    ) %>%
    arrange(match(snp, bim[["V2"]])) %>%
    mutate(MAF=if_else(MAF>0.5,1-MAF,MAF),
           varbeta = SE^2) %>%
    alignSS(.,bim)
})

#Generate a quick data summary
P01$ref_harmonised_data_overview<- lapply(sums.region,datasetSummary) %>%
  do.call(rbind,.) %>%
  cbind(traits,.)

#Return snps available after harmonising with sumstats
snpsavail<- lapply(sums.region,nrow)
sink(file = logfile, append = T)
cat("Number of SNPs in common between LD reference and each set of sumstats in tested region:\n", paste0(names(snpsavail)," = ",snpsavail,collapse = "\n"),"\n",sep='')
sink()
rm(snpsavail)

#Intersect the snp lists to find those common across traits
snplist <- Reduce(intersect,lapply(sums.region,function(x){x$snp}))
sink(file = logfile, append = T)
#Print some information to console
cat(paste0("N SNPs in common across traits after harmonising to LD reference: ", length(snplist),"\n"))
if(length(snplist)==0){cat("Analysis stopped as no SNPs remain.\n");stop("Analysis stopped as no SNPs remain.")}
sink()



#Intersect the overlapping SNPs, and assign sequence position in common snp sequence (based on rownumber).
sums.region <- lapply(sums.region, function(x){x %>%
    filter(snp %in% snplist) %>%
    mutate(position = row_number())
})

snpid<- lapply(sums.region,function(x){paste0(x$snp,"_",x$pos)})
#Return warning if the SNP position alignment seems incorrect.
#sapply compares for identical results across all vectors relative the first vector
if(!all(sapply(snpid[-1],identical, snpid[[1]]))){
  sink(file = logfile, append = T)
  cat("WARNING: SNP names and genomic positions do not match between all datasets. Please check that these have been aligned correctly.")
  sink()
}
rm(snpid)

######
# Generate LD matrix (for SuSiE)
######

if(opt$runMode %in% c("trySusie", "doBoth","finemapOnly")){
  #Extract snps and write list to file in subdirectory of results directory
  
  #Check for the expected LD matrix output; if files are absent or if --force_matrix is set, generate using PLINK
  expect <- c(".snplist",".ld")
  expect <- file.path(datadir,'LDmatrix',paste0('ld_matrix',expect))
  
  if(any(!file.exists(expect)) | opt$force_matrix){
    
    sink(file = logfile, append = T)
    cat("Computing LD matrix...")
    sink()
    
    ld_dir<- dirname(expect[1])
    
    if(!dir.exists(ld_dir)){dir.create(ld_dir,recursive = TRUE)}
    write(snplist, file=file.path(ld_dir,"snplist.txt"))
    
    #Identify LD matrixsnps from snplist present in reference
    #Syntax based on plink v 1.9
    system(paste0(opt$plink,' --bfile ', opt$LDreference,reg_range["chr"],
                  ' --extract ',ld_dir,'/snplist.txt',
                  ' --r square',
                  ' --write-snplist',
                  ' --keep-allele-order',
                  ' --out ',file.path(ld_dir,'ld_matrix')))
    
    sink(file = logfile, append = T)
    cat("Done\n")
    sink()
    
  } else {
    sink(file = logfile, append = T)
    cat("Existing LD matrix found, skipping call to PLINK.\n")
    sink()
  }
  
  #Read in the ld matrix and snp names as returned by plink 
  ld <- as.matrix(fread(expect[2]))
  ld_names<- scan(expect[1],what=character())
  dimnames(ld)<-list(ld_names, ld_names) #assign dimnames to LD object
  
  #As a sanity check, redo alignment based on plink output order to ensure snps are correctly arranged in base pair order, and assign 'positions' for coloc
  sums.region <- lapply(sums.region, function(x){
    x %>%
      filter(snp %in% ld_names) %>%
      arrange(match(snp, ld_names)) %>%
      mutate(position = row_number())
  })

} 

######
# Generate a minimal dataset across traits which can be readily used for generation of summary plots
# This will be saved to file following the SuSiE steps where a column for credible sets may be appended
######
minimal_df<- mapply(function(x,trait){
  
  mindf<- cbind(x[,c("snp","A1","A2","chr","pos","beta","SE","pvalues")],trait) #Subset to useful columns
  mindf$z <- mindf$beta/mindf$SE #Calculate z-scores
  
  return(mindf)
},sums.region,traits, SIMPLIFY = FALSE) %>%
  do.call(rbind,.) 


######
### Prepare for colocalisation analysis; convert each data.frame to list and add expected list elements
######

##Expected elements:
#LD matrix (if passing via SuSiE)
#type = "cc" or "quant" depending on trait type
#s = case control proportion if "cc"
#Adjust column names into COLOC format.
#mapply is set to SIMPLIFY=FALSE to ensure a list result
sums.region<- mapply(function(x,trait){
  
  x <- as.list(x)
  x$type <- tolower(GWASconfig[trait,on="ID"]$type)
  if(x$type=="cc"){x$s <- GWASconfig[trait,on="ID"]$prop
  } else if(x$type=="quant" && !is.na(GWASconfig[trait,on="ID"]$traitSD)){x$sdY <- GWASconfig[trait,on="ID"]$traitSD}
  if(opt$runMode %in% c("trySusie", "doBoth","finemapOnly")){x$LD <- ld}
  
  #Save the formatted dataset as a resource
  datasetpath<- file.path(datadir,"datasets")
  if(!dir.exists(datasetpath)){dir.create(datasetpath)}
  saveRDS(x,file=file.path(datasetpath,paste0("coloc_format_dataset_",trait,".Rds")))
  
  return(x)
},x=sums.region,trait=traits,SIMPLIFY=FALSE)

######
### Generate a final summary of data and write out the pre-processing report
######

#Generate a quick data summary
P01$final_data_overview<- lapply(sums.region,datasetSummary) %>%
  do.call(rbind,.) %>%
  cbind(traits,.)

if(!is.null(opt$rdsOut)) saveRDS(P01,file.path(opt$rdsOut,"P01.Rds")) #Optionally save P01 to an Rds file for global report generation


report01<- file.path(normalizePath(reportdir),paste0("01_preprocessing_report.html"))
renderReport(params="P01",
             outfile=report01,
             template=file.path(opt$scriptsDir,"rmd","preprocess_report.Rmd"))

######
### Extract genes around tested region
######

#Import gene information for the chromosome.
#Cache set to false because conflicts were encountered frequently.
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=opt$genomeAlignment)
Genes<-getBM(attributes=c('external_gene_name','chromosome_name','start_position','end_position','external_gene_source','strand'),
             filters ='chromosome_name', values =reg_range[["chr"]], mart = ensembl,useCache=FALSE)

# Use 10kb window to define gene window
gene_window<-10000
Genes$start_window<-Genes$start_position-gene_window
Genes$end_window<-Genes$end_position+gene_window

# Do main Susie steps
if(opt$runMode %in% c("trySusie", "doBoth","finemapOnly")){
  
  ## Filter genes in chromosome for matching any credible sets to nearest up/downstream genes
  genes_checkUpOrDownstream<- Genes
  if(!is.null(opt$restrict_source_for_nearest_finemap_gene)){ #Filter down the gene sources to match the option set (if any)
    src <- strsplit(opt$restrict_source_for_nearest_finemap_gene,split=",")[[1]] #String split on comma to identify sources
    
    if(sum(genes_checkUpOrDownstream$external_gene_source %in% src)==0){ warning(paste("Problem with --restrict_source_for_nearest_finemap_gene option. could not filter to the indicated source(s) since no genes remained in chromosome after restricting to:", paste0(src,collapse=", ")))
    } else { genes_checkUpOrDownstream<- genes_checkUpOrDownstream[genes_checkUpOrDownstream$external_gene_source %in% src,] }
  }
  
  sink(file = logfile, append = T)
  cat("\n######\n### SuSiE finemapping quality control\n######\n")
  sink()
  
  ## Add initial elements to the finemapping report params lists.
  
  #The number of SNPs in the datasets prior to finemapping QC, the cs coverage, genome alignment, the name of each trait.
  #To avoid stripping away trait names, assign traits within a loop rather than with mapply
  for(i in seq_along(traits)){
    assignToList(x=paste0("RMD_finemap_",traits[i]),multi=TRUE,
                 value=list(trait=traits[i],
                            initNsnp=length(sums.region[[i]]$snp), #initNsnp should always be fully harmonised across traits
                            coverage=coveragePcent,
                            alignment=opt$genomeAlignment)
    )
  }
  
  susieQC <- mapply(function(dset,trait){
      ## Compare observed z-scores with LD matrix [warning thrown when LD matrix is not semi-definite...]
      z_sc <- dset$beta/dset$SE #Get z-scores
      
      # #Estimate s with all available methods - indicates compatibility between sumstats and reference
      checkMthd<- c("null-mle", "null-partialmle", "null-pseudomle")
      checkLD<- sapply(checkMthd,estimate_s_rss,z=z_sc, R=dset$LD, n=median(dset$N), r_tol = 1e-08)
      
      #Add extra info for writing to file
      checkTofile<- data.frame(t(c("trait"= trait, "n_flips"=0,round(checkLD,5))))
      
      ### Check for allele flip issues
      #Compare observed and expected Z-scores
      z_compare<- kriging_rss(z=z_sc, R=dset$LD, n=median(dset$N), r_tol = 1e-08,s = checkLD[[1]])
      zPlot<- list("0 flips"=z_compare$plot)
      
      #Determine snps with betas that may be flipped; this check corresponds with the one implemented by internally by kriging_rss
      #possFlipped stores the number of flips detected in a single recursion of the check
      #flipIndex stores a list across all recursions which will later be reduced down to provide a final index
      #nPossFlipped stores a cumulative record of the number of SNPs recorded as flipped
      possFlipped<- with(z_compare$conditional_dist,
                         logLR>2 & abs(z)>2)
      flipIndex <- list(possFlipped)
      nPossFlipped<- sum(possFlipped)
      
      #Recursively loop (according to the setting of --finemapQC_limitCheckBetaFlips) to test for inverted summary statistic estimate issues
      #The repeated check is recommended in the Zou et al. 2022 paper to ensure all flips are identified
      if(opt$finemapQC_limitCheckBetaFlips>0){
        nloops <- 0 #Flag variable for limiting the while loop
        while(any(possFlipped)){ #If any outliers replot the figure in a circumstance where these are flipped
          nloops <- nloops+1
          
          z_sc[possFlipped]<- -z_sc[possFlipped] #Flip the flagged indices
          
          #Re-compare observed and expected Z-scores
          z_compare_flip <- kriging_rss(z=z_sc,R=dset$LD,n=median(dset$N))
          
          #Determine SNPs with betas that are currently flagged for potential flips and append to the flip index list
          possFlipped<- with(z_compare_flip$conditional_dist,
                             logLR>2 & abs(z)>2)
          flipIndex <- c(flipIndex,list(possFlipped))
          
          #Save elements (the cumulative number of flipped snps, and the new flip plot)
          #The plot name details the number of flips BEFORE plotting, nPossFlipped identifies the number as a result of the plot
          zPlot <- c(zPlot,list(z_compare_flip$plot))
          names(zPlot)[length(zPlot)] <- paste(max(nPossFlipped),"flips")
          
          nPossFlipped[length(nPossFlipped)+1]<- sum(possFlipped)+nPossFlipped[length(nPossFlipped)]
          
          if(any(possFlipped) && opt$finemapQC_limitCheckBetaFlips==nloops){
            warning("The recursion limit defined in --finemapQC_limitCheckBetaFlips was hit during the repeated observed and expected Z-score consistency check for ",trait,". If this is unexpected please review the diagnostic plots returned in the .html report or finemapQC directory.")
            break
          }
        }
      }
      
      if(any(rowSums(as.data.frame(flipIndex))>1)){
        warning("During the repeated observed and expected Z-score consistency check for ",trait,", at least one variant was flagged more than once (i.e. after the first flip it was again identified as an outlier). Please check the data carefully.")
      }
      
      #Reduce the flip index into a single logical vector
      flipIndex <- Reduce(`|`, flipIndex)
        
      #Save the plot
      z_check_outpath <- file.path(finemapQCdir,paste0("Z_score_alignment_",trait,".pdf"))
      pdf(z_check_outpath,height=6,width=6)
      for (i in seq_along(zPlot)) print(zPlot[[i]]+labs(subtitle=paste0("Check performed following a cumulative n = ",names(zPlot[i]))))
      dev.off()
      
      if(max(nPossFlipped)>0){ # If any flips repeat s-estimation after flips
        checkLD<- sapply(checkMthd,estimate_s_rss,z=z_sc, R=dset$LD, n=median(dset$N), r_tol = 1e-08)
        
        #Add extra info for writing to file
        checkTofile<- rbind(checkTofile,
                            data.frame(t(c("trait"= trait, "n_flips"=max(nPossFlipped),round(checkLD,5))))
        )
      }
      
      # #Write finemapping summary to a file and write without names if file exists already
      LDcheckFile<- file.path(finemapQCdir,"check_sumstat_LDconsistency.csv")
      write.table(checkTofile,
                  file=LDcheckFile,sep = ",",row.names=FALSE,col.names = !file.exists(LDcheckFile),append = file.exists(LDcheckFile))
      
      
      #Add multiple elements to report:
      #Z-plot 
      #LD 's' consistency check; this could have one or two rows according to whether any potential encoding flips were flagged
      #Number of possible allele flips
      assignToList(x=paste0("RMD_finemap_",trait),multi=TRUE,
                   value=list(nPossFlipped=max(nPossFlipped),
                              init_s_estimate=checkTofile,
                              zPlot=zPlot))
    
      #Report relevant messages to the logfile
      sink(file = logfile, append = T)
      cat("\nQuality control for ",trait,":\n",sep="")
      cat(max(nPossFlipped),ifelse(max(nPossFlipped)==1," SNP was"," SNPs were")," flagged for potentially flipped allele encoding", ifelse(max(nPossFlipped)>0," (marked red on the diagnostic plot returned in the finemapQC directory).\n",".\n"),sep="")
      sink()
     
      return(list(conditional_dist=z_compare$conditional_dist,toFlip=flipIndex))
      
  },sums.region,traits,SIMPLIFY = FALSE)
  
  sink(file = logfile, append = T)
  cat("\nGlobal quality control:\nTests of consistency between each set of sumstats and LD matrix by each method available in SuSiE estimate_s_rss function are written to: check_sumstat_LDconsistency.csv\n")
  cat("Comparisons between observed and expected z-scores are visualised per trait in the file(s): Z_score_alignment_<trait>.pdf. If SuSiE identifies any credible sets, these data will be replotted with colouring to mark credible sets assigned.\n\n")
  sink()
  
  #Across traits, identify if any positions betas are indicated for flipping.
  #If flipping, this will be handled per-trait, if dropping, this will be done by global index
  globalFlips <- Reduce(`|`, lapply(susieQC,function(x) x$toFlip))
  
  #If set, adjust each sumstats to remove or flip outlier records
  if(opt$finemapQC_handleBetaFlips %in% c("drop","flip") && any(globalFlips)){
    
    #Logical tests indicating how allele flipping steps should proceed (i.e. with flipping, or with dropping)
    doFlipStep <- tolower(opt$finemapQC_handleBetaFlips)=="flip" && any(globalFlips)
    doFlipDROPStep <- tolower(opt$finemapQC_handleBetaFlips)=="drop" && any(globalFlips) && any(!globalFlips)
    
    #Message about which steps are to be conducted
    sink(file = logfile, append = T)
    if(doFlipStep){ cat("Summary statistic betas will be flipped per-trait for all SNPs indicated to have reversed allele order.\n")
    } else if(doFlipDROPStep){ cat("SNPs identified as likely to have flipped test statistics in any trait will be removed...\n") }
    sink()
    
    #Across datasets flip betas/drop SNPs according to settings
    sums.region <- mapply(function(dset,QC,trait,flipDrop){
      
      #According to doFlipStep and doFlipDROPStep logicals flip or drop problematic beta coefficients.
      #Flipping is handled on a trait-by-trait basis while dropping is global to keep datasets harmonised
      if(doFlipStep){ dset$beta[QC$toFlip] <- -dset$beta[QC$toFlip]
      } else if(doFlipDROPStep){ dset <- dropSNPs(dset,flipDrop) }
      
      #Write out the dataset post-qc, overwriting the original
      #Save the formatted dataset as a resource
      datasetpath<- file.path(datadir,"datasets")
      if(!dir.exists(datasetpath)){dir.create(datasetpath)}
      saveRDS(dset,file=file.path(datasetpath,paste0("coloc_format_dataset_",trait,".Rds")))
      
      return(dset)
    }
    ,sums.region,susieQC,traits,MoreArgs = list(flipDrop=globalFlips),SIMPLIFY = FALSE)
    
    #Save IDs for SNPs that have been dropped to file (if any)
    isDropped<- !snplist %in% sums.region[[1]]$snp
    if(any(isDropped)){
      write(snplist[isDropped],
            file=file.path(finemapQCdir,"FinemapQC_droppedSNPs.txt")
      )}
    
    if(doFlipDROPStep){ #Report message indicating how subsetting
      sink(file = logfile, append = T)
      cat("A total of",length(sums.region[[1]]$snp),"SNPs remain.\n")
      if(any(isDropped)) cat("An index of removed snps can be found in the following file: FinemapQC_droppedSNPs.txt")
      if(length(sums.region[[1]]$snp)==0){cat("Analysis stopped as no SNPs remain.\n");stop("Analysis stopped as no SNPs remain.")}
      sink()
      
    }
  }
  
  #Always return the final number of SNPs in the main report
  invisible(lapply(paste0("RMD_finemap_",traits),assignToList,value=length(sums.region[[1]]$snp),name="finalNsnp"))
  
  sink(file = logfile, append = T)
  cat("\n######\n### SuSiE finemapping results\n######\n\n")
  sink()
  
  #Run SuSiE across datasets and generate report summaries
  #A guideline for N is specified since this is highly recommended
  
  SusieFail <- logical(0) #SusieFail indicates if susie was successful across traits: TRUE indicates an error was thrown
  susie <- mapply(function(dset,trait){
    L <- opt$finemap_initialL #Define the number of credible sets to initially fit
    
    #Set a list of arguments to pass to runsusie (and internal susieR functions) with do.call
    #Is programmed this way to avoid repetition in any potential while loop for increasing the L argument
    runsusieArgs <- list(d=dset,
                         n=median(dset$N),
                         estimate_residual_variance = FALSE,
                         coverage=opt$finemap_CScoverage,
                         p=opt$priors_susie,
                         refine=opt$finemap_refine,
                         maxit=ifelse(opt$finemap_refine,10000,100),
                         check_prior=TRUE
      )
    
    tryCatch({
      
      #Run Susie
      finemap<- do.call(runsusie,c(runsusieArgs,list(L=L)))
      
      #If the fit indicates a number of CS close to the maximum number of L, recursively try a larger number
      if(opt$finemap_increaseL && !is.null(finemap$sets$cs_index)){
        while(max(finemap$sets$cs_index)>=L-2){
          L=L+10 #Iteratively increase number
          finemap<- do.call(runsusie,c(runsusieArgs,list(L=L)))
          
          if(L==100) break #Set a very high limit on the loop to avoid potential infinite loop
        }
      }
      #Alternative syntax to run susie without the coloc wrapper
      #finemap<- susie_rss(dset$beta/dset$SE,R=dset$LD,n=median(dset$N),refine=TRUE)
      
      ## Some basic diagnostics
      #susie_plot(finemap,y="PIP",add_legend=TRUE) 
      
      #Save the susie results as a resource
      finemappath<- file.path(datadir,"finemapping")
      if(!dir.exists(finemappath)){dir.create(finemappath)}
      saveRDS(finemap,file=file.path(finemappath,paste0("susie_results_",trait,".Rds")))
      
      #Assign FALSE if errors weren't thrown
      assign("SusieFail",c(SusieFail,FALSE),envir=.GlobalEnv)
      
      return(finemap)},
      
      error   = function(x){ #If there is an issue, flag in global environment and print to log file
        warning(x)
        assign("SusieFail",c(SusieFail,TRUE),envir=.GlobalEnv)
        
        sink(file = logfile, append = T)
        cat("------------------------------\n")
        cat("SuSiE fine-mapping error produced for", traits[length(SusieFail)],". Please see the following:\n")
        cat("Error in", toString(last.warning),":\n",names(last.warning),"\n")
        cat("------------------------------\n")
        sink()
      })
  },sums.region,traits,SIMPLIFY = FALSE)
  
  ######
  ### Summarise susie results across datasets
  ######
  susie_rep <- mapply(function(susie.fit,dset,trait,failed){
    
    if(failed){ #For a trait with a failed finemapping step, terminate early and return empty list for downstream logic-checks
      return(list(csFound=FALSE,sets=data.frame(snp=snplist,cs=factor(NA_character_),
                                                              variable_prob=NA_real_)))
    } else { #Otherwise, produce a summary!
      
      ## Extract initial summary elements from susie.fit
      sum_susie <- summary(susie.fit)
      csFound<- !is.null(sum_susie$cs) #Logical statement to pass through indicating TRUE if a CS has been identified
      sets<- sum_susie$vars #[sum_susie$vars$cs!=-1,] #Extract all snps
    
      ##Perform a final consistency check for the final dataset
      finalConsistencyCheck<- estimate_s_rss(method="null-mle",z=dset$beta/dset$SE, R=dset$LD, n=median(dset$N), r_tol = 1e-08)
      
      ## Format the dataset subsetting to only the list elements which are snpwise values
      basedata<- as.data.frame(dset[sapply(dset,length)==length(dset$snp)])
      
      #Combine sets object with the basedata, matched by row-index
      sets<- cbind(sets,basedata[sets$variable,c("snp","chr","pos","pvalues","beta","SE")]) #Add relevant columns from GWAS sumstats
      
      #Quick check to catch a potential circumstance where CS have not matched to snp IDs correctly
      sets$setmatch <- sapply(sets$cs,function(x)ifelse(x!=-1,paste0("L",x),NA_character_))
      for(i in seq_along(susie.fit$sets$cs)){
        setName<- names(susie.fit$sets$cs)[i]
        overlaps<- intersect(sets$snp[which(sets$setmatch==setName)],names(susie.fit$sets$cs[[i]]))
        if(length(overlaps)!=length(susie.fit$sets$cs[[i]])){ stop("Fine-mapping credible sets have not been correctly matched to SNPs for the report step. Please report this error since this will be a result of an unexpected issue in the workflow.") }
      }
      sets$setmatch <- NULL
      
      #Save PIP summaries to file
      sets_outpath<- file.path(tabdir,paste0("susie_snp_summary_",trait,".csv"))
      write.table(sets,file=sets_outpath,sep = ",",row.names=FALSE)
      tab_out<- paste0("PIP summaries for SNPs for ",trait," are tabulated in:\n", basename(sets_outpath),"\n") #Info string
      
      ### Plot the pip summaries alongside snp p-values
      
      #Mutate data for plotting
      sets <- sets %>%
        mutate(cs = if_else(cs==-1,NA_character_,paste0("L",cs)),
               thresh = if_else(!is.na(cs),as.character(snp),""))
      
      if(csFound){ #if CS are found
        #Replicate Susie plot legend labelling. Order credible sets numerically for plot visual consistency (numbering is otherwise arbitrary)
        csNames<- names(susie.fit$sets$cs)
        csNames<- csNames[order(as.numeric(gsub("L","",names(susie.fit$sets$cs))))]
        
        CSlen<- sapply(susie.fit$sets$cs,length)[csNames]
        CSmin<- susie.fit$sets$purity[csNames,"min.abs.corr"] #Here matching uses rownames
        csLabs<- paste0(names(CSlen),": C=",CSlen,"/R=",round(CSmin,3))
        
        #Explicit conversion to factor for the CS to ensure correct labelling
        sets$cs <- factor(sets$cs,levels=names(CSlen),labels=csLabs)
      }
      
      pips <- ggSummaryplot(yaxis="pip",
                            bp_range=range(sets$pos),
                            dset=sets,colourMapping=switch(csFound,sym("cs"),NULL),
                            nameColourLegend="Credible set",chr=sets$chr[1])$bpfigure +
        labs(title=ifelse(is.null(names(trait)),trait,names(trait))) #Use trait name as plot title if possible
      
      ggsave(file.path(figdir,paste0("susie_PIP_",trait,".pdf")),pips,device="pdf",units="mm",width=150,height=150)
      
      if(csFound){

        ## Visualise CS-assigned SNPs against the Obs/Exp SNP matrix.
        z_compare_cs <- kriging_rss(
          z=dset$beta/dset$SE,
          R=dset$LD,
          n=median(dset$N),
          s=finalConsistencyCheck
        )
        
        #Combine the dset 'snp' vector with the z_compare results (to ensure that the row order is consistent), then add credible sets matched by snp.
        z_compare_csplot<- cbind(dset["snp"],z_compare_cs$conditional_dist) %>%
          left_join(sets[c("snp","cs")],by="snp") %>%
          mutate(alpha=if_else(!is.na(cs),1,0)) %>%
          arrange(desc(is.na(cs)),cs) %>% #Arrange to force plotting of non-NA vals (i.e. CS), atop the non-cs snps)
          ggplot(.,aes(x=condmean,y=z,colour=cs,alpha=alpha))+
          geom_point() +
          theme_bw()+
          labs(y = "Observed z scores", x = "Expected value") +
          geom_abline(intercept = 0, slope = 1) +
          scale_colour_manual(na.value = "black", values=ggpalette,breaks=levels(sets$cs))+
          scale_alpha_continuous(range=c(0.2,1))+
          guides(color=guide_legend(title="Credible set"),alpha="none")
        
        #Save to file
        ggsave(file.path(finemapQCdir,paste0("Z_score_alignment",trait,"_withCS.pdf")),z_compare_csplot,device="pdf",units="mm",width=150,height=150)
        
        ##### Create summary file
        finemap_summary <- data.frame(credible_set_bp_range = NA,
                                      LD_Zscore_consistency=finalConsistencyCheck,
                                      CScoverage=opt$finemap_CScoverage,
                                      NSNP=NA,
                                      beta_maxSNP = paste(basedata$snp[which(abs(basedata$beta)==max(abs(basedata$beta)))],collapse=", "),
                                      p_minSNP = paste(basedata$snp[which(basedata$pvalues==min(basedata$pvalues))],collapse=", "),
                                      pip_maxSNP = NA,
                                      pip_max=NA,
                                      pip_maxSNP_nearestDownstreamGene=NA,
                                      pip_maxSNP_nearestUpstreamGene=NA,
                                      sum_susie$cs,
                                      genes_plusminus10kb_window_containing_full_cs=NA,
                                      genes_plusminus10Kb_window_overlapping_cs=NA)
        
        for(j in 1:nrow(sum_susie$cs)){
          snp_index<-as.numeric(strsplit(finemap_summary$variable[[j]], ',')[[1]])
          
          finemap_summary$NSNP[j] <- length(snp_index)
          finemap_summary$pip_max[j] <- max(susie.fit$pip[snp_index])
          
          #Extract the top PIP snp(s), then identify nearest up/downstream gene to variant (intragenic variants are tagged within their respective gene)
          topPIPsnps <- names(susie.fit$pip[snp_index])[susie.fit$pip[snp_index]==max(susie.fit$pip[snp_index])]
          topPIPnearest<- lapply(topPIPsnps,getNearestGene,positions=sets[c("snp","pos")],genesInChr=genes_checkUpOrDownstream)
          
          #Record the top PIP snps and their nearest up/downstream genes
          finemap_summary$pip_maxSNP[j]<- paste(topPIPsnps,collapse=", ")
          finemap_summary$pip_maxSNP_nearestDownstreamGene[j] <- paste(sapply(topPIPnearest,function(x)x$downstream),collapse=", ")
          finemap_summary$pip_maxSNP_nearestUpstreamGene[j] <- paste(sapply(topPIPnearest,function(x)x$upstream),collapse=", ")
          
          ss_subset<-basedata[(basedata$snp %in% names(susie.fit$pip)[snp_index]),]
          bp_range <- range(ss_subset$pos)
          
          finemap_summary$credible_set_bp_range[j]<- paste0("chr",basedata$chr[1],":",bp_range[1],"-",bp_range[2]) #span of credible set
          
          #Identify all genes around the credible set region
          Genes_subset<- getNearbyGenes(bp_range=bp_range, genesInChr = Genes)

          if(nrow(Genes_subset) > 0){
            #Write the gene names
            finemap_summary$genes_plusminus10Kb_window_overlapping_cs[j]<-paste(Genes_subset$external_gene_name, collapse=', ')
            
            #Flag any high confidence genes that contain the full credible set
            highConf <- getNearbyGenes(bp_range=bp_range, genesInChr = Genes, highConfidence = TRUE)
            if(nrow(highConf) > 0) finemap_summary$genes_plusminus10kb_window_containing_full_cs[j] <- paste(highConf$external_gene_name, collapse=', ')

            #Save details for genes nearby to credible set
            sets_outpath<- file.path(tabdir,paste0("nearby_genes_",trait,"_cs",j,".csv"))
            Genes_subset %>%
              dplyr::select(-c(start_window,end_window)) %>%
              write.table(.,file=sets_outpath,sep = ",",row.names=FALSE)
          }
        }
        finemap_summary$variable <- NULL
        
        #Write out the region analysed
        finemap_summary$region_analysed <- paste0("Chr",dset$chr[[1]],":",paste0(range(dset$pos),collapse="-"))
        
        finemap_summary <- finemap_summary %>%
          dplyr::select(region_analysed,cs,!cs) #Reorder so that credible set number comes first
        
        #Write finemapping summary to file and concatenate with previous results if file already exists
        finemapSummary<- file.path(tabdir,"results_summary_finemapping.csv")
        write.table(cbind(trait=unname(trait),finemap_summary),
                    file=finemapSummary,sep = ",",row.names=FALSE,col.names = !file.exists(finemapSummary),append = file.exists(finemapSummary))
        
        ## Generate some summary plots
        
        #Select the top snps and cs index
        topsnps<- finemap_summary[c("pip_maxSNP","cs")]
        names(topsnps)[1] <- "topsnp"
        
        #Where multiple top snps have been matched, split into long format dataset
        #some plots will allow inclusion of multiple snps per CS and downsampling is handled internally
        topsnps <- separate_rows(topsnps, all_of("topsnp"), sep = ", ")
        
        #Generate numeric only credible set labels for heatmap legend
        heatmapSets<- sets[c("snp","pos","cs")] %>%
          mutate(cs=gsub("L([0-9]+):.*","\\1",cs))
          
        ## First, plot LD in the region relative to top snps from the credible set(s) identified
        
        #Visualise the LD between Top PIP SNPs and other SNPs in the dataset / assigned to CS in heatmap
        #The plotting function expects an index of the top snps to plot ("topsnp") and the credible set to which they correspond ("cs")
        ldHeatmap<- susie_cs_ld(sets=heatmapSets,R=dset$LD,topsnps=topsnps,heatmapPalette="OrRd",discretePalette=ggpalette,
                                plotR2=TRUE,chr=sets$chr[1])
        
        #If 2+ CS, visualise correlations between different credible sets
        if(length(susie.fit$sets$cs)>1){
          cs_heatmap <- susieCScorrs(susie.fit,R=dset$LD,inc_Z=TRUE,sets=sets,topsnps =topsnps)
          ggsave(file.path(finemapQCdir,paste0("cs_correlations_",trait,".pdf")),cs_heatmap$heatmap,device="pdf",units="mm",width=150,height=150)
          
          #Save heatmaps to list which will pass to a Rmd report
          assignToList(x=paste0("RMD_finemap_",trait),value=cs_heatmap$heatmap,name="CS_corrmap")
          
        }
        
        #Assign elements relevant to having 1+ cs to report
        assignToList(x=paste0("RMD_finemap_",trait),multi=TRUE,
                     value=list(LDheatmaps=ldHeatmap,
                                zPlot_withCS=z_compare_csplot))
        
      } else {
        finemap_summary <- paste("No",coveragePcent,"credible sets could be identified")
      }
      
      #Append a series of report
      assignToList(x=paste0("RMD_finemap_",trait),multi=TRUE,
                   value=list(niter=susie.fit$niter,
                              converged=susie.fit$converged,
                              final_s_estimate=finalConsistencyCheck,
                              sets=sets,
                              nCS=length(susie.fit$sets$cs),
                              coloc_format_data=dset,
                              finemapSummary=finemap_summary
                              ))
      
      #Sink directly to file
      sink(file = logfile, append = T)
      cat("------------------------------\n")
      cat("SuSiE finemap result for ", trait,":\n",sep="")
      cat("Model fitted using", susie.fit$niter,"iterations, and converged =",susie.fit$converged,"\n\n")
      print(finemap_summary)
      cat("\n",tab_out)
      cat("------------------------------\n")
      sink()
      
      #Return results just in case
      return(list(finemap_summary=finemap_summary,tab_out=tab_out,csFound=csFound, sets=sets))
    }
  },susie,dset=sums.region,trait=traits,failed=SusieFail,SIMPLIFY = FALSE)
  
  
  #Identify the objects storing files to write to a report
  finemapReports<- grep("RMD_finemap_",ls(envir=.GlobalEnv),value=TRUE)
  
  #Render the QC report per-trait
  invisible(mapply(renderReport,
         params=finemapReports,
         outfile=file.path(normalizePath(reportdir),paste0("02_",gsub("RMD_","",finemapReports),"_report.html")),
         MoreArgs = list(template=file.path(opt$scriptsDir,"rmd","finemap_report.Rmd"))))
  
  #Optionally prepare for global report; since Finemapping reports are handled per-trait, combine them into a single list to recursively generate reports
  if(!is.null(opt$rdsOut)){  
    P02<- lapply(finemapReports,get)
    names(P02) <- lapply(P02,function(x) unlist(x$trait))
    saveRDS(P02,file.path(opt$rdsOut,"P02.Rds"))
  }
  
  if(file.exists(file.path(tabdir,"results_summary_finemapping.csv"))){
    sink(file = logfile, append = T)
    cat("\nSummaries of credible sets identified by susie across traits are all available in the file: results_summary_finemapping.csv\n")
    sink()
  }
  
  
  #Flag either of the SuSiE calls failed
  if(any(SusieFail)){
    
    ######
    ### Check alignment of Beta and LD, may lead to issues with LD matrix convergence if not aligned in the same direction
    ######
    pdf(file=file.path(finemapQCdir,"Beta_LD_alignment.pdf"),width=7,height=7)
    lapply(sums.region,check_alignment)
    dev.off()
    
    sink(file = logfile, append = T)
    cat("The SuSiE finemapping step returned an error for at least one trait (see above).\n")
    cat("Please check the resources returned in the finemapQC directory to evaluate the  alignment of summary statistic test statistics against those expected given the LD matrix.\nSee also https://chr1swallace.github.io/coloc/articles/a02_data.html for details.\n")
    sink()
    
  } else if(!all(sapply(susie_rep[1:ifelse(length(susie_rep)>2,2,length(susie_rep))],function(x)x$csFound)) && opt$runMode!="finemapOnly"){ #Return this message on the basis of only the first two traits
    sink(file = logfile, append = T)
    cat(coveragePcent," credible sets were not identified for at least one of ",paste0(traits[1:ifelse(length(traits)>2,2,length(traits))],collapse=" & "),". Therefore, coloc.susie was not used.\nAdjusting the susie_rss coverage parameter using the --finemap_CScoverage option may allow lower coverage credible sets to be identified but these should be treated with caution (See: https://chr1swallace.github.io/coloc/articles/a06_SuSiE.html).\n",sep='')
    sink()
  }
  
  
  #Extract PIPs for each trait 
  snp_PIP <- mapply(function(x,trait,keepcols){cbind(x$sets[,c("snp","variable_prob")],trait)},
                    susie_rep, traits, SIMPLIFY = FALSE) %>%
    do.call(rbind,.) 
  
  #Combine SuSiE PIP results with the minimal dataset
  minimal_df <- right_join(minimal_df,snp_PIP,by=c("snp","trait"))  #add each PIP; but retain only snps passing QC
  
  #If any CS were identified, extract these across traits  and add to the minimal DF
  if(any(sapply(susie_rep,function(x)x$csFound))){
    
    #Extract credible set summaries per-snp across traits.
    #Subset to columns, adjust credible set labels, then reduce list across matched DFs and unite cols
    snp_CS<- mapply(function(x,trait){ 
      y <- x$sets[,c("snp","cs")]
      if(length(levels(y$cs))>0){levels(y$cs) <- paste0(trait,":",gsub("^L([0-9]+)\\:.*","\\1",levels(y$cs)))}
      return(y)
    },susie_rep,traits,SIMPLIFY = FALSE) %>%
      reduce(full_join,by="snp") %>%
      tidyr::unite(.,col=cs, starts_with("cs"), sep = " & ", remove = TRUE, na.rm = TRUE) %>%
      mutate(cs=as.factor(if_else(cs=="",NA_character_,cs)))
    
    #Combine SuSiE credible sets with the minimal dataset
    minimal_df <- right_join(minimal_df,snp_CS,by="snp")
  }
  
  minDFtext<- " with SNPwise finemapping information "
  
} else {
  #If susie is skipped entirely, produce some dummy values required for subsequent logic checks
  susie_rep <- lapply(seq_along(sums.region),function(x)list(csFound=FALSE))
    
  SusieFail <- sapply(seq_along(sums.region),function(x) FALSE)
  
  minDFtext<- " "
  
}

## Save the minimal dataset to file, for later plotting
sink(file = logfile, append = T)
cat("---------\nSaving file containing harmonised summary statistics",minDFtext,"across traits to the file: data/datasets/harmonised_sumstats.csv\n",sep="")
sink()

write.table(minimal_df,file=file.path(datadir,"datasets","harmonised_sumstats.csv"),sep = ",",row.names=FALSE,col.names = TRUE,quote = FALSE)

if(opt$runMode=="finemapOnly"){
  sink(file = logfile, append = T)
  cat("Skipping colocalisation analysis because the 'finemapOnly' setting is declared in the --runMode option.\n")
  sink()
  
} else if(length(traits)==1){
    sink(file = logfile, append = T)
    cat("Skipping colocalisation analysis because only one trait has been provided.\n")
    sink()
    
} else {
  P03 <- list() #Prepare list in which to store coloc report parameters
  
  #Subset to traits 1 and 2, if more have been supplied for prior steps.
  if(length(traits)>2){
    sink(file = logfile, append = T)
    cat("Colocalisation analysis can only be performed for pairs of traits. Proceeding with the first two traits only.\n")
    sink()
    
    traits<- traits[1:2]
    SusieFail<- SusieFail[1:2]
  }
  
  #Separate into distinct objects consistent with initial script config for downstream analysis.
  #Ensure only the first two datasets are passed onward
  sums1.region <- sums.region[[1]]
  sums2.region <- sums.region[[2]]

######
### Proceed with colocalisation step
######
sink(file = logfile, append = T)
cat("\n######\n### Colocalisation results\n######\n\n")
sink()

## Extract priors for coloc steps [coloc.abf and coloc.susie defaults unless options are otherwise modified]
colocPriors<- strsplit(c(opt$priors_coloc.abf,
                         opt$priors_coloc.susie),split=",")
names(colocPriors) <- c("coloc.abf","coloc.susie") #set names respective to the relevant function

P03$colocPriors <- colocPriors <- lapply(colocPriors,function(x){
  x <- as.numeric(x)              #Set to numeric
  names(x) <- c("p1","p2","p12")  #Declare relevant argument name
  x<- as.list(x)                  #Convert to list
  return(x)})

#If subsetting to named gene sources, generate the vector, and cat message to file if no rows would remain
if(!is.null(opt$restrict_nearby_gene_plotting_source)){geneSources <-strsplit(opt$restrict_nearby_gene_plotting_source,split=",")[[1]] #String split on comma to identify sources
  } else { geneSources <- NULL }

#If either SuSiE call failed, a credible set has not been found for both traits, or if analysis is passed direct to coloc, run coloc.abf
if(any(SusieFail) || !all(sapply(susie_rep[1:ifelse(length(traits)>2,2,length(traits))],function(x)x$csFound)) || opt$runMode %in% c("doBoth","skipSusie")){  
  
  #Run coloc.abf via do.call to allow passing of the priors list
  clc.abf<- do.call(coloc.abf,c(colocPriors$coloc.abf,list(dataset1=sums1.region, dataset2=sums2.region)))
  
  P03$coloc_abf_summary <- clc.abf$summary #Add coloc.abf summary to report
  
  #Save the coloc abf results as an Rds resource
  colocpath<- file.path(datadir,"colocalisation")
  if(!dir.exists(colocpath)){dir.create(colocpath)}
  saveRDS(clc.abf,file=file.path(colocpath,paste0("coloc_abf_",paste0(traits,collapse="_"),".Rds")))

  sink(file = logfile, append = T)
  cat("Colocalisation will now be performed without passing first to SuSiE [see coloc::coloc.abf].\n")
  cat("Thus, the single causal variant assumption has not been relaxed.\n")
  cat("coloc.abf was performed with the following priors:\n")
  print(clc.abf$priors)
  cat("coloc.abf results summary:\n")
  print(clc.abf$summary)
  cat("This summary is also returned in the file: results_summary_coloc_abf.csv\n")
  sink()
  
  ## Save the results summary in a table which can be readily combined with other outputs
  Sum_abf<- t(enframe(c(traits,reg_range,region=paste0("Chr", reg_range["chr"],":",reg_range["start"],"-",reg_range["stop"]),clc.abf$priors,clc.abf$summary)))
  Sum_abf[1,1:length(traits)] <- paste0("TraitID_",1:length(traits))
  
  write.table(Sum_abf,file=file.path(tabdir,"results_summary_coloc_abf.csv"),sep = ",",row.names=FALSE,col.names = FALSE,quote = FALSE)
  
  
  if(clc.abf$summary[["PP.H4.abf"]]>0.2){
    ## Combine with positional data, and label the 95% credible SNPs, saving per-SNP statistics to file
    snpwiseSummary <- sums1.region %>%
      as.data.frame() %>%
      dplyr::select(snp,chr,pos) %>%
      right_join(clc.abf$results,by="snp") %>%
      arrange(desc(SNP.PP.H4)) %>%
      mutate(cs_coloc=cumsum(SNP.PP.H4),                                         #Establish colocalisation credible sets
             cs_coloc=case_when(row_number() == 1 ~ "within set",                #Top snp is always within set
                                lag(cs_coloc) < 0.95 ~ "within set",             #lag checks previous record. If this value is <0.95 then SNP is in set
                                TRUE ~ "outside set"))
    
    write.table(snpwiseSummary,file=file.path(tabdir,"coloc.abf_snpwise_PP_H4_abf.csv"),sep = ",",row.names=FALSE)
    
    ### Prep data for the SNP pp plot
    plt <- snpwiseSummary %>%
      rename(SNP.PP=SNP.PP.H4) %>%
      mutate(cs_coloc= factor(cs_coloc,levels=c("within set","outside set")),
             lab = if_else(row_number()<11 & SNP.PP > 0.1 & cs_coloc=="within set",as.character(snp),""))
    
    ## Generate the plot and save to report
    P03$abf.PP.plot <- snpPP <- ggSummaryplot(dset=plt,
                                              yaxis="SNP.PP",
                                              bp_range=range(plt$pos),
                                              chr=reg_range[["chr"]],
                                              shapeMapping="cs_coloc",
                                              nameShapeLegend="Colocalisation\n(95% credible SNPs)")$bpfigure +
      geom_text_repel(aes(label=lab), max.overlaps = 20, na.rm=TRUE, show.legend = FALSE)
    
    
    
    ## Stack the plot above nearby genes
    P03$abf.PP.plot.geneNear <- snpPPstackGenes <- plotGenesStack(snpPP,
                                                                  genesInChr=Genes,
                                                                  bp_range=range(snpPP$data$pos[snpPP$data$cs_coloc=="within set"]),
                                                                  geneSources=geneSources,
                                                                  gene_tracks=opt$gene_tracks,
                                                                  nudge_y=0.3)
    
    #Reformat the Y axis of the with-gene plot to ensure no label clashes
    P03$abf.PP.plot.geneNear$patches$plots[[1]] <- P03$abf.PP.plot.geneNear$patches$plots[[1]]+labs(y="Posterior probability\nfor being shared variant")
      
    
    #Save the snp-wise posterior probability plot, and genes nearby to CS if any were found.
    prefix <- paste0(c("coloc_abf",traits,"SNPposteriorProbs"),collapse="_")
    file <- ggsave(file.path(figdir,paste0(prefix,".pdf")),snpPP,device="pdf",units="mm",width=150,height=175)
    if(!"character" %in% class(snpPPstackGenes)){file <- c(file,ggsave(file.path(figdir,paste0(prefix,"_nearbygenes.pdf")),snpPPstackGenes,device="pdf",units="mm",width=150,height=175))}
    
    #Identify all genes nearby to the credible snps range for each row
    nearbyGenes<- getNearbyGenes(bp_range=range(snpPP$data$pos[snpPP$data$cs_coloc=="within set"]),
                                 genesInChr=Genes)
    
    #If any genes have been matched, then write to a summary table
    if(nrow(nearbyGenes)>0){
      nearbyGenes %>%
        dplyr::select(-c(start_window,end_window)) %>%
        write.table(.,file=file.path(tabdir,paste0("nearby_genes_coloc_abf.tsv")),sep = "\t",row.names=FALSE)
    }
    
    sink(file = logfile, append = T)
    cat(sep='',"\nPer-SNP posterior probabilities for being a shared variant (assuming H4=TRUE) are tabulated in: coloc.abf_snpwise_PP_H4_abf.csv\nSee also: ",basename(file[1]),"\n")
    if(nrow(nearbyGenes)>0) cat(sep='',"\nGenes within a 10Kb window around the 95% credible SNPs from the analysis are tabulated in: nearby_genes_coloc_abf.tsv\n",ifelse(length(file)==2,paste0("See also: ",basename(file[2])),""),"\n")
    sink()
  } else {
    sink(file = logfile, append = T)
    cat("\nFurther summaries have not been returned because negligible support was found for the shared variant hypothesis (H4).\n")
    sink() 
  }
}

#If both susie calls are successful and credible sets identified, run coloc.susie.
#Note that plots will overwrite any from Coloc.abf.
if(!any(SusieFail) && all(sapply(susie_rep[1:2],function(x)x$csFound)) && opt$runMode %in% c("doBoth","trySusie")){
  
  ###run coloc.susie based on susie outputs and priors set
  clc<- do.call(coloc.susie,c(colocPriors$coloc.susie,list(dataset1=susie[[1]], dataset2=susie[[2]])))
  
  P03$coloc_susie_summary <- clc$summary #Add coloc.susie summary to report
  
  #Save the coloc susie results as an Rds resource
  colocpath<- file.path(datadir,"colocalisation")
  if(!dir.exists(colocpath)){dir.create(colocpath)}
  saveRDS(clc,file=file.path(colocpath,paste0("coloc_susie_",paste0(traits,collapse="_"),".Rds")))

  sink(file = logfile, append = T)
  cat("------------------------------\n")
  cat("coloc.susie was performed with the following priors:\n")
  print(clc$priors)
  invisible(print(clc$summary)) #Call summary invisibly first to avoid triggering bug where no output is printed
  cat("coloc.susie results summary:\n")
  print(clc$summary)
  cat("This summary is also returned in the file: results_summary_coloc_susie.csv.\n")
  sink()
  
  ## Save the results summary in a table which can be readily combined with other outputs
  configSummary <- c(traits,reg_range,region=paste0("Chr", reg_range["chr"],":",reg_range["start"],"-",reg_range["stop"]),clc$priors)
  if(nrow(clc$summary)>1){ Sum_clcsusie<- t(configSummary) %>% .[rep(1,nrow(clc$summary)),] %>% cbind(.,clc$summary)
  } else { Sum_clcsusie<- t(c(configSummary,clc$summary)) }
  colnames(Sum_clcsusie)[1:length(traits)] <- paste0("TraitID_",1:length(traits))
  write.table(Sum_clcsusie,file=file.path(tabdir,"results_summary_coloc_susie.csv"),sep = ",",row.names=FALSE,col.names = TRUE,quote = FALSE)
  
  possColoc<- clc$summary$PP.H4.abf>0.2
  if(any(possColoc)){
    
    #Identify 95% credible SNPs for each comparison row, splitting each row results into a list and then recombining into a single data frame
    annot <- list()
    for(i in 2:ncol(clc$results)){
      annot[[length(annot)+1]]  <- clc$results[,c(1,i),with=FALSE] %>%
        arrange(desc(.[[2]])) %>%
        mutate(cs_coloc=cumsum(.[[2]]),                                         #Establish colocalisation credible sets
               cs_coloc=case_when(row_number() == 1 ~ "within set",             #Top snp is always within set
                                  lag(cs_coloc) < 0.95 ~ "within set",          #lag checks previous record. If this value is <0.95 then SNP is in set
                                  TRUE ~ "outside set")) %>%
        rename_with(.cols="cs_coloc",function(x)paste0("cs_coloc.row",i-1))
    }
    annot <- Reduce(function(x,y)full_join(x,y,by="snp"),annot)
    
    #Add genomic position and susie credible set annotations
    snpwiseSummary<- sums1.region %>%
      as.data.frame() %>%
      dplyr::select(snp,chr,pos) %>%
      full_join(snp_CS,by="snp") %>%
      rename(cs_susie=cs) %>%
      right_join(annot,by="snp")
    
    ## Save PP.H4 summaries, with additional details, to file
    write.table(snpwiseSummary,file=file.path(tabdir,"coloc.susie_snpwise_PP_H4_abf.csv"),sep = ",",row.names=FALSE)
    
    sink(file = logfile, append = T)
    cat("\nPer-SNP posterior probabilities for being a shared variant (assuming H4=TRUE) from each pairwise comparison of fine-mapping credible sets across traits are tabulated in: coloc.susie_snpwise_PP_H4_abf.csv\n")
    cat("\nGenerating additional summaries for pairwise comparisons with PP.H4>0.2 (a low threshold to ensure outputs are produced for every reasonable PP.H4).\n\n")
    sink()
    
    ## Generate the snp PP plot for each row of coloc.susie results that are suggestive of colocalisation
    possColocRows<-paste0("row",which(possColoc))
    names(possColocRows) <- possColocRows
    snp_PPplots <- lapply(possColocRows,function(rows,snpwiseSummary){
      
      cols <- c("snp","pos","cs_susie",grep(paste0("(",rows,"$|abf)"),colnames(snpwiseSummary),value = TRUE))
      
      #Prepare data for plotting
      plt <- snpwiseSummary[cols] %>% 
        rename(SNP.PP=matches("SNP.PP.H4"),
               cs_coloc=matches("cs_coloc")) %>%
        arrange(desc(SNP.PP)) %>%
        mutate(cs_coloc= factor(cs_coloc,levels=c("within set","outside set")),
               lab = if_else(row_number()<11 & SNP.PP > 0.1 & cs_coloc=="within set",as.character(snp),""))
      
      ## Generate the plot
      snpPP <- ggSummaryplot(dset=plt,
                             yaxis="SNP.PP",
                             bp_range=range(plt$pos),
                             chr=plt$chr[1],
                             colourMapping="cs_susie",
                             nameColourLegend=paste0("Fine-mapping\n(Trait: ",coveragePcent," credible set)"),
                             shapeMapping="cs_coloc",
                             nameShapeLegend="Colocalisation\n(95% credible SNPs)")$bpfigure +
        geom_text_repel(aes(label=lab), max.overlaps = 20, na.rm=TRUE, show.legend = FALSE)
      
      return(snpPP)
    },snpwiseSummary)
    
    #Identify credible SNP range for each snpPP plot and then mapply across the original plots and the range to visualise nearby genes
    credibleSnpRange<-lapply(snp_PPplots,function(snpPP)range(snpPP$data$pos[snpPP$data$cs_coloc=="within set"]))
    
    ## Stack each plot above nearby genes
    snpPPstackGenes <- mapply(plotGenesStack,snp_PPplots,credibleSnpRange,
                              SIMPLIFY = FALSE, MoreArgs = list(
                                genesInChr=Genes,
                                geneSources=geneSources,
                                gene_tracks=opt$gene_tracks,
                                nudge_y=0.3)
    )
    
    #Reformat the Y axis of the with-gene plot to ensure no label clashes
    snpPPstackGenes <- lapply(snpPPstackGenes,function(x){
      x$patches$plots[[1]] <- x$patches$plots[[1]]+labs(y="Posterior probability\nfor being shared variant")
      return(x)
    })
    
    #Assign the per-SNP plots and the plots with nearby genes to report list
    P03$susie.PP.plots <- list()
    for(i in seq_along(snp_PPplots)){  
      P03$susie.PP.plots[[length(P03$susie.PP.plots)+1]] <- list(ALLsnps=snp_PPplots[[i]],snpsWgenes=snpPPstackGenes[[i]])
    }
    names(P03$susie.PP.plots) <- names(snp_PPplots)
    
    
    ### Save plots to file
    
    #Generate a prefix for plots of each coloc.susie comparison row
    rowCompares<- paste0(traits[1],"cs",clc$summary$idx1[possColoc],"_",traits[2],"cs",clc$summary$idx2[possColoc])
    prefix <- sapply(rowCompares,function(x)paste0(c("coloc_susie",x,"SNPposteriorProbs"),collapse="_"))
    
    invisible(mapply(function(snp_PPplots,snpPPstackGenes,prefix){
      ggsave(file.path(figdir,paste0(prefix,".pdf")),snp_PPplots,device="pdf",units="mm",width=150,height=100)
      if(!"character" %in% class(snpPPstackGenes)){ggsave(file.path(figdir,paste0(prefix,"_nearbygenes.pdf")),snpPPstackGenes,device="pdf",units="mm",width=150,height=175) }
    },snp_PPplots,snpPPstackGenes,prefix))
    
    sink(file = logfile, append = T)
    cat("Per-SNP posterior probabilities for being a shared variant are plotted in .pdf files, named with the format:",paste0(basename(gsub("[0-9]+","<number>",prefix[1])),collapse="\n"))
    sink()
    
    ## Write out all nearby genes
    #Identify all genes nearby to the credible snps range for each row
    nearbyGenes<- lapply(credibleSnpRange,getNearbyGenes,genesInChr=Genes)
    
    #If any genes have been matched, then write to a summary table
    genesMatched<- sapply(nearbyGenes,function(x)nrow(x)>0)
    if(any(genesMatched)){
      
      #Annotate nearby genes lists with comparison rows, and concatenate into a single table
      allnearbygenes<- mapply(cbind,susie_CS_comparisons=rowCompares[genesMatched],nearbyGenes[genesMatched],SIMPLIFY = FALSE) %>%
        Reduce(rbind,.)
      
      #Check for any top SNPs tagged tagged multiple times, if present concatenate their credible set assignment column and deduplicate the dataset
      allnearbygenes$checkString <- with(allnearbygenes,paste0(external_gene_name,chromosome_name,start_position,end_position))
      tab<-table(allnearbygenes$checkString)
      if(any(tab>1)){
        dups <- names(tab)[tab>1]
        for(i in seq_along(dups)){ 
          duprows<- allnearbygenes$checkString==dups[i]
          allnearbygenes$susie_CS_comparisons[duprows] <- paste0(allnearbygenes$susie_CS_comparisons[duprows],collapse = ",") 
        }
        allnearbygenes <- unique.data.frame(allnearbygenes)
      }
      allnearbygenes$checkString <- NULL
      
      allnearbygenes %>%
        dplyr::select(-c(start_window,end_window)) %>%
        write.table(.,file=file.path(tabdir,paste0("nearby_genes_coloc_susie.tsv")),sep = "\t",row.names=FALSE)
      
      sink(file = logfile, append = T)
      cat("\nGenes within a 10Kb window around the 95% credible SNPs from pairwise comparisons with PP.H4>0.2 are tabulated in: nearby_genes_coloc_susie.tsv\n")
      sink()
    }
  } else {
    sink(file = logfile, append = T)
    cat("\nFurther summaries have not been returned because negligible support was found for the shared variant hypothesis (H4).\n")
    sink() 
  }
}

######
### Generate the colocalisation report
######
if(!is.null(opt$rdsOut)) saveRDS(P03,file.path(opt$rdsOut,"P03.Rds")) #Optionally save P03 to an Rds file for global report generation

report03<- file.path(normalizePath(reportdir),paste0("03_colocalisation_report.html"))
renderReport(params="P03",
             outfile=report03,
             template=file.path(opt$scriptsDir,"rmd","coloc_report.Rmd"))

    
} #Bracket indicating the end of the else condition performed when opt$runMode!="finemapOnly"

sink(file = logfile, append = T)
cat("------------------------------\n")
cat("An html report overviewing the analyses performed can be found in the main results directory.\n")
sink()
