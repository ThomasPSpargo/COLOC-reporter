#####
# Author: Thomas Spargo (thomas.spargo@kcl.ac.uk)
# Obtained from GitHub repository: https://github.com/ThomasPSpargo/COLOC-reporter
#
# This script returns figures which compare LD between credible sets identified across multiple traits
#####


## Setup

suppressPackageStartupMessages(library(optparse))

option_list = list(

make_option("--susieFitsDir", action="store", default=NULL, type='character',
            help="Path to directory containing Rds files storing fine-mapping results. Note that at least 1 credible set is expected across traits."),
make_option("--LDmatrix", action="store", default=NULL, type='character',
            help="Path and prefix for ld matrix files containing all SNPs used within SuSiE fine-mapping analysis; expects files ending with the suffixes: '.snplist','.ld' [optional; required if --LDreference is not provided]"),
make_option("--plink", action="store", default='plink', type='character',
            help="Path to PLINK executable (syntax written for PLINK 1.9). By default, has the value 'plink' [only required if using --LDreference option"),
make_option("--LDreference", action="store", default=NA, type='character',
            help="Path to, and prefix for, per-chromosome PLINK binary files used to compute LD matrix for SNPs analysed. Expected format of <prefix>i[.bim/.bed/.fam], where i is the chromosome number [optional; required if --LDmatrix is not provided]"),
make_option("--harmonisedSumstats", action="store", default=NULL, type='character',
            help="A file containing harmonised summary statistic information relating to the files from susieFitsDir"),
make_option("--helperFunsDir", action="store", default=NULL, type='character',
            help="Filepath to directory containing helper functions used within the script"),
make_option("--stackHeatmaps", action="store", default=FALSE, type='logical',
            help="logical, defaulting to FALSE. Set TRUE indicate that heatmap figures for the entire region of SNPs analysed and for the snp credible sets only should be stacked vertically into a single figure."),
make_option("--rdsOut",action="store", default=NULL, type='character',
            help="Specify a file path and prefix for an Rds file containing figures generated within this script '.Rds' will be appended [optional]"),
make_option("--rdsOnly", action="store", default=FALSE, type='logical',
            help="Logical, defaulting to FALSE, if set to TRUE, plots will only output as .Rds files; see the --rdsOut option."),
make_option("--outdir", action="store", default=NULL, type='character',
            help="Filepath to directory within which figure files will be written.")
)
opt = parse_args(OptionParser(option_list=option_list))

suppressPackageStartupMessages({
  library(susieR)
  library(tidyverse)
  library(data.table)
  library(patchwork)
})

if(is.null(opt$outdir) && is.null(opt$rdsOut)) { cat("ERROR: Please specify arguments to either or both of --rdsOut and --outdir. Terminating script\n"); q("no") }
if(opt$rdsOnly && is.null(opt$rdsOut)) warning("Rds files will not be saved the --rdsOut option has not been set.")

test <- FALSE # test <- TRUE            
if(test){
  opt <- list()
  opt$susieFitsDir <- "/Users/tom/OneDrive - King's College London/PhD/PhD project/COLOC/git.local.COLOC-reporter/testing/tidy_processingTEST_PD.SZ.chr17_coloc/data/finemapping"
  opt$LDmatrix <- "/Users/tom/OneDrive - King's College London/PhD/PhD project/COLOC/git.local.COLOC-reporter/testing/tidy_processingTEST_PD.SZ.chr17_coloc/data/LDmatrix/ld_matrix"
  opt$harmonisedSumstats <- "/Users/tom/OneDrive - King's College London/PhD/PhD project/COLOC/git.local.COLOC-reporter/testing/tidy_processingTEST_PD.SZ.chr17_coloc/data/datasets/harmonised_sumstats.csv"
  
  opt$helperFunsDir <- "~/OneDrive - King's College London/PhD/PhD project/COLOC/git.local.COLOC-reporter/scripts/helper_functions"
  opt$stackHeatmaps <- FALSE
  opt$rdsOut <- "~/Downloads/testfile"
  opt$outdir <- "~/Downloads"
}

#Setup a colour palette to be used across ggplots (colours taken from plot SuSiE)
ggpalette = c("dodgerblue2", "green4", "#6A3D9A", "#FF7F00", 
              "gold1", "skyblue2", "#FB9A99", "palegreen2", "#CAB2D6", 
              "#FDBF6F", "gray70", "khaki2", "maroon", "orchid1", "deeppink1", 
              "blue1", "steelblue4", "darkturquoise", "green1", "yellow4", 
              "yellow3", "darkorange4", "brown")

#Read-in custom helper functions stored in directory specified by opt$scriptsDir
list.files(opt$helperFunsDir,full.names = TRUE,pattern=".R") %>%
  lapply(.,source) %>%
  invisible(.)

### RUN

#Read in the susie results for each trait and Identify trait IDs from Rds filenames
fits<- list.files(opt$susieFitsDir,full.names=TRUE,pattern = ".Rds")
susie<- lapply(fits,readRDS)
traits <- gsub(".*_(.*).Rds$","\\1",fits)

#Check that there is at least 1 credible set present
csFound<- any(sapply(susie,function(x)!is.null(summary(x)$cs))) #Logical check whether any CS are present in the data
if(!csFound) { cat("Terminating script since no fine-mapping credible sets have been identified across the SuSiE models supplied.\n"); q("no") }

## Read in the harmonised sumstats file, which gives positional and credible set information
sets <- tibble(fread(opt$harmonisedSumstats))
heatmapSets<- unique.data.frame(sets[c("snp","cs","pos")]) #Extract the positional information

######
### Read in or generate the LD matrix
######

suffix <- c(".snplist",".ld")
if(!is.null(opt$LDmatrix)){
  LDpath <- paste0(opt$LDmatrix,suffix)
} else {
  LDpath <- paste0("ld_matrix",suffix)
}

if(any(!file.exists(LDpath))){
  cat("Computing LD matrix...")
  
  #Write out sampled SNPs and generate LD matrix
  write(heatmapSets$snp, file="temp_snplist_forLDmatrix.txt")
  
  #Syntax based on plink v 1.9
  system(paste0(opt$plink,' --bfile ', opt$LDreference,sets$chr[1],
                ' --extract temp_snplist_forLDmatrix.txt',
                ' --r square',
                ' --write-snplist',
                ' --keep-allele-order',
                ' --out ld_matrix'))
  
  cat("Done\n")
  
  #Logical to indicate to remove any LD files if they were newly generated
  rmLDfiles <- TRUE
} else {
  rmLDfiles <- FALSE
  cat("Existing LD matrix found, skipping call to PLINK.\n")
}

#Read in the ld matrix and snp names as returned by plink 
ld <- as.matrix(fread(LDpath[2]))
ld_names<- scan(LDpath[1],what=character())
dimnames(ld)<-list(ld_names, ld_names) #assign dimnames to LD object

#Remove any LD files if newly generated
if(rmLDfiles) system("rm ld_matrix.ld ld_matrix.log ld_matrix.snplist temp_snplist_forLDmatrix.txt")


## Ensure that, minimally, all credible set SNPs are represented in the LD matrix
#Extract Ids for all snps assigned to credible sets
csSNPs<- lapply(susie,function(x)lapply(x$sets$cs,names))
csSNPs<- unique(unlist(csSNPs))

#Perform check
if(!all(csSNPs %in% ld_names)){
  warning("Not all SNPs assigned to credible sets were identified in the LD matrix snplist")
  if(!any(csSNPs %in% ld_names)) { cat("Terminating script because no snps assigned to credible sets were identified in the LD matrix snplist.\n"); q("no") }
}

## Filter to only SNPs for which positional information are provided
if(!all(ld_names %in% heatmapSets$snp)){
  if(!any(ld_names %in% heatmapSets$snp)){ cat("Terminating script because no snps were provided positional information within the file supplied to --harmonisedSumstats.\nIf this is unknown, a sequence of relative positions could be supplied.\n"); q("no") }else{
    cat("Filtering down LD matrix to include only SNPs included within the file supplied to --harmonisedSumstats.\n")
    ld <- ld[dimnames(ld)[[1]] %in% heatmapSets$snp,dimnames(ld)[[2]] %in% heatmapSets$snp]      
  }
}

#Generate the LD heatmaps
ldHeatmap<- susie_cs_ld(susie.fit=susie ,R=ld,sets=heatmapSets,discretePalette=ggpalette,plotR2=TRUE,trait=traits,
                        chr=switch("chr" %in% names(sets),sets$chr[1],NULL), #Include chromosome number if provided in the sets information
                        separateLegends=opt$stackHeatmaps #Do NOT collect legends when heatmaps are to be stacked; this will be done at the final stacking stage
                        )

#Write plots to Rds file
if(!is.null(opt$rdsOut)){ 
  if(!dir.exists(dirname(opt$rdsOut))) dir.create(dirname(opt$rdsOut),recursive = TRUE)
  saveRDS(ldHeatmap,paste0(opt$rdsOut,".Rds"))
  
  if(opt$rdsOnly) q(save="no")
}

if(opt$stackHeatmaps){

  #Collect these into a single 'tall figure
  heatmapStack<- wrap_plots(ldHeatmap[[1]],ncol=1)/wrap_plots(ldHeatmap[[2]],ncol=1)+
    plot_layout(nrow=2,guides="collect",heights = c("1","1"))
  
  ggsave(file.path(opt$outdir,paste0("LD_CScorrs_stacked_",paste0(traits,collapse="_"),".pdf")),
                   device="pdf",units="mm",width=200,height=120)
} else {

  ggsave(file.path(opt$outdir,paste0("LD_CScorrs_allSNPs_",paste0(traits,collapse="_"),".pdf")),
         ldHeatmap$heatmap_allSNPs,device="pdf",units="mm",width=200,height=120)
  
  ggsave(file.path(opt$outdir,paste0("LD_CScorrs_csSNPs_",paste0(traits,collapse="_"),".pdf")),
         ldHeatmap$heatmap_CSonly,device="pdf",units="mm",width=200,height=120)
}

