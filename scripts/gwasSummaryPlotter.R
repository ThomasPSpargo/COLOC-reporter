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
  make_option("--harmonisedSumstats", action="store", default=NULL, type='character',
              help="A file containing harmonised summary statistics across two or more traits"),
  make_option("--traits", action="store", default=NULL, type='character',
              help="Comma separated list of two trait ID for traits to analyse (e.g. 'P1,P2'), as relevant to the --harmonisedSumstats input."),
  make_option("--GWASconfig", action="store", default=NULL, type='character',
              help="[optional] Path to file specifying configuration of GWAS sumstats. This is used only for identifying trait labels associated with traits in the --traits option."),
  make_option("--GWASsumplots", action="store", default="p", type='character',
              help="Comma separated string declaring variables to plot; defaults to 'p'. The preset plot formatting is compatable with any combination of 'PIP','p','z','abs_z','beta'. Other strings which name numeric columns in the supplied data can also be plotted. PIP refers to posterior inclusion probabilities for snps as indicated from SuSiE in a column called 'variable_prob'; 'p' expects the column 'pvalues' which will be converted into -log10(p); 'z' or 'abs_z' requires either the column 'z' (case insensitive)  or 'beta' and 'SE' (from which Z-scores can be determined) and from which absolute Z-scores are determined if supplying 'abs_z' ; 'beta' requires the column 'beta'"),
  make_option("--GWASsumplots_onefile", action="store", default=FALSE, type='logical',
              help="Logical, defaults to FALSE. When --GWASsumplots has more then one element, indicate TRUE to return all comparisons bound into a single figure; plots will be stacked vertically, in the order of the --GWASsumplots vector."),
  make_option("--GWASsumplots_incfinemapping", action="store", default=TRUE, type='logical',
              help="Logical, defaults to TRUE. If TRUE, and if any credible sets are identified across the two traits, plots returned by --GWASsumplots will include colouring indicating (if any) the finemapping credible sets assigned to each snp."),
  make_option("--genomeAlignment", action="store", default=37, type='numeric',
              help="Indicate the reference genome to which summary statistics and LD reference are aligned. Defaults to 37."),
  make_option("--helperFunsDir", action="store", default=NULL, type='character',
              help="Filepath to directory containing helper functions used within the script"),
  make_option("--rdsOut",action="store", default=NULL, type='character',
              help="Specify a directory in which to return figures in an Rds file [optional]"),
  make_option("--outdir", action="store", default=NULL, type='character',
              help="Filepath to directory to which files will be written")
)

opt = parse_args(OptionParser(option_list=option_list))

#######
### Begin script
#######

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)    #For reading-in Gwas config
  library(patchwork)     #arranging summary plot
})


test <- FALSE # test <- TRUE
if(test==TRUE){
  #####
  setwd("/Users/tom/OneDrive - King's College London/PhD/PhD project/COLOC/git.local.COLOC-reporter/testing")
  opt <- list()
  
  opt$harmonisedSumstats <- "/Users/tom/OneDrive - King's College London/PhD/PhD project/COLOC/git.local.COLOC-reporter/testing/tidy_processingTEST_PD.SZ.chr17_coloc/data/datasets/harmonised_sumstats.csv"
  opt$traits <- "PD,SZ"
  
  opt$GWASconfig <- "./GWAS_samples_testing.txt"
  opt$genomeAlignment <- 37
  
  opt$helperFunsDir <- "/Users/tom/OneDrive - King's College London/PhD/PhD project/COLOC/git.local.COLOC-reporter/scripts/helper_functions"
  
  opt$GWASsumplots <- c("PIP,p,beta")
  opt$GWASsumplots_onefile <- FALSE
  opt$GWASsumplots_incfinemapping <- TRUE
  
  opt$rdsOut <- "~/Downloads"
  
  opt$outdir <- "/Users/tom/OneDrive - King's College London/PhD/PhD project/COLOC/git.local.COLOC-reporter/testing/tidy_processingTEST_PD.SZ.chr17_coloc/plots"
  
}

#Declare trait IDs and, if provided, name for each trait

#Extract names of traits compared, first dropping file path, and then any prefixes indicated by a preceding underscore
traits <- strsplit(opt$traits,",")[[1]][1:2]
cat("Plotting summary plots for trait IDs: ", paste0(1:length(traits)," = ",traits,collapse = ", "),".\n",sep='')

if(!is.null(opt$GWASconfig)){
  #Read in the configuration options
  GWASconfig<- fread(opt$GWASconfig)
  
  #Assign names to traits
  names(traits) <- GWASconfig[traits,on="ID"]$traitLabel
}

#Read in helper functions
list.files(opt$helperFunsDir,full.names = TRUE,pattern=".R") %>%
  lapply(.,source) %>%
  invisible(.)


#Read in dataset
data<- read.table(opt$harmonisedSumstats,sep=",",header=TRUE)

#Filter traits to match only those indicated in the --traits option
data<- data[data$trait %in% traits,]
if(nrow(data)==0){
  stop("No records remain in the data supplied to --harmonisedSumstats after filtering to records where the 'trait' column corresponds to the first two traits indicated in the --traits option")
}

######
### Prepare for, and then generate Summary plots comparing the region-results for the two GWAS
######

#Parse the comma delimited list
opt$GWASsumplots <- strsplit(opt$GWASsumplots,",")[[1]]

#Check for appropriate format across non-predetermined options
checksumplots<- tolower(opt$GWASsumplots) %in% c('pip','p','z','beta')
if(any(!checksumplots)){
  customStrings<- sapply(opt$GWASsumplots[!checksumplots],function(x,data){
    !is.null(data[x]) && is.numeric(data[[x]])  #Return true if the custom string corresponds to a numeric column - otherwise declare FALSE, to be removed
  },data)
  
  if(any(!customStrings)){
    warning("Summary plots cannot be produced for the custom string(s) ", paste0(names(customStrings)[!customStrings],collapse=", ") ," since corresponding numeric columns were not detected within the dataset.")
  }
  
  #Create a new list only retaining the recognised strings and any others which match to a column
  opt$GWASsumplots <- c(opt$GWASsumplots[checksumplots],opt$GWASsumplots[!checksumplots][customStrings])
}

## Check for presence of required parameters
if("pip" %in% tolower(opt$GWASsumplots) && is.null(data$variable_prob)){
  warning("Finemapping PIPs cannot be plotted because the 'variable_prob' column output by SuSiE is absent from the dataset")
  opt$GWASsumplots <- opt$GWASsumplots[!grepl("pip",tolower(opt$GWASsumplots))]
}

if("beta" %in% tolower(opt$GWASsumplots) && is.null(data$beta)){
  warning("Summary statistics betas cannot be plotted because the expected 'beta' column is absent from the dataset")
  opt$GWASsumplots <- opt$GWASsumplots[!grepl("beta",tolower(opt$GWASsumplots))]
}

if("p" %in% tolower(opt$GWASsumplots) && is.null(data$pvalues)){
  warning("Summary statistics p values cannot be plotted because the expected 'pvalues' column is absent from the dataset")
  opt$GWASsumplots <- opt$GWASsumplots[!grepl("p",tolower(opt$GWASsumplots))]
}


if("z" %in% tolower(opt$GWASsumplots) && any(is.null(data$beta),is.null(data$SE)) || sum(grepl("^(z|Z)$",colnames(data)))!=1){
  warning("Summary statistics Z-scores cannot be plotted because the required columns were missing from the dataset. Please provide either the column 'Z' (case insensitive) or both of 'beta' and 'SE'.")
  opt$GWASsumplots <- opt$GWASsumplots[!grepl("z",tolower(opt$GWASsumplots))]
}

if(length(opt$GWASsumplots)==0){
  stop("No valid options for summary plotting remain. Terminating script")
}



#Set colour attribute 
if(!is.null(data$cs) && opt$GWASsumplots_incfinemapping){
  colourMapping <- "cs" 
  data$cs <- factor(data$cs)
} else if (is.null(data$cs) && opt$GWASsumplots_incfinemapping){
  warning("Credible sets cannot be shown in the figure since the 'cs' column is missing from the dataset.")
  colourMapping = NULL #Set colour attribute 
} else {
  colourMapping = NULL #Set colour attribute 
}

#Setup a colour palette to be used across ggplots (colours taken from plot SuSiE)
ggpalette = c("dodgerblue2", "green4", "#6A3D9A", "#FF7F00", 
              "gold1", "skyblue2", "#FB9A99", "palegreen2", "#CAB2D6", 
              "#FDBF6F", "gray70", "khaki2", "maroon", "orchid1", "deeppink1", 
              "blue1", "steelblue4", "darkturquoise", "green1", "yellow4", 
              "yellow3", "darkorange4", "brown")

#Dynamically plot gwas stat summary figures, varying the y-axis accordingly. This is subsequently called within lapply
#Internally, recode preset options according to the element specified to 'yaxis'.
ggsummaryplot <- lapply(opt$GWASsumplots,ggSummaryplot,
                        dset=data,
                        xlim=range(data$pos),
                        chr=data$chr[1],
                        colourMapping=colourMapping,
                        traits=traits,
                        facetNrow=ifelse(opt$GWASsumplots_onefile,2,1),
                        nameColourLegend="Trait: Credible set",
                        facetTraits=TRUE,
                        build=paste0("GRCh",opt$genomeAlignment),
                        compareTraits=TRUE)
names(ggsummaryplot)<- opt$GWASsumplots

#Write individual plots to rds file
if(!is.null(opt$rdsOut)){ saveRDS(ggsummaryplot,file.path(opt$rdsOut,"SummaryFigs.Rds")) }

#If more than one file requested, concatenate
if(opt$GWASsumplots_onefile && length(ggsummaryplot)>1){
  
  #Drop repeated x-axis text from the bp figures
  index=1:(length(ggsummaryplot)-1)
  ggsummaryplot[index] <-lapply(ggsummaryplot[index],function(x){
    x$bpfigure <-x$bpfigure+
      theme(axis.title.x = element_blank())
    return(x)
  })
  
  #Align all the bp figures and traitxy figures column-wise to ensure x-axes panels are aligned correctly, then align the two columns
  wrapOneFile <- wrap_plots(
    
    wrap_plots(lapply(ggsummaryplot,function(x)x$bpfigure),ncol=1),
    wrap_plots(lapply(ggsummaryplot,function(x)x$traitxy_figure),ncol=1),
    ncol = 2)+
    plot_layout(guides = 'collect')
  
  #Save the combined plots to file
  height=75*length(ggsummaryplot)
  ggsummaryfile <- ggsave(file.path(opt$outdir,paste0(paste0(traits,collapse="_"),"_summary_plots.pdf")),
                          wrapOneFile,device="pdf",units="mm",width=200,height=height)
  
  
} else {
  
  #Patchwork wrap the plots individually
  wrapped <- lapply(ggsummaryplot,function(x){
    wrap<- x$bpfigure / x$traitxy_figure+
      plot_layout(guides = 'collect',heights = c(1,2))
    
    return(wrap)
  })
  
  for(i in 1:length(wrapped)){
    ggsave(file.path(opt$outdir,paste0(paste0(traits,collapse="_"),"_",names(wrapped)[i],"_summaryplot.pdf")),
           wrapped[[i]],device="pdf",units="mm",width=175,height=175)
  }
}

