#####
# Author: Thomas Spargo (thomas.spargo@kcl.ac.uk)
# Obtained from GitHub repository: https://github.com/ThomasPSpargo/COLOC-reporter
#
# The purpose of this script is to generate ggplots visualising the output of genetic correlation analysis using LAVA
# The script expects input files to follow with the format 'Trait1.Trait2.<ext>' where <ext> is 'univ' for a univariate and bivar for a bivariate analysis.
#####


#Rscript to generate ggplots for lava correlations
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(ggrepel))
suppressMessages(library(optparse))

option_list = list(
  make_option("--infile", action="store", default=NULL, type='character',
              help="The file path to the LAVA output file, in the format 'Trait1.Trait2.<ext>' where <ext> is 'univ' for a univariate and bivar for a bivariate analysis. and Trait1/2 are IDs for the traits compared.
              "),
  make_option("--outdir", action="store", default=NULL, type='character',
              help="Directory in which to return output plots"),
  make_option("--global_rg", action="store", default=NULL, type='numeric',
              help="Numeric indicating the global RG for the current comparison"),
  make_option("--pthresh", action="store", default=NULL, type='numeric',
              help="Manually assign p-value significance thresholds for univariate and bivariate analyses as a comma separated list (e.g. '0.0004,0.001'); the first numeric will be applied to univariate analyses, while the second is applied to bivariate analyses."),
  make_option("--suppressLabels", action="store", default=FALSE, type='logical',
              help="Logical to indicate whether or not to supress ggrepel label plotting. TRUE will produce plots with no labels")
)

opt = parse_args(OptionParser(option_list=option_list))

options(echo=FALSE)
#Plot adapted from
#https://danielroelfs.com/blog/how-i-create-manhattan-plots-using-ggplot/

#Create directory if non-existent - This should already be created when originally calling the script
if(!dir.exists(opt$outdir)){dir.create(opt$outdir,recursive = TRUE)}

#Define the 'type' object based on whether input is bivar or univ which is used to determine which plots are drawn
if(grepl("univ",opt$infile)){
  type <- "univ"
} else if(grepl("bivar",opt$infile)){
  type <- "bivar"
} else {
  error("The input was not recognised as being either 'univ' or 'bivar' please supply either a univariate or bivariate analysis, indicating the analysis type in the file extension.")
}

#Gsub out filepath and identify phenotypes compared in the plot
phenos<- gsub(".*/","",gsub("\\.(univ|bivar)","",opt$infile))

#Import dataset 
data <- read.delim(opt$infile,sep="\t",header = TRUE)


#Calculate x-axis positions relevant for each chromosome
axis_set <- data %>% 
  group_by(chr) %>% 
  summarize(center = mean(locus))

if(!is.null(opt$pthresh)){ #Split the P-threshold string and assign to univ and bivar analyses accordingly
  pthresh <- as.numeric(strsplit(opt$pthresh,split=",")[[1]]) 
  
  if(type=="univ"){ 
    sig <- pthresh[1]
  } else if(type=="bivar"){
    sig <- pthresh[2]
  }
} else {
  
  if(type=="univ"){
    ncompar <- 2495
  } else if(type=="bivar"){
    ncompar <- nrow(data)
  }
  sig <- 0.05/ncompar #Identify significance threshold
  
}
#Add labels to loci smaller than sig threshold
data <- data %>%
  mutate(thresh= if_else(p < sig, as.character(locus),"")) 

#Plot pvalues in manhattan plot
manhp <- ggplot(data,aes(x=locus,y=-log10(p),color=chr %% 2==0))+
  geom_point()+
  geom_hline(yintercept = -log10(sig), linetype="dashed")+  #plot significant threshold
  scale_color_manual(values = rep(c("#2FD095", "#D02F6A"), length(unique(data$chr))),guide="none")+
  scale_x_continuous(label = axis_set$chr, breaks = axis_set$center) +
  labs(x = "Chromosome", 
       y = bquote(-log[10]~p))+
  theme_bw()

if(opt$suppressLabels==FALSE){
  manhp <- manhp+
    geom_text_repel(aes(label=thresh))
}

if(type=="univ"){
  #Facet univariate plots by phenotype
  manhp <- manhp+
    facet_wrap(~phen,ncol=1)
} 

ggsave(manhp, 
       filename = paste0(opt$outdir,phenos,".",type,".manhat.pdf"),
       units="mm",width=175,height=150)

if(type=="bivar"){
  
  ####Adjust locus scaling to create larger x-axis gaps when moving to a new chromosome
  point <- 1 #Set first x axis value
  for(i in 1:length(data$locus)){
    
    #If the new value of i is a new chromosome, add spacing to the locus, else make the locus value 1+
    if(data$chr[i]!=data$chr[i-1] && i!=1){
      point <- point+5
    } else {
      point <- point+1
    }
    
    #Assign the value of 'point' to the locus vector 
    data$locus[i] <- point
    
  }
  
  #Calculate x-axis positions relevant for each chromosome
  axis_set <- data %>% 
    group_by(chr) %>% 
    summarize(center = mean(locus))
  
  #Plot row
  rhop <- ggplot(data,aes(x=locus,y=rho,color=chr %% 2==0))+
    geom_point()+
    geom_errorbar(aes(ymin=rho.lower,ymax=rho.upper),width=0,
                  position=position_dodge(width = 0.6))+
    #geom_hline(yintercept=0)+  #Plot a 'no association line
    geom_hline(aes(yintercept = mean(rho),linetype="LAVA average"), linetype="dotted", color = "grey40",show.legend =TRUE)+  #plots of mean association
    scale_color_manual(values = rep(c("#2FD095", "#D02F6A"), length(unique(data$chr))),guide="none")+
    scale_x_continuous(label = axis_set$chr, breaks = axis_set$center)+
    geom_text_repel(aes(label=thresh))+
    labs(x = "Chromosome", 
         y = bquote(Genetic~correlation~"[95% CI]"))+
    ylim(c(-1,1))+
    theme_bw() 
  
  if(opt$suppressLabels==FALSE){
    rhop <- rhop +
      geom_text_repel(aes(label=thresh))
  }
  
  #If a global genetic correlation has been specified, plot this as a h-line
  if(!is.null(opt$global_rg)){
    rhop <- rhop +
      geom_hline(aes(yintercept=opt$global_rg,linetype="LDSC"), linetype="longdash", color = "grey40",show.legend =TRUE) #plots of LDSC global association
  }
  
  
  ggsave(rhop, 
         filename = paste0(opt$outdir,phenos,".rgcor.pdf"),
         units="mm",width=300,height=150)
  
}



