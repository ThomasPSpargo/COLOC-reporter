#####
# Author: Thomas Spargo (thomas.spargo@kcl.ac.uk)
# Obtained from GitHub repository: https://github.com/ThomasPSpargo/COLOC-reporter
#
# This script visualises sensitivity analysis of COLOC results
# It is a wrapper around the coloc::sensitivity function and generates ggplots as opposed to base R graphics for the posterior probability plot
#####

suppressPackageStartupMessages(library(optparse))

option_list = list(

  #General options
  make_option("--colocResultsDir", action="store", default=NULL, type='character',
              help="Directory containing .Rds files which store outputs from coloc.abf and/or coloc.susie, across which analysis sensitivity can be plotted. Traits compared are identified based on file name, expecting the format *susie_<Trait1>_<Trait2>.Rds for a coloc.susie output and *abf_<Trait1>_<Trait2>.Rds for a coloc.abf output."),
  make_option("--rule", action="store", default=NULL, type='character',
              help="From the coloc::sensitivity function documentation: 'a decision rule. This states what values of posterior probabilities \"pass\" some threshold. This is a string which will be parsed and evaluated, better explained by examples. \"H4 > 0.5\" says post prob of H4 > 0.5 is a pass. \"H4 > 0.9 & H4/H3 > 3\" says post prob of H4 must be > 0.9 AND it must be at least 3 times the post prob of H3.\"'"),
  make_option("--drawHline", action="store", default=NULL, type='numeric',
              help="Draw a horizontal line on ggplots (e.g. to indicate the decision rule boundary) [optional]"),
  make_option("--rdsOut",action="store", default=NULL, type='character',
              help="Specify a file path and prefix for an Rds file containing figures generated within this script '.Rds' will be appended [optional]"),
  make_option("--rdsOnly", action="store", default=FALSE, type='logical',
              help="Logical, defaulting to FALSE, if set to TRUE, plots will only output as .Rds files; see the --rdsOut option."),
  make_option("--out", action="store", default=NULL, type='character',
              help="Path and prefix for output figure.")
)

opt = parse_args(OptionParser(option_list=option_list))


suppressPackageStartupMessages({
  library(coloc)
  library(tidyverse)
  library(patchwork)
})

if(is.null(opt$out) && is.null(opt$rdsOut)) { cat("ERROR: Please specify arguments to either or both of --rdsOut and --out. Terminating script\n"); q("no") }
if(opt$rdsOnly && is.null(opt$rdsOut)) warning("Rds files will not be saved the --rdsOut option has not been set.")


######
### Read-in
######
test <- FALSE
if(test){
  opt <- list()
  opt$colocResultsDir <-  "~/Downloads/colocalisation"
  opt$rule <-"H4 > 0.8"
  opt$drawHline <- 0.8
}
  
######
### Setup
######

#Identify and read-in the Rds files
files <- list.files(opt$colocResultsDir,full.names = TRUE,pattern=".Rds$")
colocRes <- lapply(files,readRDS)

#Identify the trait-pair analysed per-file
traits<- gsub(".*(abf|susie)_(.*).Rds","\\2",basename(files)) %>%
  strsplit(.,split="_")

#Define labels for all hypotheses
hypotheses <- paste0("H",0:4)
names(hypotheses)<- paste0("PP.H",0:4,".abf")

#coloc::sensitivity function rule checking method:
# rule <- gsub("(H.)","PP.\\1.abf",opt$rule,perl=TRUE)
# check <- function(pp) { with(as.list(pp),eval(parse(text=rule))) }
# check(x)

plotlist <- list()
for(i in seq_along(colocRes)){
  
  nrow<- ifelse("coloc_abf" %in% class(colocRes[[i]]),1,nrow(colocRes[[i]]$summary))
  
  for(j in 1:nrow){
    #Run sensitivity analysis, saving the posterior probability matrix 
    x <- sensitivity(colocRes[[i]],rule=opt$rule,row=j,
                     doplot = FALSE,plot.manhattans=FALSE)
    
    #Generate the plot
    plotlist[[length(plotlist)+1]] <- x %>%
      pivot_longer(cols=!all_of(c("p12","pass")),names_to = "Hypothesis",values_to = "PP") %>%
      mutate(Hypothesis=recode(Hypothesis,!!!hypotheses)) %>% #Prepare pretty hypothesis names
      ggplot(.,aes(x=p12,y=PP,colour=Hypothesis))+
      {if(any(x$pass)) geom_rect(aes(xmin=min(p12[pass]),xmax=max(p12[pass])),ymin=-Inf,ymax=Inf,alpha=0.1,fill="wheat1",colour="wheat1") }+
      {if(!is.null(opt$drawHline)) geom_hline(yintercept = opt$drawHline,colour="black") }+
      geom_vline(xintercept = colocRes[[i]]$priors["p12"],colour="grey70",lty=2)+
      geom_point(na.rm=TRUE)+
      ylim(c(0,1))+
      scale_colour_brewer(palette = "Set1",labels=scales::parse_format())+
      scale_x_continuous(trans="log10",
                         expand=expansion(mult = 0),
                         breaks=scales::trans_breaks('log10', function(x) 10^x),
                         labels = function(x) parse(text=gsub("e", "%*%10^",format(x, scientific = TRUE)))
      )+
      theme_bw()+
      theme(panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank(),
            axis.text.x = element_text(hjust=0.75),
            plot.subtitle = element_text(hjust = 1))+
      labs(y="Posterior probability for hypothesis",x="Prior probability of shared variant",
           subtitle = ifelse("coloc_abf" %in% class(colocRes[[i]]),
                             paste0("coloc.abf (", paste0(traits[[i]],collapse= " & "),")"),
                             paste0("coloc.susie (", traits[[i]][1],":", colocRes[[i]]$summary$idx1[j], " & ", traits[[i]][2],":",colocRes[[i]]$summary$idx2[j],")")
           )
      )
    
  }
}

#Write plots to Rds file
if(!is.null(opt$rdsOut)){ 
  if(!dir.exists(dirname(opt$rdsOut))) dir.create(dirname(opt$rdsOut),recursive = TRUE)
  saveRDS(plotlist,paste0(opt$rdsOut,".Rds"))
  
  if(opt$rdsOnly) q(save="no")
}

if(length(plotlist)>1){
  #Combine into a single figure
  plots<- Reduce("+",plotlist)+
    plot_layout(guides="collect",ncol=switch(length(plotlist) %in% 2:4,2,NULL))
  
  plots<- plots & theme(
    plot.margin = unit(c(1,3,1,1), "mm"),
    legend.margin = margin(0,0,0,0, "mm")
  ) #Add space to the right of panels to avoid cropping the x-axis text
  
  width=200
  height=ifelse(length(plotlist)==2,100,200)
  
} else {
  plots <- plotlist[[1]]
  
  width=120
  height=100
}

ggsave(paste0(opt$out,".pdf"),plots,device="pdf",units="mm",width=width,height=height)
