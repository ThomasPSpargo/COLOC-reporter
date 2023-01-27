#####
# Author: Thomas Spargo (thomas.spargo@kcl.ac.uk)
# Obtained from GitHub repository: https://github.com/ThomasPSpargo/COLOC-reporter
#
# Defines a custom plotting function to be recycled for (optional) summary plots and later for plotting snp pp's
#####

#Define a custom plotting function to be recycled for (optional) summary plots and later for plotting snp pp's
ggSummaryplot <- function(yaxis,xstring="pos",xlim,dset=minimal_P1P2_df,colourMapping=NULL,shapeMapping=NULL,figdir,traits,returnplot=TRUE,
                          facetTraits=FALSE,facetNrow=1,nameColourLegend="Trait: Credible set",nameShapeLegend=NULL,alignment){
  
  #Setup y-axis parameters (and where required adjust data)
  if(tolower(yaxis)=="p"){
    ylab <- bquote(-log[10]~p)
    ystring <- "pvalues"
    dset$pvalues <- -log10(dset$pvalues) #Convert p-values
  } else if (tolower(yaxis)=="z"){
    dset$z <- abs(dset$beta/dset$SE)
    ylab <- "Absolute Z-score"
    ystring <- "z"
  } else if(tolower(yaxis)=="pip"){
    ylab <-  "PIP"
    ystring <- "variable_prob"
  } else if(tolower(yaxis)=="snp.pp"){
    ylab <- "Posterior probability of being a shared variant\n"
    ystring <- "SNP.PP"
  } else {
    ylab <- yaxis
    ystring <- tolower(yaxis)
  }
  
  if(facetTraits==TRUE){
    dset$trait<- factor(dset$trait,levels=traits,labels = names(traits))
  }
  
  #Rescale x-axis into MBs or KBs
  if((xlim[2]-xlim[1])>100000){
    xscale <- 1000000
    xscale_magnitude <- " [Mb]"
  } else if ((xlim[2]-xlim[1])>100) {
    xscale <- 1000
    xscale_magnitude <- " [Kb]"
  } else {
    xscale <- 1
    xscale_magnitude <- ""
  }
  
  #Supply string-based aesthetic elements via list
  aesthetics<- list(x=xstring,y=ystring,colour=colourMapping,shape=shapeMapping) %>%
    lapply(., function(x) if(!is.null(x)){sym(x)})
  
  #Generate plot
  plot_root <- ggplot(dset,aes(!!!aesthetics))+
    geom_point()+
    theme_bw()+
    scale_x_continuous(limits = xlim,
                       breaks = scales::breaks_extended(n=4), 
                       labels = scales::label_number(scale = 1 / xscale,accuracy = 0.1) #Scale axis magnitude dynamically
    )+ 
    labs(y=ylab,x=paste0("GRCh",alignment," genomic position (",target_region,")",xscale_magnitude))
  
  if(facetTraits==TRUE){
    plot_root <- plot_root +
      facet_wrap(~trait, nrow = facetNrow)
  }
  
  #Add Colouring if specified
  if(!is.null(colourMapping)){
    plot_root <- plot_root +
      scale_colour_manual(na.value = "black", values=ggpalette,
                          breaks=levels(dset[[which(colnames(dset)==colourMapping)]])
      )+
      guides(color=guide_legend(title=nameColourLegend,order=1))
  }
  
  #Add shapes if specified
  if(!is.null(shapeMapping)){
    plot_root <- plot_root +
      #scale_shape_manual(na.value = 19, values=17,breaks=levels(dset[[which(colnames(dset)==shapeMapping)]])
      scale_shape_manual(values=c(17,19),breaks=levels(dset[[which(colnames(dset)==shapeMapping)]])
      )+
      guides(shape=guide_legend(title=nameShapeLegend,order=2))
  }
  
  #Return either the plot itself or ggsave directly, according to return plot option
  if(returnplot){
    return(plot_root)
  } else {
    #GGsave returns character string indicating where file was written
    ggsave(file.path(figdir,paste0(paste0(traits,collapse="_"),"_",tolower(yaxis),"_summaryplot.pdf")),
           plot_root,device="pdf",units="mm",width=175,height=75)
  }
}