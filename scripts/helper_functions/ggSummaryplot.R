#####
# Author: Thomas Spargo (thomas.spargo@kcl.ac.uk)
# Obtained from GitHub repository: https://github.com/ThomasPSpargo/COLOC-reporter
#
# Defines a custom plotting function to be recycled for (optional) summary plots and later for plotting snp pp's
#####

#Inputs
#'chr' refers to chromosome and should be provided to get a complete x-axis title
#These should be provided to get a complete x-axis title. The function will fall back upon basepair positions

#Define a custom plotting function to be recycled for (optional) summary plots and later for plotting snp pp's
ggSummaryplot <- function(yaxis,bp_range,dset,colourMapping=NULL,shapeMapping=NULL,traits=NULL,
                          facetTraits=FALSE,facetNrow=NULL,nameColourLegend=NULL,nameShapeLegend=NULL,chr=NULL,compareTraits=FALSE){
  
  #Setup y-axis parameters (and where required adjust data)
  if(tolower(yaxis)=="p"){
    ylab <- -log[10]~p
    ystring <- "pvalues"
    dset$pvalues <- -log10(dset$pvalues) #Convert p-values
  } else if (tolower(yaxis) %in% c("z","abs_z")){
    ystring <- "z"
    
    #Either calculate Z from beta and se, or grep the column name
    if(!"z" %in% tolower(colnames(dset))){
      dset$z <- dset$beta/dset$SE
    } else {
      dset$z <- dset[[grep("^(z|Z)$",colnames(dset))]]
    }
    
    if(tolower(yaxis)=="abs_z"){ #Conditionally convert to absolute z-scores.
      ylab <- "Absolute Z-score"
      dset$z <- abs(dset$z)
    } else {
      ylab <- "Z-score"
    }
  } else if(tolower(yaxis)=="pip"){
    ylab <-  "PIP"
    ystring <- "variable_prob"
  } else if(tolower(yaxis)=="snp.pp"){
    ylab <- "Posterior probability for being shared variant"
    ystring <- "SNP.PP"
  } else if(tolower(yaxis)=="beta"){
    ylab <- ~beta
    ystring <- "beta"
  } else {
    ylab <- yaxis
    ystring <- yaxis
  }
  
  if(facetTraits && !is.null(traits)){
    dset$trait<- factor(dset$trait,levels=traits)
    if(!is.null(names(traits))) levels(dset$trait) <- names(traits)
  } 
  
  #Rescale x-axis into MBs or KBs
  if((max(bp_range)-min(bp_range))>100000){
    xscale <- 1000000
    xscale_magnitude <- " [Mb]"
  } else if ((max(bp_range)-min(bp_range))>100) {
    xscale <- 1000
    xscale_magnitude <- " [Kb]"
  } else {
    xscale <- 1
    xscale_magnitude <- ""
  }
  
  #Setup x-axis label. If Chr is provided try to give a complete legend, otherwise go just by x-limits
  xrange<- paste0(min(bp_range),"-",max(bp_range))
  if(!is.null(chr)){
    xlab<- paste0("Genomic position",xscale_magnitude,"\n(Chr",chr,":",xrange,")")
  } else {
    xlab<- paste0("Position",xscale_magnitude," (",xrange,")")
  }
  
  if(!is.null(colourMapping)){ #Arrange data so that coloured points will always be plotted atop the NA values
    dset <- dset %>%
      arrange(desc(is.na(.data[[colourMapping]])),.data[[colourMapping]])
  }
  
  #Logical to indicate whether the plot is for a measure with negative values - if yes draw a hatched line along 0
  incNegvals=min(dset[[ystring]])<0
  
  
  
  
  #Supply string-based aesthetic elements via list
  aesthetics<- list(x="pos",y=ystring,colour=colourMapping,shape=shapeMapping) %>%
    lapply(., function(x) if(!is.null(x)){sym(x)})
  
  #Generate plot
  pos_fig <- ggplot(dset,aes(!!!aesthetics))+
    geom_point(na.rm=TRUE)+
    theme_bw()+
    scale_x_continuous(limits = bp_range,
                       breaks = scales::breaks_extended(n=4), 
                       labels = scales::label_number(scale = 1 / xscale,accuracy = 0.01) #Scale axis magnitude dynamically
    )+ 
    labs(y=bquote(.(ylab)),x=xlab)+
    { if(incNegvals) geom_hline(yintercept = 0,lty=2) }+
    { if(facetTraits) facet_wrap(~trait, nrow = facetNrow) }+ #Optionally add faceting
    
    { if(!is.null(colourMapping)){  #Optionally add colour aesthetic
      list(scale_colour_manual(na.value = "black", values=ggpalette,
                               breaks=levels(dset[[which(colnames(dset)==colourMapping)]])
      ),
      guides(color=guide_legend(title=nameColourLegend,order=1))) } } +
    
    { if(!is.null(shapeMapping)){ #Optionally add shape aesthetic 
      list(scale_shape_manual(values=c(17,19),breaks=levels(dset[[which(colnames(dset)==shapeMapping)]])),
                              guides(shape=guide_legend(title=nameShapeLegend,order=2))) } }
  
  
  ## Optionally derive a further plot comparing a pair of traits along the x and y axis
  if(compareTraits){
    if(class(dset$trait)=="factor") traits<- levels(dset$trait) else traits<- unique(dset$trait)
    aesthetics$x <- sym(traits[1])
    aesthetics$y <- sym(traits[2])
    
    lims <- range(dset[[ystring]])  
    
    traitXY_fig<- dset %>%
      dplyr::select(all_of(c("snp","trait",ystring,colourMapping))) %>%
      pivot_wider(values_from = all_of(ystring),names_from = "trait") %>%
      ggplot(.,aes(!!!aesthetics))+
      lims(x=lims,y=lims)+
      labs(x=bquote({.(ylab)[~(.(traits[1]))]},splice=TRUE),y=bquote({.(ylab)[~(.(traits[2]))]},splice=TRUE))+
      geom_point(na.rm = TRUE)+
      theme_bw()+
      { if(incNegvals) list(geom_hline(yintercept = 0,lty=2),
                            geom_vline(xintercept = 0,lty=2))}+
      { if(!is.null(colourMapping)){  #Optionally add colour aesthetic 
        list(scale_colour_manual(na.value = "black", values=ggpalette,
                                 breaks=levels(dset[[which(colnames(dset)==colourMapping)]])
        ),
        guides(color=guide_legend(title=nameColourLegend,order=1))) } } +
      { if(!is.null(shapeMapping)){ #Optionally add shape aesthetic 
        list(scale_shape_manual(values=c(17,19),breaks=levels(dset[[which(colnames(dset)==shapeMapping)]]),
                                guides(shape=guide_legend(title=nameShapeLegend,order=2)))) } }
    
  
    return(list(bpfigure=pos_fig,traitxy_figure=traitXY_fig))
  } else {
    return(list(bpfigure=pos_fig))
  }
}