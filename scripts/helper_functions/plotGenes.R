#####
# Author: Thomas Spargo (thomas.spargo@kcl.ac.uk)
# Obtained from GitHub repository: https://github.com/ThomasPSpargo/COLOC-reporter
#
# This script defines two functions.
# The first visualises genomic position for genes occurring across a given region. 
# The second is a wrapper around the function and uses patchwork to stack the plot output by the first function beneath a second panel with an x-axis encompassing the genes plot range
# Attempts to plot a large number of genes may not be suitable for this plot - the gene_tracks option may help 
#####

#inputs
#genes: a list of genes including the columns "external_gene_name", "chromosome_name","start_position","end_position", "external_gene_source", "start_window", "end_window"
#start_window and end_window define slight extensions from start_position and end_position to ensure nearby genes are correctly captured
#bp_range - optionally, a range of base pair values that approximately correspond to the positions of the genes list, used as a variable for setting plot limits, taking into consideration also the start and end windows for the supplied gene list
#geneSources: Restrict plotting to certain sources by supplying a vector of strings (e.g. "HGNC Symbol") which matches with the "external_gene_source" column of genes The setting is ignored if nothing remains for plotting after filtering
#gene_tracks: indicate an integer which acts as a threshold for declaring that, if plotting more genes than this threshold, the y axis should no longer include a row per gene and instead attempt to plot multiple genes across 'tracks'. When using tracks, gene names are plotted as labels in the panel with nudge_y controlling their positioning. Without tracks, y-axis labels declare genes
#nudge_y: controls the y-axis position of text labels if using the gene_tracks option

plotGenes <- function(genes,bp_range=NULL,geneSources=NULL,gene_tracks=NULL,nudge_y=0.2){
  
  if(nrow(genes)==0) stop("The list supplied to the 'genes' argument is empty; therefore nothing remains to be plotted.")
  
  #Arrange the genes list
  genelist <- genes %>%
    arrange(start_position) %>%                      #Sort by start position
    mutate(midpoint=(start_position+end_position)/2, #Determine middle of gene position
           external_gene_name=factor(external_gene_name,levels = rev(external_gene_name[!duplicated(external_gene_name)]))) #Save as factor variable, to maintain sort order
  
  #If subsetting to named gene sources
  if(!is.null(geneSources)){
    sourceSubset<- which(genelist$external_gene_source %in% geneSources)        #Find matching rows
    if(length(sourceSubset)>0){ genelist<- genelist[sourceSubset,]              #If at least one gene remains, subset
    } else { warning("Could not subset nearby gene plotting using the geneSources argument since no genes remained after restricting to ", paste0(geneSources,collapse=", "), ". Plotting with all available gene sources.\n") }
  }
  
  #Establish x-axis for the plot  
  if(!is.null(bp_range)){ plot_rng <- c(min(c(min(bp_range),genelist$start_position)),max(c(max(bp_range),genelist$end_position)))
  } else { plot_rng <- c(min(genelist$start_position),max(genelist$end_position)) }
  
  #Determine x-scale magnitude conversion
  xScale<-xScaler(plot_rng[2]-plot_rng[1])
  
  if(!is.null(gene_tracks) && nrow(genelist)>=gene_tracks){
    
    #Create y axis 'tracks' along which genes are plotted. This recurses across gene start and end positions aiming to plot genes with >20% of the xlim space between each gene
    genespacing<- (plot_rng[2]-plot_rng[1])/5
    
    genelist$track <- 1           #Assign all genes to track 1 initially
    for(t in 2:nrow(genelist)){   #For each gene, shift tracks upwards to avoid overlap with those in current track
      prior <- genelist[1:(t-1),]
      
      r=1                             #Start from track number 1, and while overlapping, increase tracks until adequate plotting space is found
      inTrack<- which(prior$track==r) #identify rows in current track
      while(length(inTrack)>0){       #While in an existing track
        if(genelist$start_position[t]>=(max(prior$end_position[inTrack])+genespacing)){ #If the start position is sufficiently apart from the previous gene in track, then assign gene to track
          genelist$track[t] <- r
          break #break loop
          
        } else { #increase track number
          r <- r+1
          inTrack<- which(prior$track==r) #identify rows in current track
          if(length(inTrack)==0){ #If this is a new track, assign final track number and loop will break
            genelist$track[t] <- r
          }
        }
      }
    } #End tracking loop
    
    #Convert to factor in order to display lower tracks higher
    tracklevels<- rev(sort(unique(genelist$track)))
    genelist$track <- factor(genelist$track,levels=tracklevels) 
    
    useTracks=TRUE
  } else {
    useTracks=FALSE
  }
      
  #Generate the plot, with conditional layers according to whether plotting is performed based on tracks or otherwise
  fig <- ggplot(genelist,aes(x=midpoint,y=!!sym(ifelse(useTracks,"track","external_gene_name")))) +
    geom_errorbar(aes(xmin=start_position,xmax=end_position),width=0)+
    theme_bw()+
    scale_x_continuous(limits = plot_rng,
                       breaks = scales::breaks_extended(n=4), 
                       labels=scales::label_number(scale = 1 / xScale$xscale,accuracy = 0.01))+ #Scale axis magnitude dynamically
    labs(x=paste0("Genomic position",xScale$xscale_magnitude,"\n(Chr",genelist$chromosome_name[1],":",plot_rng[1],"-",plot_rng[2],")"),
         y="Nearby genes")+
    { if(useTracks){
      list(geom_text(aes(label=external_gene_name),fontface="italic",nudge_y = nudge_y,size=3),
           theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
      ) } else { theme(axis.text.y = element_text(face="italic")) } }
      
  
  return(list(plot=fig,xlims=plot_rng))
}

#Plot genes nearby to a range of positions, then stack the original plot above the nearby gene plot
#Note that genesInChr is expected to already have filtered to the correct chromosome
#highConfidence declares how nearby genes filtering should be handled, set TRUE to only return genes containing the full bp region
#... passes options onto plotGenes
plotGenesStack <- function(snpPP,genesInChr,bp_range,highConfidence=FALSE,...){
  
  #Filter to genes nearby to bp_range
  Genes_subset <- getNearbyGenes(bp_range=bp_range,genesInChr=genesInChr,highConfidence=highConfidence)
  
  if(nrow(Genes_subset)>0){
    
    #Plot nearby genes
    gene_near <- plotGenes(bp_range=range(bp_range),genes=Genes_subset,...)
    
    #Narrow the snp plot xlim to match gene_near
    suppressMessages({ #Suppress warning about 'double-setting' xlim
      snpPP_range <- snpPP +
        scale_x_continuous(limits = gene_near$xlim,
                           breaks = scales::breaks_extended(n=4))+
        theme(axis.title.x = element_blank(),
              axis.text.x = element_blank())
    })
    
    #Update any colour-mapping variable to reflect region visualised
    if(!is.null(snpPP$mapping$colour)){
      #To ensure consistent colouring when reassigning the colour factor, assign factor levels to specific colours
      levels<- levels(snpPP$data[[snpPP$labels$colour]])
      col_levels <-ggpalette[1:length(levels)]
      names(col_levels)<- levels
      
      suppressMessages({ #Suppress warning about redoing existing legend
        snpPP_range <- snpPP_range +
          scale_colour_manual(na.value = "black", values=col_levels,
                              breaks = levels(droplevels(snpPP$data[[snpPP$labels$colour]][with(snpPP$data,which(pos > gene_near$xlims[1] & pos < gene_near$xlims[2]))])))
      })
    }
    
    #Stack the snp plot above the nearby gene plot, adjust the plot stack proportions according to n genes plotted and the way that labels are plotted
    if(is.null(gene_near$plot$data$track)){ propor <- 0.05*nrow(gene_near$plot$data)
    } else { propor <- 0.2*length(levels(gene_near$plot$data$track)) }
    if(propor>1){propor <- 1}
    
    stackPlot <- (snpPP_range/gene_near$plot)+
      plot_layout(guides = 'collect',heights=c(1,propor))
    
    return(stackPlot)
  } else {
    return("No nearby genes could be plotted")
  } 
}

    
