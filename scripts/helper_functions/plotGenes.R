#####
# Author: Thomas Spargo (thomas.spargo@kcl.ac.uk)
# Obtained from GitHub repository: https://github.com/ThomasPSpargo/COLOC-reporter
#
# This script defines a function for visualising genomic position for genes occurring across a given region. 
# Attempts to plot a large number of genes may not be suitable for this plot - the gene_tracks option may help 
#####

#inputs
#bp_range - the range of base pair values that approximately correspond to the positions of the genes list, used to set plot limits
#chr - the chromosome to which bp_range corresponds
#genes: a list of genes including the columns "external_gene_name", "chromosome_name","start_position","end_position", "external_gene_source", "start_window", "end_window"
  #start_window and end_window define slight extensions from start_position and end_position to ensure nearby genes are correctly captured
#geneSources: Restrict plotting to certain sources by supplying a vector of strings (e.g. "HGNC Symbol") which matches with the "external_gene_source" column of genes The setting is ignored if nothing remains for plotting after filtering
#gene_tracks: indicate an integer which acts as a threshold for declaring that, if plotting more genes than this threshold, the y axis should no longer include a row per gene and instead attempt to plot multiple genes across 'tracks'. When using tracks, gene names are plotted as labels in the panel with nudge_y controlling their positioning. Without tracks, y-axis labels declare genes
#nudge_y: controls the y-axis position of text labels if using the gene_tracks option

plotGenes <- function(bp_range,chr,genes,geneSources=NULL,gene_tracks=NULL,nudge_y=0.2){
  
  if(nrow(genelist)==0){
    stop("The list supplied to the 'genes' argument is empty; therefore nothing remains to be plotted.")
  }

  #Arrange the genes list
  genelist <- genes %>%
    arrange(start_position) %>%                      #Sort by start position
    mutate(midpoint=(start_position+end_position)/2, #Determine middle of gene position
           external_gene_name=factor(external_gene_name,levels = rev(external_gene_name[!duplicated(external_gene_name)]))) #Save as factor variable, to maintain sort order
  
  #If subsetting to named gene sources
  if(!is.null(geneSources)){
    
    sourceSubset<- which(genelist$external_gene_source %in% geneSources)        #Find matching rows
    if(length(sourceSubset)>0){                                                     #If at least one gene remains, subset
      genelist<- genelist[sourceSubset,]
    } else {
      warning("Could not subset nearby gene plotting by source since no genes remained after restricting to ", paste0(geneSources,collapse=", "), "\n")
    }
  }
  
  #Establish x-axis for the plot  
  plot_rng <- c(min(c(min(bp_range),genelist$start_window)),max(c(max(bp_range),genelist$end_window)))
  
  #Determine x-scale magnitude conversion
  xScale<-xScaler(plot_rng[2]-plot_rng[1])
  
  if(!is.null(gene_tracks) && nrow(genelist)>=gene_tracks){
    
    #Create y axis 'tracks' along which genes are plotted. This recurses across gene start and end positions aiming to plot genes with >5% of the xlim space between each gene
    genespacing<- (plot_rng[2]-plot_rng[1])/5
    
    genelist$track <- 1           #Assign all genes to track 1 initially
    for(t in 2:nrow(genelist)){   #For each gene, shift tracks upwards to avoid overlap with those in current track
      prior <- genelist[1:(t-1),]
      
      r=1                             #Start from track number 1, and while overlapping, increase tracks until adequate plotting space is found
      inTrack<- which(prior$track==r) #identify rows in current track
      while(length(inTrack)>0){       #While in an existing track
        if(genelist$start_position[t]>=(max(prior$end_position[inTrack])+genespacing)){ #If the start position is more than 100Kb apart from the previous gene in track, then assign gene to track
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
    labs(x=paste0("GRCh",opt$genomeAlignment," genomic position",xScale$xscale_magnitude,"\n(Chr", reg_range["chr"],":",plot_rng[1],"-",plot_rng[2],")"),
         y="Nearby genes")+
    { if(useTracks){
      list(geom_text(aes(label=external_gene_name),fontface="italic",nudge_y = nudge_y,size=3),
           theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
      ) } else { theme(axis.text.y = element_text(face="italic")) } }
      
  
  return(list(plot=fig,xlims=plot_rng))
}
    
