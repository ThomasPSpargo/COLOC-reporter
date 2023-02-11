#####
# Author: Thomas Spargo (thomas.spargo@kcl.ac.uk)
# Obtained from GitHub repository: https://github.com/ThomasPSpargo/COLOC-reporter
#
# This script generates heatmaps visualizing LD between credible sets relative to the top snps passed to the function
#####


#dependencies:
# tidyverse
# patchwork


#Inputs
#sets - a dataframe containing information about the credible sets identified by susie. the columns "snp" - ID, corresponding to the topsnps,"cs" - credible set assignments,"pos" - genomic position
#R - an LD correlation matrix associated with susie.fit
#palette - the palette to use on the ggplot heatmap

#topsnps - a dataframe declaring the top snp for each credible set in the column "topsnp" - used as an anchor point and the credible set to which this corresponds in the column 'cs'
#heatmapPalette="OrRd" - specify the colour palette to use for the heatmap
#discretePalette=NULL  - specify the colour palette to use for visualising different credible sets in a bar that corresponds to the heatmap
#alignment=37          - Indicate the genome alingnment. 37 is assumed.
#plotR2=TRUE           - TRUE returns an R2 plot, FALSE returns correlations
#trait                 - trait string id for the finemapped trait
#maxPerCS              - indicate a numeric limit for the number of top 'snps' that will be plotted along the y-axis for each credible set. This is included to avoid overplotting.

susie_cs_ld <- function(sets,R,topsnps,heatmapPalette="OrRd",discretePalette=NULL,alignment=37,plotR2=TRUE,trait,maxPerCS=10){
  
  #Filter step to reduce overplotting if a large number of snps are identified as 'top'
  if(any(tabulate(topsnps$cs)>maxPerCS)){
    message("Randomly downsampling credible sets with >",maxPerCS," top snps identified")
    
    incSets<- unique(topsnps$cs)
    for(i in seq_along(incSets)){
      logiIndex<- topsnps$cs==incSets[i] #Get a logical index of snps in and outside the current CS
      
      #If the cs includes snps above the threshold, filter down
      if(sum(logiIndex)>maxPerCS){
        topsnps<- rbind(topsnps[!logiIndex,],
                 topsnps[sample(which(logiIndex),size=maxPerCS),])
      }
    }
  }
  
  #Extract LD estimates relative to the top SNPs
  topSNPLD <- R[,dimnames(R)[[1]] %in% topsnps$topsnp,drop=FALSE] %>%
    as.data.frame() %>%
    rownames_to_column("snp")
  
  matchCheck <- !topsnps$topsnp %in% colnames(topSNPLD)
  if(any(matchCheck)){
    warning("Some topsnps couldnt be matched to columns in the LD matrix. Matching failed for row(s): ", paste0(which(matchCheck), " ('",topsnps$topsnp[matchCheck],"')",collapse=","))
  }
  
  
  #Label the credible set to which each SNP belongs
  names(topSNPLD)[-1] <- sapply(names(topSNPLD)[-1],function(x){paste0(x," (",trait,":",topsnps$cs[topsnps$topsnp==x],")")})
  
  #Generate a list containing two datasets suitable for ggplot.
  #The first contains all snps in the region and second is a subset of these, including only those assigned to credible sets
  minimalSET<- list(allsnps=NULL,csonly=NULL)
  
  minimalSET$allsnps<- sets[c("snp","cs","pos")] %>%
    mutate(cs=gsub("L([0-9]+):.*",paste0(trait,":\\1"),cs)) %>%
    right_join(topSNPLD,by = "snp") %>%
    pivot_longer(cols=!all_of(c("snp","cs","pos")), names_to = "leadSNP",values_to = "corrs") %>%
    mutate(xval=ordered(pos,levels=sort(pos),labels=snp[order(pos)])) #Fix order to correspond with genomic position
  
  #Square the corrs if plotting R2 and declare the colourbar lower limit
  if(plotR2){
    minimalSET$allsnps$corrs <- minimalSET$allsnps$corrs^2
    min <- 0
    legendname <- ~r^2
  } else {
    min <- -1
    legendname <- ~r
  }
  
  minimalSET$csonly<- minimalSET$allsnps[!is.na(minimalSET$allsnps$cs),]
  
  #Generate LD heatmaps for both data subsets across the region relative to the top SNPs
  LDtopsnps<- lapply(minimalSET,function(x,breaks){
    
    ggplot(x,aes(y=leadSNP,x=xval,fill=corrs))+
      geom_tile()+
      scale_x_discrete(breaks=breaks)+ #On the x-axis include only positions for the top snps
      scale_fill_distiller(palette = heatmapPalette, guide="colourbar",
                           limit = c(min,1), direction=1,
                           na.value = "white",name=bquote(.(legendname))
                           )+
      coord_cartesian(expand=FALSE)+
      theme(axis.text.x = element_text(angle=315, vjust=.5, hjust=0))+
      labs(y= "Top SNP (Credible set)",x=paste0("GRCh",alignment," Genomic position (Chr",sets$chr[1],":",min(x$pos),"-",max(x$pos),")"))
    
  },breaks=topsnps$topsnp)
  
  #Generate a second figure matching the credible set only heatmap and generate a discrete colour scale above the plot
  csassign<- ggplot(minimalSET$csonly,aes(y=1,x=xval,fill=cs))+
    geom_tile()+
    {if(!is.null(discretePalette)) scale_fill_manual(values = discretePalette) }+
    labs(fill="Credible set")+
    coord_cartesian(expand=FALSE)+
    theme(axis.text = element_blank(), axis.ticks = element_blank(),axis.title = element_blank(),
          plot.margin = unit(c(0,0,0,0), "cm"))
  
  #With patchwork, combine the CS only heatmap with the assignments
  csOnlyHeatmap<- csassign/
    LDtopsnps$csonly+plot_layout(guides = 'collect',heights=c(1,10))  
  
  return(list(heatmap_allSNPs=LDtopsnps$allsnps,heatmap_CSonly=csOnlyHeatmap))
}