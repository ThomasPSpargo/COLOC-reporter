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
#topsnps - optionally, a dataframe declaring the top snp for each credible set in the column "topsnp" - used as the reference point in the heatmap. Expects the columns "cs" (the credible set label) and "topsnp" (the snp to which this label relates)
#susie.fit - optionally, and if not supplying information to 'topsnps', this is either a model output by susie which contains at least one credible set or a list of SuSiE models, across which at least one credible set occurs. The snps with the top PIP from each credible set will be extracted from each model.
#sets - a dataframe containing information about the credible sets identified by susie and their associated SNPs. Expects the columns "snp" - ID, corresponding to the topsnps,"cs" - credible set assignments per SNP, and "pos" - genomic position
#R - an LD correlation matrix for all SNPs supplied to susie.

#heatmapPalette="OrRd" - specify the colour palette to use for the heatmap
#discretePalette=NULL  - specify the colour palette to use for visualising different credible sets in a bar marking CS assignments relevant to each snp
#plotR2=TRUE           - TRUE returns an R2 plot, FALSE returns correlations
#trait                 - optionally, trait string id for the finemapped trait. This is ecc If supplying a list of models to susie.fit, a vector of model trait IDs must be supplied in the order to which the list corresponds.
#maxPerCS              - indicate a numeric limit for the number of 'top' snps that will be plotted along the y-axis for each credible set. This is included to avoid overplotting.
#chr                   - optionally indicate the chromosome 
#separateLegends=FALSE - logical, defaulting to FALSE. Set TRUE indicate that patchwork should NOT collect legends across the colourbar and heatmap plot. This is required for any further collection of patchwork the figures into higher level patchworks


susie_cs_ld <- function(topsnps=NULL,susie.fit=NULL,R,sets,heatmapPalette="OrRd",discretePalette=NULL,plotR2=TRUE,trait=NULL,maxPerCS=10,chr=NULL,separateLegends=FALSE){
  
  ## If top snps not supplied directly, determine top PIP snps from the SuSiE fit supplied
  if(is.null(topsnps)){ 
    
    ## Internal function for extracting top snps from each cs (return an empty dataframe if none are present)
    extractTopPIPs<- function(susie.fit,trait=NULL){
      sum_susie <- summary(susie.fit)
      if(!is.null(sum_susie$cs)){
        topsnp <- data.frame(cs=as.character(sum_susie$cs$cs),topsnp=rep(NA_real_,length(sum_susie$cs$cs)))
        if(!is.null(trait)) topsnp$cs <- paste0(trait,":",topsnp$cs)
        for(j in 1:nrow(sum_susie$cs)){
          snp_index <- as.numeric(unlist(strsplit(sum_susie$cs$variable[j], ',')))
          topsnp$topsnp[j]<- paste(names(susie.fit$pip[snp_index])[susie.fit$pip[snp_index]==max(susie.fit$pip[snp_index])],collapse=", ")
        }
      } else {
        topsnp <- data.frame(cs=character(0),topsnp=numeric(0))
      }
      return(topsnp)
    }
    
    if(class(susie.fit)=="list" && !is.null(trait)){   #Concatenate across multiple SuSiE fits
      topsnps<- mapply(extractTopPIPs,susie.fit,trait,SIMPLIFY=FALSE)
      topsnps <- do.call(rbind,topsnps)
      
    } else if (class(susie.fit)=="list" && is.null(trait)){ #Only allow list imports when trait IDs are specified
      stop("Please supply trait IDs required for labelling which top SNPs are associated with each credible set. This should be provided to the 'trait' argument in the same order as the list of susie fits.")
      
    } else { #Import from a single susie fit
      topsnps <- extractTopPIPs(susie.fit,trait)
    }
    
    #Split into long-format any instances where multiple 'top' snps are matched for a cs
    topsnps <- separate_rows(topsnps, all_of("topsnp"), sep = ", ")
  }
  
  #Check for any top SNPs tagged for more than one CS, if present, concatenate their credible set assignments and deduplicate the dataset
  tab<-table(topsnps$topsnp)
  if(any(tab>1)){
    dups <- names(tab)[tab>1]
    for(i in seq_along(dups)){ topsnps$cs[topsnps$topsnp==dups[i]] <- paste0(topsnps$cs[topsnps$topsnp==dups[i]],collapse = ", ") }
    topsnps <- unique.data.frame(topsnps)
  }
  
  ## Filter step to reduce overplotting if a large number of snps are identified as 'top'
  if(any(table(topsnps$cs)>maxPerCS)){ message("Randomly downsampling credible sets with >",maxPerCS," top snps identified")
    
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
  
  #Warn if any snps aren't matched
  matchCheck <- !topsnps$topsnp %in% colnames(topSNPLD)
  if(any(matchCheck)) warning("Some topsnps couldnt be matched to columns in the LD matrix. Matching failed for row(s): ", paste0(which(matchCheck), " ('",topsnps$topsnp[matchCheck],"')",collapse=","))
  
  #Label the credible set to which each top SNP belongs [the columns of the list]
  names(topSNPLD)[-1] <- sapply(names(topSNPLD)[-1],function(x){paste0(x," (",topsnps$cs[topsnps$topsnp==x],")")})
  
  #Generate a list containing two datasets suitable for ggplot.
  #The first contains all snps in the region and second is a subset of these, including only those assigned to credible sets
  minimalSET<- list(allsnps=NULL,csonly=NULL)
  
  minimalSET$allsnps<- sets[c("snp","cs","pos")] %>%
    left_join(topSNPLD,by = "snp") %>%
    pivot_longer(cols=!all_of(c("snp","cs","pos")), names_to = "leadSNP",values_to = "corrs") %>%
    arrange(pos) %>% #Fix order to correspond with genomic position
    mutate(xval=ordered(pos,levels=pos,labels=snp))
  
  #Square the corrs if plotting R2 and declare the colourbar lower limit
  if(plotR2){
    minimalSET$allsnps$corrs <- minimalSET$allsnps$corrs^2
    min <- 0
    legendname <- ~r^2
  } else {
    min <- -1
    legendname <- ~r
  }
  
  #Create a data subset with only the SNPs assigned to credible sets
  minimalSET$csonly<- minimalSET$allsnps[!is.na(minimalSET$allsnps$cs),]
  
  #Generate LD heatmaps for both data subsets across the region relative to the top SNPs
  LDtopsnps<- lapply(minimalSET,function(x,breaks){
    
    ggplot(x,aes(y=leadSNP,x=xval,fill=corrs,colour=corrs))+
      geom_tile()+
      scale_x_discrete(breaks=breaks,
                       guide = guide_axis(n.dodge=2)
                       )+ #On the x-axis include only positions for the top snps
      scale_fill_distiller(palette = heatmapPalette, guide="colourbar",
                           limit = c(min,1), direction=1,
                           na.value = "white",name=bquote(.(legendname))
                           )+
      scale_colour_distiller(palette = heatmapPalette, guide="none",
                           limit = c(min,1), direction=1,
                           na.value = "white")+
      coord_cartesian(expand=FALSE)+
      theme(panel.grid = element_blank())+
      labs(y= "Top SNP (Credible set)", x=ifelse(!is.null(chr),
                    paste0("Genomic position (Chr",chr,":",min(x$pos),"-",max(x$pos),")"),
                    paste0("Position (", min(x$pos),"-",max(x$pos),")")
                    )
           )
        
    
  },breaks=topsnps$topsnp)
  
  csassign<- lapply(minimalSET,function(x){
    #Generate a second figure matching the dataset displaying a discrete colour bar to mark credible sets above the plot
    ggplot(x,aes(y=1,x=xval,fill=cs,colour=cs))+
      geom_tile()+
      {if(!is.null(discretePalette)) list(scale_fill_manual(values = discretePalette,na.value = "white",breaks=levels(factor(minimalSET$allsnps$cs))),
                                          scale_colour_manual(values = discretePalette,na.value = "white",breaks=levels(factor(minimalSET$allsnps$cs)),guide="none")
                                          ) }+
      labs(fill="Credible set")+
      coord_cartesian(expand=FALSE)+
      theme(axis.text = element_blank(),
            axis.ticks = element_blank(),axis.title = element_blank(),
            plot.margin = unit(c(0,0,0,0), "cm"))
    
  })

  #With patchwork, combine the heatmap(s) with the CS assignments
  heatmapWbar<-  mapply(function(heatmap,csBar){ csBar/heatmap+plot_layout(guides=switch(!separateLegends,'collect',NULL),
                                                                           heights=c(1,10),ncol=1,nrow=2) },
                        heatmap=LDtopsnps,csBar=csassign,SIMPLIFY=FALSE)
  
  return(list(heatmap_allSNPs=heatmapWbar$allsnps,heatmap_CSonly=heatmapWbar$csonly))
}
