#####
# Author: Thomas Spargo (thomas.spargo@kcl.ac.uk)
# Obtained from GitHub repository: https://github.com/ThomasPSpargo/COLOC-reporter
#
# This script generates a heatmap to visualise correlations between credible sets identified by susie estimated via the get_cs_correlation() function.
#####

#Inputs
#susie.fit the results of a susie finemapping analysis
#R - an LD matrix associated with susie.fit
#palette - the palette to use on the ggplot
#inc_Z - Indicate TRUE to include Z-scores for top PIP snps in the plot This does requireme 
#sets - a data.frame containing at least the column: "cs" - credible set along with either "z" or"beta" and "SE" (for deriving Z)
#topsnps - a data.frame containing at least the columns: "cs" - credible set,"topsnp" - id for the top SNP in the cs

susieCScorrs <- function(susie.fit,R,palette="PuOr",inc_Z=TRUE,sets=NULL,topsnps=NULL){
  
  #Obtain correlations between credible sets
  cs_correlations<- get_cs_correlation(susie.fit,Xcorr = R)
  ## Note that: cs_correlations^2 is equivalent to cs LD r^2  [which is calculated internally within coloc::runsusie wrapper function]
  
  # Format for heatmap plotting
  # 'corrlabs' will be plotted on the upper tri of the heatmap and the diagonal will either be blank or contain the CS names
  plot_prep<- as.data.frame(cs_correlations) %>%
    rownames_to_column("rows") %>%
    pivot_longer(cols=!rows, names_to = "cols",values_to = "corrs") %>%
    mutate(rows=as.numeric(gsub("L","",rows)),
           cols=as.numeric(gsub("L","",cols)),
           corrs=if_else(rows<=cols,NA_real_,corrs), #Change these signs to alter top and bottom tri plotting
           corrlabs=if_else(rows<=cols,NA_real_,corrs)
    )
  
  #Determine the number of SNPs in each credible set; this will be plotted along the diagonal
  nsnps <- enframe(sapply(susie.fit$sets$cs,length)) %>%
    mutate(cs=as.numeric(gsub("L","",name)))
  
  #Prepare a column to store values which will go along the diagonal
  plot_prep$diag <- NA
  
  
  if(!inc_Z){
    
    #Filter down to only upper tri elements
    plot_prep <- plot_prep %>%
      mutate(corrs=if_else(rows>=cols,NA_real_,corrs))
    
    #Obtain N for the diagonal
    for(i in 1:nrow(plot_prep)){
      plot_prep[i,"diag"] <- ifelse(plot_prep[i,"rows"]==plot_prep[i,"cols"],paste0("N = ",nsnps$value[nsnps$cs %in% plot_prep[i,"rows"]]),NA_character_)
    }
    
  } else {
    
    if(is.null(sets) || is.null(topsnps)){
      stop("Please provide data to the sets and topsnps arguments in order to obtain the Z statistic.")
    }
    
    #Filter step to ensure an index of only one top snp is used per CS
    if(any(tabulate(topsnps$cs)>1)){
      message("Randomly downsampling credible sets with >",1," top snps provided")
      
      incSets<- unique(topsnps$cs)
      for(i in seq_along(incSets)){
        logiIndex<- topsnps$cs==incSets[i] #Get a logical index of snps in and outside the current CS
        
        #If the cs includes snps above the threshold, filter down
        if(sum(logiIndex)>1){
          topsnps<- rbind(topsnps[!logiIndex,],
                          topsnps[sample(which(logiIndex),size=1),])
        }
      }
    }
    
    #Subset sets to the columns of interest
    cols <- c("beta","SE","z","snp") 
    cols<- cols[cols %in% colnames(sets)]
    
    #Extract the beta for the top PIP SNP in each cs
    snpgrid<- sets[cols] %>%
      rename(topsnp=snp) %>%
      right_join(topsnps[c("cs","topsnp")],by = "topsnp")
    
    #If not provided, calculate the z column
    if(!"z" %in% colnames(snpgrid)){ 
      if(!all(c('beta',"SE") %in% colnames(snpgrid))) stop("Please provide either the column 'z' or the columns 'beta' and 'SE' from which 'z' can be derived.")
      snpgrid$z<- snpgrid$beta/snpgrid$SE 
    }
  
    tsnps<-full_join(snpgrid,nsnps,by="cs")
    for(i in 1:nrow(plot_prep)){
      index <- tsnps$cs %in% plot_prep[i,"rows"]
      plot_prep[i,"diag"] <- ifelse(plot_prep[i,"rows"]==plot_prep[i,"cols"],
                                    paste0("N = ",tsnps$value[index],",\nZ = ",signif(tsnps$z[index],3)),NA_character_)
    }
    
    
    # # Obtain Z for every pair of CS
    # cs_betacompare<- expand.grid(list(rows=snpgrid$cs,cols=snpgrid$cs)) %>%
    #   rowwise() %>%
    #   #mutate(stat= paste0("frac(",signif(snpgrid$z[snpgrid$cs==cols],3),",",signif(snpgrid$z[snpgrid$cs==rows],3),")")) %>%
    #   ungroup() %>%
    #   mutate(stat=if_else(rows<=cols,"",stat))
    
    
  }
  
  
  ## Return a heatmap of Credible set correlations
  #Plot rows and cols as factors for complete CS labelling, and sort y-axis ordering to present the heatmap like a table
  cs_heatmap<- ggplot(plot_prep,aes(
    y=factor(rows,levels=rev(sort(snpgrid$cs))),
    x=factor(cols,levels=sort(snpgrid$cs))
  ))+
    geom_tile(aes(fill=corrs))+
    geom_text(aes(label=round(corrlabs,3)),na.rm=TRUE)+
    geom_text(aes(label=diag),na.rm = TRUE)+
    scale_fill_distiller(palette = palette, guide="colourbar",
                         limit = c(-1,1),
                         na.value = "white",name="Correlation")+
    coord_cartesian(expand=FALSE)+
    theme(axis.title = element_blank()
          #axis.ticks = element_blank(),axis.text = element_blank()
    )# +
    # {if(inc_Z) { geom_text(data=cs_betacompare,aes(label=stat),parse=TRUE)
    # } else {
    # theme(plot.caption = element_blank(),
    #       axis.ticks = element_blank(),axis.text = element_blank(),
    #       legend.position = c(0.1,0.5)
    #       )
    #} }
  
  
  return(list(heatmap=cs_heatmap,corrs=cs_correlations))
}
