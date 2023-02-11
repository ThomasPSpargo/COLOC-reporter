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
#return_beta_comparisons
#sets - a data.frame containing at least the columns: "cs" - credible set,"beta" - an effect estimate for the direction of effect associated with the SNP
#topsnps - a data.frame containing at least the columns: "cs" - credible set,"topsnp" - id for the top SNP in the cs

susieCScorrs <- function(susie.fit,R,palette="PuOr",return_beta_comparisons=TRUE,sets=NULL,topsnps=NULL){
  
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
           corrs=if_else(rows==cols,NA_real_,corrs),
           corrlabs=if_else(rows>=cols,NA_real_,corrs)
           )
  
  #Determine the number of SNPs in each credible set; this will be plotted along the diagonal
  nsnps <- enframe(sapply(susie.fit$sets$cs,length)) %>%
    mutate(cs=as.numeric(gsub("L","",name),.keep="unused"))
  plot_prep$diag <- NA
  for(i in 1:nrow(plot_prep)){
   plot_prep[i,"diag"] <- ifelse(plot_prep[i,"rows"]==plot_prep[i,"cols"],nsnps$value[nsnps$cs %in% plot_prep[i,"rows"]],NA_real_)
  }
  
  if(return_beta_comparisons){
    
    if(is.null(sets) || is.null(topsnps)){
      stop("Please provide data to the sets and topsnps arguments in order to compare test statistics.")
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
    
    
    #Extract the beta for the top PIP SNP in each cs
    snpgrid<- sets[c("beta","snp")] %>%
      rename(topsnp=snp) %>%
      right_join(topsnps[c("cs","topsnp")],by = "topsnp")
    
    # Compare across every pair of CS, labelling '+' if effect estimate valence matches and '-' if they opposed
    # to be plotted on lower tri of subsequent heatmap
    cs_betacompare<- expand.grid(list(rows=snpgrid$cs,cols=snpgrid$cs)) %>%
      rowwise() %>%
      mutate(eff_dir=if_else(snpgrid$beta[snpgrid$cs==rows]>0 && snpgrid$beta[snpgrid$cs==cols]<0,"-","+")) %>%
      ungroup() %>%
      mutate(eff_dir=if_else(rows<=cols,"",eff_dir))
    
    
  } else {
    #Filter down to only upper tri elements
    plot_prep <- plot_prep %>%
      mutate(corrs=if_else(rows>=cols,NA_real_,corrs))
    
  }
  
  
  ## Return a heatmap of Credible set correlations
  #Plot rows and cols as factors for complete CS labelling, and sort y-axis ordering to present the heatmap like a table
  #The upper tri contains correlation coefficients between CS
  #The lower tri labels whether effect estimate direction for the top SNPs in each CS-pair is in the same or opposite direction
  cs_heatmap<- ggplot(plot_prep,aes(
    y=factor(rows,levels=rev(sort(snpgrid$cs))),
    x=factor(cols,levels=sort(snpgrid$cs))
  ))+
    geom_tile(aes(fill=corrs))+
    geom_text(aes(label=ifelse(!is.na(corrlabs),round(corrlabs,3),"")))+
    geom_text(aes(label=ifelse(!is.na(diag),diag,"")))+
    scale_fill_distiller(palette = palette, guide="colourbar",
                         limit = c(-1,1),
                         na.value = "white",name="Correlation")+
    coord_cartesian(expand=FALSE)+
    theme(axis.title = element_blank()
          #axis.ticks = element_blank(),axis.text = element_blank()
    ) +
    {if(return_beta_comparisons) {
      list(geom_text(data=cs_betacompare,aes(label=eff_dir)),
           labs(caption=str_wrap("The diagonal specifies the number of SNPs in a given credible set. The upper tri displays the visualised correlation estimate and the lower triangle indicates if the effect estimates for snps with the top PIP across two CS have the same valence (+) or not (-).",width = 100)),
           theme(plot.caption = element_text(hjust = 0)))
    } else {
           theme(plot.caption = element_blank(),
                 axis.ticks = element_blank(),axis.text = element_blank(),
                 legend.position = c(0.1,0.5)
                 )
    }
    }
  
  
  return(list(heatmap=cs_heatmap,corrs=cs_correlations))
}