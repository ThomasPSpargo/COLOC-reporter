library(data.table)
library(ggplot2)
library(patchwork)
ggpalette = c("dodgerblue2", "green4", "#6A3D9A", "#FF7F00",
              "gold1", "skyblue2", "#FB9A99", "palegreen2", "#CAB2D6",
              "#FDBF6F", "gray70", "khaki2", "maroon", "orchid1", "deeppink1",
              "blue1", "steelblue4", "darkturquoise", "green1", "yellow4",
              "yellow3", "darkorange4", "brown")

tableS1<- data.table(openxlsx::read.xlsx("/Users/tom/Downloads/media-1 (5).xlsx",sheet="Table S1"))
tableS2<- data.table(openxlsx::read.xlsx("/Users/tom/Downloads/media-1 (5).xlsx",sheet="Table S2"))


## Create additional / recode columns for table S1 used in the script
alphabetical<- c("AD","ALS","FTD","PD","SZ")
tableS1[,`:=`(compar=factor(paste0(phen1," : ",phen2)),
              phen1=factor(phen1,levels=alphabetical),
              phen2=factor(phen2,levels=alphabetical),
              rho_direction=ifelse(rho>0,"pos","neg"),
              abs_rho=abs(rho),
              locusLabel=paste0(phen1,":", phen2," ",chr,":",start,"-",stop)
              )]

#Extract the more and less significant GWAS pval into their own column
tableS1[,`:=`(smallestP=apply(.SD,1,min),
              largestP=apply(.SD,1,max)),.SDcols=c("min_snp_pval_phen1","min_snp_pval_phen2")]


######
### Visualise correlation results overall as a volcano plot
######
volcano<- ggplot(tableS1,aes(x=rho,y=-log10(p),shape=rg_sigthresh,
                              # colour=fcase(rho>0 & p.fdr<0.05,"Positive",
                              #              rho<0 & p.fdr<0.05,"Negative")
                             colour=compar
                              ))+
  geom_point()+
  geom_hline(yintercept = -log10(max(tableS1[p.fdr<0.05,p]))-0.01,lty=2)+ # FDR line
  lims(x=c(-1,1))+
  # geom_errorbar(aes(xmin=ifelse(p.fdr<0.05,rho.lower,NA_real_),
  #                   xmax=ifelse(p.fdr<0.05,rho.upper,NA_real_)))+
  #scale_colour_discrete(na.value="black",breaks=c("Negative","Positive"),name="Direction of correlation in locus")+
  scale_colour_manual(values=ggpalette,breaks=levels(tableS1$compar),name="Trait pair")+
  scale_shape_manual(values=c("Bonf_sig"=15,"FDR_sig"=17),na.value = 19,
                     breaks=c("Bonf_sig","FDR_sig"),
                     labels=c(parse(text=paste0("P<",gsub("e", "%*%10^", signif(0.05/605,3)))),
                              bquote(~P[fdr]<0.05)),name="Significance")+
  labs(x=bquote(Local~genetic~correlation~(r[g])),
       y=bquote(-log[10]~p))+
  theme_minimal()+
  theme(legend.margin = margin(t=0.2,r=0.2,b=0.2,l=0.2,"mm"))

#ggsave("~/Downloads/lavaVolcano.pdf",device=cairo_pdf,width=150 ,height=125,units="mm")


#Ongoing checks will focus on the fdr significant correlations. Therefore, append the TableS1 columns of interest to the Table S2 data.
#Use non-equi join for start and stop to match the locus with the 10kb widened boundaries
S1keepCols<- c("rho","p","p.fdr","min_snp_pval_phen2","min_snp_pval_phen1","rho_direction","abs_rho","largestP","smallestP","locusLabel")
tableS2[tableS1,on=.(Trait_1=phen1,Trait_2=phen2,chr,start<=start,stop>=stop),c(S1keepCols):=mget(paste0("i.",S1keepCols))]

#Compute how many traits in pair have been assigned at least one credible set
tableS2[,nTraitsWithSet:=factor(
  fcase(Trait1_n_finemapping_credible_Sets>0 & Trait2_n_finemapping_credible_Sets>0,"Both",
        Trait1_n_finemapping_credible_Sets>0 | Trait2_n_finemapping_credible_Sets>0,"One",
        Trait1_n_finemapping_credible_Sets==0 & Trait2_n_finemapping_credible_Sets==0,"Neither"),
  levels=c("Neither","One","Both"))
  ]

######
### Compare the pvals
######

#Compute -log10(p) for largest and smallest P value
tableS2[,c("-log10_largestP","-log10_smallestP"):=lapply(.SD,function(x)-log10(x)),.SDcols=c("largestP","smallestP")]

metricplot<- lapply(c("-log10_largestP","-log10_smallestP","abs_rho"),function(x){
  #Generate y axis labels to be parsed
  # -log10(p) for A is the most sig snp from the trait with the MORE significant minimum p
  # -log10(p) for B is the most sig snp from the trait with the LESS significant minimum p
  ylab=setNames(nm=c("-log[10]~p[(A)]","-log[10]~p[(B)]","Absolute~r[g]"),c("-log10_smallestP","-log10_largestP","abs_rho"))
  
  compare_pval<- ggplot(tableS2[!is.na(rho_direction)],aes(y=get(x),x=rho_direction,fill=rho_direction))+#,x=locusLabel))+
    geom_boxplot(na.rm=TRUE,outlier.shape = NA)+
    geom_point(aes(group=locusLabel,shape=nTraitsWithSet),position=position_dodge(width=0.75))+
    scale_colour_manual(values=c("Negative"="red","Positive"="blue"),na.value = "black",
                        name="Direction of correlation in locus")+
    guides(fill="none")+
    scale_x_discrete(labels=c("neg"="Negative","pos"="Positive"))+
    scale_y_continuous(limits = c(0,NA),expand = expansion(mult=c(0.01,0.03)))+
    #ylim(c(0,NA))+
    labs(x="Direction of correlation in locus",y=eval(parse(text=names(ylab)[ylab==x])),
         shape="Traits with ≥1\ncredible set")+
    theme_bw()+
    theme(legend.position = "top",axis.title.x = element_blank(),
          panel.grid.major.x = element_blank())
  
  #Split by number of credible sets found in trait pair
  compare_pval_wfacets<- compare_pval+
    facet_grid(cols=vars(nTraitsWithSet),labeller = label_parsed,scales="free_y")+#,switch="x")+
    theme(axis.text.y=element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y=element_blank(),
          panel.spacing.x =unit(0.4,"mm"),
          #strip.background.x =element_blank(),strip.text.x = element_blank()
          #strip.placement = "outside",
          strip.background.y =element_blank(),strip.text.y = element_blank()
    )
  
  #Add faceting
  compare_pval_totFacets<- compare_pval+
    facet_grid(cols=vars("Total"),scales="free_y",labeller = label_parsed,switch="y")+
    theme(strip.background.y = element_blank(),
          #axis.title.y = element_blank(),
          strip.placement = "outside")
  #theme(strip.background.y =element_blank(),strip.text.y = element_blank())
  
  return(list(tot=compare_pval_totFacets,
              fullfacet=compare_pval_wfacets))
  
})

#Generate panel spacing
metricplot[[1]]$tot <- metricplot[[1]]$tot+theme(plot.margin = margin(t=0.5,r=0.1,b=0.1,l=0.1,unit = "mm"))
metricplot[[1]]$fullfacet <- metricplot[[1]]$fullfacet+theme(plot.margin = margin(t=0.5,r=0.1,b=0.1,l=0.1,unit = "mm"))

metricplot[[2]]$tot <- metricplot[[2]]$tot+theme(plot.margin = margin(t=0.1,r=0.1,b=0.1,l=0.1,unit = "mm"))
metricplot[[2]]$fullfacet <- metricplot[[2]]$fullfacet+theme(plot.margin = margin(t=0.1,r=0.1,b=0.1,l=0.1,unit = "mm"))
  
metricplot[[3]]$tot <- metricplot[[3]]$tot+theme(plot.margin = margin(t=0.1,r=0.1,b=0.5,l=0.1,unit = "mm"))
metricplot[[3]]$fullfacet <- metricplot[[3]]$fullfacet+theme(plot.margin = margin(t=0.1,r=0.1,b=0.5,l=0.1,unit = "mm"))

metricplot[1:2] <- lapply(metricplot[1:2],function(x){
  x$tot <- x$tot+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x =element_blank()
    )
  x$fullfacet <- x$fullfacet+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x =element_blank()
    )
  
  return(x)
})

metricplot[2:3] <- lapply(metricplot[2:3],function(x){
  x$tot <- x$tot+
    theme(strip.background =element_blank(),
          strip.text = element_blank()
    )
  x$fullfacet <- x$fullfacet+
    theme(strip.background =element_blank(),
          strip.text = element_blank()
    )
  
  return(x)
})

#Combine into into a patchwork
layout <- c(
  "ABBB
    CDDD
    EFFF
    GGGG"
)

pvalPlots <- metricplot[[1]]$tot+metricplot[[1]]$fullfacet+
  metricplot[[2]]$tot+metricplot[[2]]$fullfacet+
  metricplot[[3]]$tot+metricplot[[3]]$fullfacet+
  (ggplot(data.frame(),aes(x="test"))+theme_void()+labs(x="Direction of correlation in locus")+theme(axis.title.x = element_text()))+ #Create dummy plot to provide shared x-axis
  plot_layout(design=layout,guides="collect",heights=c(1,1,1,0.01)
              ) & 
  theme(legend.position = "top")
  #theme(legend.position = "right",legend.direction = "vertical")
  

volcBox<- free(volcano)/
  pvalPlots+
  plot_annotation(tag_levels = list(c("A","B"))) & theme(plot.tag = element_text(face = "bold"))# &
  #theme(legend.position)

ggsave(volcBox,filename = "~/Downloads/LAVAcorrelation_volcano_boxplots.pdf",device=cairo_pdf,width=150,height=210,units="mm")

#Save final figure
#ggsave("~/Downloads/LAVAcorrelation_finemap_comparison.pdf",device=cairo_pdf,width=150,height=125,units="mm")

