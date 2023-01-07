#####
# Author: Thomas Spargo (thomas.spargo@kcl.ac.uk)
# Obtained from GitHub repository: https://github.com/ThomasPSpargo/COLOC-reporter
#
# Generate plots in a column with a shared legend, using custom function which wraps around egg::ggarrange and draws upon grid functionality
#####

#x - a list of ggplots to be stacked vertically and assigned a shared figure legend
#position - place legend at either "right" or "bottom"
# ... - arguments passed to egg::ggarrange

egg_ggarr_wlegend <- function(p,position="right",...){
  #Syntax developed based on
  #https://github.com/tidyverse/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
  
  # #Extract legend from first figure as a grob, then determine parameters
  g <- ggplotGrob(p[[1]] + theme(legend.position = position))$grobs
  is.legend <- which(sapply(g, function(x) x$name) == "guide-box")
  if(length(is.legend)>0){
    legend <- g[[is.legend]]
    lheight <- sum(legend$height)
    lwidth <- sum(legend$width)
    #Drop legend from individual plots
    p <-  lapply(p, function(x) x + theme(legend.position="none"))
  }
  
  
  #Always drop duplicated x-axis labeling elements, since columns are always stacked vertically
  p[1:(length(p)-1)] <-  lapply(p[1:(length(p)-1)], function(x) x + theme(axis.title.x = element_blank(),
                                                                          axis.line.x = element_blank(),
                                                                          axis.text.x = element_blank())
  )
  
  #Using egg_ggarrange, combine the plots with legend removed; egg::ggarrange ensures combined plots are well formatted
  eggP <- egg::ggarrange(plots=p,...)
  
  if(length(is.legend)>0){
    #Then, using grid functionality (which should be installed already, as a dependency of egg)
    eggP <- switch(position,
                   "bottom" = arrangeGrob(eggP,
                                          legend,
                                          ncol = 1,
                                          heights = grid::unit.c(grid::unit(1, "npc") - lheight, lheight)
                   ),
                   "right" = arrangeGrob(eggP,
                                         legend,
                                         ncol = 2,
                                         widths = grid::unit.c(grid::unit(1, "npc") - lwidth, lwidth)
                   ))
    
    grid::grid.newpage()
    grid::grid.draw(eggP)
 
    # return gtable invisibly
    invisible(eggP)
  }
  return(eggP)
}