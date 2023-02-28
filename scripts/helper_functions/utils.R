#####
# Author: Thomas Spargo (thomas.spargo@kcl.ac.uk)
# Obtained from GitHub repository: https://github.com/ThomasPSpargo/COLOC-reporter
#
# This script contains a series of small helper functions which do not warrant their own script individually but clutter the colocalisation workflow otherwise
#####

#Function to scale x-axis based on plot range
xScaler <- function(range){
  #Rescale x-axis into MBs or KBs
  if(range>100000){
    xscale <- 1000000
    xscale_magnitude <- " [Mb]"
  } else if (range>100) {
    xscale <- 1000
    xscale_magnitude <- " [Kb]"
  } else {
    xscale <- 1
    xscale_magnitude <- ""
  }
  return(list(xscale=xscale,xscale_magnitude=xscale_magnitude))
}


# #Small wrapper function for automating assignment to global environment list variables
# #Used for the purpose of passing lists which contain params for Rmd documents
# # x - a string referring to a list object in the global environment and to which which a new element will be appended
# # value - the new element that should be appended to the list
# multi - set to TRUE to supply an existing list unaltered
assignToList <- function(x,value,name=NULL,envir = .GlobalEnv,multi=FALSE){
  
  if(multi){
    newElem <- value
  } else {
    #Create the new list element
    newElem <- list(value)
    if(!is.null(name)) names(newElem) <- name
  }
  if(!x %in% ls(envir=envir)){
    #Assign to a NEW list
    assign(x, value=newElem, envir = envir)
    
  } else {
    #Append to the existing list, called via get
    assign(x, value=c(get(x,envir=envir),newElem)
           , envir = envir)
  }
}

#Simple matching function to obtain genes within a chromosome nearby to a base pair window
# set highConfidence to TRUE to only return genes which encompass the full BP range
getNearbyGenes <- function(bp_range,genesInChr,highConfidence=FALSE){
  if(!highConfidence){
    nearbyGenes<-  genesInChr[with(genesInChr,
                                   (start_window > min(bp_range) & start_window < max(bp_range)) | # Test if start position is within range
                                     (end_window > min(bp_range) & end_window < max(bp_range)) |   # or if end position is within range
                                     (start_window < min(bp_range) & end_window > max(bp_range))   # or if start and end positions are both outside range [gene straddles window]
    ),]
  } else {
    nearbyGenes<-  genesInChr[genesInChr$start_window < min(bp_range) & genesInChr$end_window > max(bp_range),]   # Only return instances where the gene straddles the window
  }
  return(nearbyGenes)
}

#This is a wrapper around render to simplify report writing
renderReport<- function(params,outfile,template){
  rmarkdown::render(template,
                    output_file = outfile,
                    quiet=TRUE,
                    params = get(params,envir=.GlobalEnv),
                    envir = new.env(),
                    intermediates_dir = tempdir())
}
