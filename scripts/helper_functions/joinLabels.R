#####
# Author: Thomas Spargo (thomas.spargo@kcl.ac.uk)
# Obtained from GitHub repository: https://github.com/ThomasPSpargo/COLOC-reporter
#
# This function generates a character vector according to the distribution of variants across credible sets from analysis of two traits
#####

# x and y are dataframes with 'snp' column used as an ID and 'cs' column indicating credible set assignments as a factor. nrow(x) should equal nrow(y)
# traits is character vector of length 2 which identifies each of the traits 'x' and 'y'
joinLabels <- function(x, y,traits){
  
  # z is a factor class vector, trait is the associated trait
  refactor <- function(z,trait){
    #Reformat factor levels if not all NA
    len<- length(levels(z))
    if(len>0){levels(z) <- c(paste0(trait,":",1:len))}
    z <- addNA(z)                      #Add a NA factor level
    levels(z)[length(levels(z))] <- "" #Set NA level as empty string
    z
  }
  
  #Refactor both x and y traits
  x$cs<- refactor(x$cs,traits[1])
  y$cs<- refactor(y$cs,traits[2])
  
  #Combine into a snp and credible set data frame
  xy <- full_join(x,y,by="snp") %>%
    mutate(cs=case_when(cs.x == "" & cs.y == "" ~ NA_character_,
                        cs.x != "" & cs.y == "" ~ paste0(cs.x),
                        cs.x == "" & cs.y != "" ~ paste0(cs.y),
                        cs.x != "" & cs.y != "" ~ paste0(cs.x," & ",cs.y)),
           cs=as.factor(cs)
    ) %>%
    dplyr::select(snp,cs)
  
  xy
}
