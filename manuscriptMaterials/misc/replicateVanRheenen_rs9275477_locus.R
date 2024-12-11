library(coloc)
library(data.table)
### Replicate van Rheenen colocalisation analysis in AD:ALS locus around rs9275477

#Declare the region for analysis
reg_range <-c(chr=6,start=32572641,stop=32772641)


#Read in the summary stats,filtering to the region of interest.
#'Clean' versions of sumstats reflect data after filtering variants in the 200kb locus defined by van Rheenen using the protocol followed elswhere in the current study.
# Versions without the sumstat filtering protocol are 'raw' in the sense that we perform colocalisation analysis across the full set of SNPs provided for this locus across the original studies.

sumstatFile<- list.files("/Users/tom/Downloads/vanRheenenreplicate",full.names=TRUE)
# [1] "/Users/tom/Downloads/vanRheenenreplicate/AD_sumstats_Kunkle.clean.pt05.gz"     
# [2] "/Users/tom/Downloads/vanRheenenreplicate/AD_sumstats_Kunkle.txt"               
# [3] "/Users/tom/Downloads/vanRheenenreplicate/ALS_sumstats_vanRheenen.clean.pt05.gz"
# [4] "/Users/tom/Downloads/vanRheenenreplicate/ALS_sumstats_vanRheenen.tsv" 

traits<- lapply(sumstatFile,fread)
for(i in seq_along(traits)) traits[[i]] <- traits[[i]][CHR==reg_range["chr"] & BP>=reg_range["start"] & BP<=reg_range["stop"]]
names(traits) <- paste0(gsub("(AD|ALS).*","\\1",basename(sumstatFile)),ifelse(grepl("clean",basename(sumstatFile)),"","_raw"))


#Restrict to shared SNPs only
traits$ALS_raw <- traits$ALS_raw[SNP %in% traits$AD_raw$SNP]
traits$AD_raw <- traits$AD_raw[SNP %in% traits$ALS_raw$SNP]

traits$ALS <- traits$ALS[SNP %in% traits$AD$SNP][order(BP)]
traits$AD <- traits$AD[SNP %in% traits$ALS$SNP][order(BP)]

#sapply(traits,nrow) - 1563 snps pass our filtering protocol. 2342 were analysed originally
#   AD  AD_raw     ALS ALS_raw 
# 1563    2342    1563    2342 


#Append allele freqs from the ALS dataset for the 'raw' dataset analysis (these are missing from the AD data)
traits$AD[traits$ALS,on=.(SNP),FREQ:=i.FREQ]
traits$AD_raw[traits$ALS_raw,on=.(SNP),FREQ:=i.FREQ]

#Vidsally compare SNPS
compareP <- traits$ALS[,.(SNP,P,FREQ)][traits$AD[,.(SNP,P)],on=.(SNP)]
compareP[,i.P:=as.numeric(i.P)]
ggplot(compareP,aes(x=-log10(i.P),y=-log10(P),colour=FREQ<=0.4 | FREQ>=0.6))+
  geom_point()

#setequal(traits$ALS$SNP,test$SNP)

#Prepare columns used by coloc
traits <- lapply(traits,function(x){
  
  x <- x[!is.na(SNP)]
  
  cols<- c("pvalues"="P","snp"="SNP","beta"="BETA")
  setnames(x,cols,names(cols))
  x$varbeta <-x$SE^2
  x$A1 <-toupper(x$A1)
  x$A2 <-toupper(x$A2)
  
  x$MAF <- NA
  x$MAF[x$FREQ<=0.5] <- x$FREQ[x$FREQ<=0.5]
  x$MAF[x$FREQ>0.5] <- 1-x$FREQ[x$FREQ>0.5]
  
  x<- as.list(x)
  x$type="cc"
  x$s <- 0.5
  
  return(x)
})

#Run coloc.abf and write results:

coloc.abf(traits$ALS,traits$AD)
# nsnps           H0           H1           H2           H3           H4 
# 1.563000e+03 5.148664e-09 1.719301e-03 1.614886e-07 5.298081e-02 9.452997e-01 

coloc.abf(traits$ALS_raw,traits$AD_raw)
# nsnps           H0           H1           H2           H3           H4 
# 2.342000e+03 2.448401e-09 8.915116e-04 1.355858e-06 4.931890e-01 5.059182e-01 


