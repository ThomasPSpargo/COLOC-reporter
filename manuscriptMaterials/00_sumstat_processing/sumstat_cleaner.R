#!/usr/bin/Rscript
# This script was written by Oliver Pain whilst at King's College London University.
start.time <- Sys.time()
suppressMessages(library("optparse"))

option_list = list(
  make_option("--sumstats", action="store", default=NA, type='character',
              help="Path to summary statistics file [required]"),
  make_option("--ref_plink_chr", action="store", default=NA, type='character',
              help="Path to per chromosome PLINK files [required]"),
  make_option("--ref_freq_chr", action="store", default=NA, type='character',
              help="Path to per chromosome PLINK frequency (.frq) files [optional]"),
  make_option("--info", action="store", default=0.9, type='numeric',
              help="INFO threshold [optional]"),
  make_option("--maf", action="store", default=0.01, type='numeric',
              help="MAF threshold [optional]"),
  make_option("--exclude_region", action="store", default=NULL, type='character',
              help="Plain text file with the header CHR, START, and STOP. File rows indicate genomic regions to exclude [optional]"),
  make_option("--maf_diff", action="store", default=0.2, type='numeric',
              help="Difference between reference and reported MAF threshold [optional]"),
  make_option("--insert_ref_maf", action="store", default=T, type='logical',
              help="Set to T to insert reference allele frequency [optional]"),
  make_option("--SDcheck", action="store", default="none", type='character',
              help="Specify 'quant' to perform check of whether SD is reasonable for all snps in a quantitative trait, or 'binary' to check for a binary trait. see doi: 10.1093/bioinformatics/btaa1029, section 3.4"),
  make_option("--gz", action="store", default=T, type='logical',
              help="Set to T to gzip summary statistics [optional]"),
  make_option("--output", action="store", default='./Output', type='character',
              help="Path for output files [optional]")
)

opt = parse_args(OptionParser(option_list=option_list))

library(data.table)

opt$output_dir<-paste0(dirname(opt$output),'/')
system(paste0('mkdir -p ',opt$output_dir))

sink(file = paste(opt$output,'.log',sep=''), append = F)
cat(
  '#################################################################
# sumstat_cleaner.R
# For questions contact Oliver Pain (oliver.pain@kcl.ac.uk)
#################################################################
Analysis started at',as.character(start.time),'
Options are:\n')

cat('Options are:\n')
print(opt)
cat('Analysis started at',as.character(start.time),'\n')
sink()

#####
# Read in sumstats
#####

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Reading in GWAS sumstats.\n')
sink()

if(substr(opt$sumstats,(nchar(opt$sumstats)+1)-3,nchar(opt$sumstats)) == '.gz'){
  GWAS<-data.frame(fread(cmd=paste0('zcat ',opt$sumstats)))
} else {	
  GWAS<-data.frame(fread(opt$sumstats))
}

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('GWAS contains',dim(GWAS)[1],'variants.\n')
sink()

# Restrict sumstats to optional or required columns
GWAS<-GWAS[,(names(GWAS) %in% c('SNP','A1','A2','P','OR','BETA','SE','N','FREQ','INFO'))]

#####
# Read in ref_bim
#####

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Reading in ref_plink_chr\n')
sink()

#Not run alternative syntax for reading in all autosomes bim file
# #If the string 'chr' is detected, then expect per-chromosome reference files, otherwise, read from a single file
# if(grepl("chr$",basename(opt$ref_plink_chr))){
ref_bim<-NULL
for(i in 1:22){
  bim<-fread(paste0(opt$ref_plink_chr,i,'.bim'))
  ref_bim<-rbind(ref_bim, bim)
}
# } else {
#   ref_bim<-fread(opt$ref_plink_chr)
# }

names(ref_bim)<-c('CHR','SNP','POS','BP','A1','A2')
ref_bim$POS<-NULL

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('ref_plink_chr contains',dim(ref_bim)[1],'variants.\n')
sink()


#####
# Remove SNPs that are not in ref_bim
#####

GWAS<-GWAS[(GWAS$SNP %in% ref_bim$SNP),]

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('After removal of variants that are not in ref_bim,',dim(GWAS)[1],'variants remain.\n')
sink()

#####
# Remove SNPs within genomic regions indicated by file pointed to within --exclude_region option
#####
if(!is.null(opt$exclude_region)){
  
  exclregions <- fread(opt$exclude_region)
  
  for(i in 1:nrow(exclregions)){
    
    sink(file = paste(opt$output,'.log',sep=''), append = T)
    cat(paste0("Excluding variants within the region chr",exclregions[i,"CHR"],":",exclregions[i,"START"],"-",exclregions[i,"STOP"],"...\n"))
    sink()
    
    #Identify SNPs in region from reference dataset
    ref_bim_excl <- ref_bim[ref_bim$BP >= exclregions[i,"START"] & ref_bim$BP <= exclregions[i,"STOP"] & ref_bim$CHR == exclregions[i,"CHR"],]
    
    #Identify any GWAS SNPs in the exclusion region and remove
    region_overlap <- which(GWAS$SNP %in% ref_bim_excl$SNP)
    if(length(region_overlap)>0) GWAS <- GWAS[-region_overlap,]
    
    #Print number dropped
    sink(file = paste(opt$output,'.log',sep=''), append = T)
    cat(length(region_overlap),'variants this region have been dropped,',dim(GWAS)[1],'variants remain.\n')
    sink()
  }
}

#####
# Remove variants that are not SNPs or were strand-ambiguous.
#####

GWAS<-GWAS[nchar(GWAS$A1) == 1 & nchar(GWAS$A2) == 1,]

#Convert SNP character strings for A1 and A2 into upper case
GWAS$A1 <- toupper(GWAS$A1)
GWAS$A2 <- toupper(GWAS$A2)

#IUPAC codes: https://www.snp-nexus.org/v4/guide/
GWAS$IUPAC[GWAS$A1 == 'A' & GWAS$A2 =='T' | GWAS$A1 == 'T' & GWAS$A2 =='A']<-'W'
GWAS$IUPAC[GWAS$A1 == 'C' & GWAS$A2 =='G' | GWAS$A1 == 'G' & GWAS$A2 =='C']<-'S'
GWAS$IUPAC[GWAS$A1 == 'A' & GWAS$A2 =='G' | GWAS$A1 == 'G' & GWAS$A2 =='A']<-'R'
GWAS$IUPAC[GWAS$A1 == 'C' & GWAS$A2 =='T' | GWAS$A1 == 'T' & GWAS$A2 =='C']<-'Y'
GWAS$IUPAC[GWAS$A1 == 'G' & GWAS$A2 =='T' | GWAS$A1 == 'T' & GWAS$A2 =='G']<-'K'
GWAS$IUPAC[GWAS$A1 == 'A' & GWAS$A2 =='C' | GWAS$A1 == 'C' & GWAS$A2 =='A']<-'M'

GWAS<-GWAS[(GWAS$IUPAC %in% c('R', 'Y', 'K', 'M')),]

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('After removal of variants that are not SNPs or are ambiguous,',dim(GWAS)[1],'variants remain.\n')
sink()

#####
# Check allele match reference and insert CHR & BP from reference
#####

if(sum(names(GWAS) == 'CHR') == 1){
  GWAS$CHR<-NULL
}

if(sum(names(GWAS) == 'ORIGBP') == 1){
  GWAS$ORIGBP<-NULL
}

ref_bim$IUPAC[ref_bim$A1 == 'A' & ref_bim$A2 =='G' | ref_bim$A1 == 'G' & ref_bim$A2 =='A']<-'R'
ref_bim$IUPAC[ref_bim$A1 == 'C' & ref_bim$A2 =='T' | ref_bim$A1 == 'T' & ref_bim$A2 =='C']<-'Y'
ref_bim$IUPAC[ref_bim$A1 == 'G' & ref_bim$A2 =='T' | ref_bim$A1 == 'T' & ref_bim$A2 =='G']<-'K'
ref_bim$IUPAC[ref_bim$A1 == 'A' & ref_bim$A2 =='C' | ref_bim$A1 == 'C' & ref_bim$A2 =='A']<-'M'

ref_bim_GWAS<-merge(ref_bim,GWAS, by='SNP')

ref_bim_GWAS<-ref_bim_GWAS[ref_bim_GWAS$IUPAC.x == ref_bim_GWAS$IUPAC.y | ref_bim_GWAS$IUPAC.x == 'R' & ref_bim_GWAS$IUPAC.y == 'Y' | ref_bim_GWAS$IUPAC.x == 'Y' & ref_bim_GWAS$IUPAC.y == 'R' | ref_bim_GWAS$IUPAC.x == 'K' & ref_bim_GWAS$IUPAC.y == 'M' | ref_bim_GWAS$IUPAC.x == 'M' & ref_bim_GWAS$IUPAC.y == 'K',]

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('After removal of SNPs with alleles that do not match ref_bim,',dim(ref_bim_GWAS)[1],'variants remain.\n')
sink()

#####
# Flip strand to match reference
#####

GWAS_noflip<-ref_bim_GWAS[ref_bim_GWAS$IUPAC.x == ref_bim_GWAS$IUPAC.y,]

GWAS_flip<-ref_bim_GWAS[ref_bim_GWAS$IUPAC.x == 'R' & ref_bim_GWAS$IUPAC.y == 'Y' | ref_bim_GWAS$IUPAC.x == 'Y' & ref_bim_GWAS$IUPAC.y == 'R' | ref_bim_GWAS$IUPAC.x == 'K' & ref_bim_GWAS$IUPAC.y == 'M' | ref_bim_GWAS$IUPAC.x == 'M' & ref_bim_GWAS$IUPAC.y == 'K',]

GWAS_flip$A1.y.flip<-GWAS_flip$A1.y
GWAS_flip$A1.y.flip[GWAS_flip$A1.y == 'A']<-'T'
GWAS_flip$A1.y.flip[GWAS_flip$A1.y == 'T']<-'A'
GWAS_flip$A1.y.flip[GWAS_flip$A1.y == 'C']<-'G'
GWAS_flip$A1.y.flip[GWAS_flip$A1.y == 'G']<-'C'
GWAS_flip$A1.y<-GWAS_flip$A1.y.flip
GWAS_flip$A1.y.flip<-NULL

GWAS_flip$A2.y.flip<-GWAS_flip$A2.y
GWAS_flip$A2.y.flip[GWAS_flip$A2.y == 'A']<-'T'
GWAS_flip$A2.y.flip[GWAS_flip$A2.y == 'T']<-'A'
GWAS_flip$A2.y.flip[GWAS_flip$A2.y == 'C']<-'G'
GWAS_flip$A2.y.flip[GWAS_flip$A2.y == 'G']<-'C'
GWAS_flip$A2.y<-GWAS_flip$A2.y.flip
GWAS_flip$A2.y.flip<-NULL

GWAS_flip$IUPAC.y[GWAS_flip$A1.y == 'A' & GWAS_flip$A2.y =='G' | GWAS_flip$A1.y == 'G' & GWAS_flip$A2.y =='A']<-'R'
GWAS_flip$IUPAC.y[GWAS_flip$A1.y == 'C' & GWAS_flip$A2.y =='T' | GWAS_flip$A1.y == 'T' & GWAS_flip$A2.y =='C']<-'Y'
GWAS_flip$IUPAC.y[GWAS_flip$A1.y == 'G' & GWAS_flip$A2.y =='T' | GWAS_flip$A1.y == 'T' & GWAS_flip$A2.y =='G']<-'K'
GWAS_flip$IUPAC.y[GWAS_flip$A1.y == 'A' & GWAS_flip$A2.y =='C' | GWAS_flip$A1.y == 'C' & GWAS_flip$A2.y =='A']<-'M'

GWAS_flip_clean<-GWAS_flip[GWAS_flip$IUPAC.x == GWAS_flip$IUPAC.y,]

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Strand was flipped for',dim(GWAS_flip_clean)[1],'variants to match ref_bim.\n')
sink()

GWAS<-rbind(GWAS_noflip,GWAS_flip_clean)

rm(GWAS_flip_clean, GWAS_flip, GWAS_noflip)

names(GWAS)[names(GWAS) == 'A1.y']<-'A1'
names(GWAS)[names(GWAS) == 'A2.y']<-'A2'
GWAS$A1.x<-NULL
GWAS$A2.x<-NULL
GWAS$IUPAC.x<-NULL
GWAS$IUPAC.y<-NULL

#####
# Remove SNPs with missing values
#####

GWAS$CHR<-as.numeric(GWAS$CHR)
GWAS$A1<-as.character(GWAS$A1)
GWAS$A2<-as.character(GWAS$A2)
GWAS$BP<-as.numeric(GWAS$BP)
GWAS$N<-as.numeric(GWAS$N)
GWAS$P<-as.numeric(GWAS$P)

if(sum(names(GWAS) == 'FREQ') == 1){
  GWAS$FREQ<-as.numeric(GWAS$FREQ)
}

if(sum(names(GWAS) == 'OR') == 1){
  GWAS$OR<-as.numeric(GWAS$OR)
}

if(sum(names(GWAS) == 'BETA') == 1){
  GWAS$BETA<-as.numeric(GWAS$BETA)
}

if(sum(names(GWAS) == 'SE') == 1){
  GWAS$SE<-as.numeric(GWAS$SE)
}

GWAS<-GWAS[complete.cases(GWAS),]

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('After removal of SNPs with missing data,',dim(GWAS)[1],' variants remain.\n')
sink()

#####
# Remove SNPs with INFO < opt$info
#####

if(sum(names(GWAS) == 'INFO') == 1){
  GWAS<-GWAS[GWAS$INFO >= opt$info,]
  
  sink(file = paste(opt$output,'.log',sep=''), append = T)
  cat('After removal of SNPs with INFO < ',opt$info,', ',dim(GWAS)[1],' variants remain.\n', sep='')
  sink()
} else {
  sink(file = paste(opt$output,'.log',sep=''), append = T)
  cat('INFO column is not present\n', sep='')
  sink()
}

#####
# Remove SNPs with reported MAF < opt$maf
#####

if(sum(names(GWAS) == 'FREQ') == 1){
  GWAS<-GWAS[GWAS$FREQ >= opt$maf & GWAS$FREQ <= (1-opt$maf),]
  
  sink(file = paste(opt$output,'.log',sep=''), append = T)
  cat('After removal of SNPs with reported MAF < ',opt$maf,', ',dim(GWAS)[1],' variants remain.\n', sep='')
  sink()
} else {
  sink(file = paste(opt$output,'.log',sep=''), append = T)
  cat('Reported MAF column is not present\n', sep='')
  sink()
}

#####
# Remove SNPs with reference MAF < opt$maf
#####

if(!is.na(opt$ref_freq_chr)){
  # Read in reference frequency data
  freq_all<-NULL
  for(i in 1:22){
    freq_i<-fread(paste0(opt$ref_freq_chr,i,'.frq'))
    freq_all<-rbind(freq_all, freq_i)
  }
  
  sink(file = paste(opt$output,'.log',sep=''), append = T)
  cat('ref_freq_chr contains',dim(freq_all)[1],'variants.\n')
  sink()
  
  freq_all<-freq_all[freq_all$MAF > opt$maf & freq_all$MAF < (1-opt$maf),]
  GWAS<-GWAS[(GWAS$SNP %in% freq_all$SNP),]
  
  sink(file = paste(opt$output,'.log',sep=''), append = T)
  cat('After removal of SNPs with reference MAF < ',opt$maf,', ',dim(GWAS)[1],' variants remain.\n', sep='')
  sink()
}

#####
# Remove SNPs with discordant MAF
#####

if(sum(names(GWAS) == 'FREQ') == 1 & !is.na(opt$ref_freq_chr)){
  GWAS_freq<-GWAS[,c('SNP','A1','A2','FREQ')]
  freq_all<-freq_all[,c('SNP','A1','A2','MAF')]
  
  gwas_ref_freq_match<-merge(GWAS_freq,freq_all,by=c('SNP','A1','A2'))
  gwas_ref_freq_switch<-merge(GWAS_freq,freq_all,by.x=c('SNP','A1','A2'),by.y=c('SNP','A2','A1'))
  gwas_ref_freq_switch$FREQ<-1-gwas_ref_freq_switch$FREQ
  gwas_ref_freq<-rbind(gwas_ref_freq_match[,c('SNP','FREQ','MAF')], gwas_ref_freq_switch[,c('SNP','FREQ','MAF')]) 
  gwas_ref_freq$diff<-abs(gwas_ref_freq$MAF-gwas_ref_freq$FREQ)
  
  bitmap(paste0(opt$output,'.MAF_plot.png'), unit='px', res=300, width=1200, height=1200)
    plot(gwas_ref_freq$MAF[gwas_ref_freq$diff > opt$maf_diff],gwas_ref_freq$FREQ[gwas_ref_freq$diff > opt$maf_diff], xlim=c(0,1), ylim=c(0,1), xlab='Reference Allele Frequency', ylab='Sumstat Allele Frequency')
    abline(coef = c(0,1))
  dev.off()
  
  GWAS<-GWAS[(GWAS$SNP %in% gwas_ref_freq$SNP[gwas_ref_freq$diff < opt$maf_diff])]
  
  sink(file = paste(opt$output,'.log',sep=''), append = T)
  cat('After removal of SNPs with absolute MAF difference of < ',opt$maf_diff,', ',dim(GWAS)[1],' variants remain.\n', sep='')
  sink()
} else {
  sink(file = paste(opt$output,'.log',sep=''), append = T)
  cat('Reported MAF column is not present, so discordance with reference cannot be determined.\n', sep='')
  sink()
}

#####
# Insert reference MAF [and check for discordant SD]
#####

if(!is.na(opt$ref_freq_chr)){
  freq_all<-freq_all[,c('SNP','A1','A2','MAF')]
  gwas_ref_freq_match<-merge(GWAS,freq_all,by=c('SNP','A1','A2'))
  gwas_ref_freq_switch<-merge(GWAS,freq_all,by.x=c('SNP','A1','A2'),by.y=c('SNP','A2','A1'))
  gwas_ref_freq_switch$MAF<-1-gwas_ref_freq_switch$MAF
  gwas_ref_freq<-rbind(gwas_ref_freq_match, gwas_ref_freq_switch)
  names(gwas_ref_freq)[names(gwas_ref_freq) == 'MAF']<-'REF.FREQ'
  GWAS<-gwas_ref_freq
  rm(gwas_ref_freq,gwas_ref_freq_switch,gwas_ref_freq_match)

  if(opt$insert_ref_maf == T & sum(names(GWAS) == 'FREQ') == 0){
    keeprefmaf <- TRUE
    sink(file = paste(opt$output,'.log',sep=''), append = T)
    cat('Reference allele frequency inserted.\n', sep='')
    sink()
  } else {
    keeprefmaf <- FALSE
  }
  
  ######
  ### Check for discordant SD in the sumstats
  ### Requires reference allele frequency
  ### c.f. https://doi.org/10.1093/bioinformatics/btaa1029
  ### https://privefl.github.io/bigsnpr/articles/LDpred2.html
  ######
  if(tolower(opt$SDcheck) %in% c("binary","quant")){

    #Compute (validation) SD from reference FREQ:
    #https://privefl.github.io/bigsnpr/articles/LDpred2.html
    GWAS$sd_af <- with(GWAS,
                       sqrt(2 * REF.FREQ * (1 - REF.FREQ)))
    
    #For binary traits eqn numerator is 2, and for quant traits, standard deviation for phenotype (sdY)
     if(tolower(opt$SDcheck)=="binary"){
       numerator <- 2
     } else {
       #Estimate sdY via one of the two methods proposed:
       #https://doi.org/10.1093/bioinformatics/btaa1029
       numerator <-with(GWAS,
                        min(sqrt(0.5)*SE*sqrt(N)))
       # numerator <-with(GWAS,
       #                  median(sd_af*SE*sqrt(N)))
     }
    
    #Denominators for eqns are equivalent for binary and quant phenotypes
    #The code used is equation 1/2 from this paper:
    #https://www.sciencedirect.com/science/article/pii/S2666247722000525
    GWAS$sd_ss <- with(GWAS,
                       numerator / sqrt(N * SE^2 + BETA^2))
    
    #however, an alternate version is proposed in eqn 3 here:
    #https://doi.org/10.1093/bioinformatics/btaa1029
    # GWAS$sd_ss <- with(GWAS,
    #                    numerator /(SE*sqrt(N)))

    #Check for snps with greatly deviating SD:
    #https://doi.org/10.1093/bioinformatics/btaa1029
    GWAS$is_bad <- with(GWAS,
                        sd_ss < (0.5 * sd_af) | sd_ss > (sd_af + 0.1) |
                        sd_ss < 0.1 | sd_af < 0.05)
    
    # #### If GWAS allele frequency provided, repeat SD check with in-sample variant frequency and an imputed N
    # Compute and check Allele frequency validation from the FREQ column provided [If same population ≈ REF.FREQ]
    if("FREQ" %in% colnames(GWAS)){
      sink(file = paste(opt$output,'.log',sep=''), append = T)
      cat('Imputing N for sample based on GWAS-reported allele frequencies, BETA, and SE.\n', sep='')
      sink()
      
      #Standard deviation for GWAS based on allele frequencies
      GWAS$sd_sumstataf <- with(GWAS,
                                sqrt(2 * FREQ * (1 - FREQ)))
      
      #Impute N and repeat SD check relative to the new N value
      GWAS$impN<-with(GWAS,
                      (4/(sd_sumstataf^2)-(BETA^2))/(SE^2))

      #Calculate the in-sample SD based on imputed N
      GWAS$sd_ss_impN <- with(GWAS,
                              numerator / sqrt(impN * SE^2 + BETA^2))
      
      #Run the SD check
      GWAS$is_bad_impN <- with(GWAS,
                               sd_ss_impN < (0.5 * sd_af) | sd_ss_impN > (sd_af + 0.1) |
                               sd_ss_impN < 0.1 | sd_af < 0.05)
      
      sink(file = paste(opt$output,'.log',sep=''), append = T)
      cat('SD check with imputed N suggests retention of:', sum(!GWAS$is_bad_impN), "variants. Check this and plots for compararability to the main SD check.\n", sep='')
      sink()
      
      #Prepare to plot the main SD check and the imputed-N version
      pltwidth<- 2400
    } else {
      #If FREQ is not in the in-sample dataset, skip imputation of N and plot main SD check with a single panel
      #Declare width for subsequent plot
      pltwidth<- 1200
      
      sink(file = paste(opt$output,'.log',sep=''), append = T)
      cat('SD check could not be repeated using an imputed N because FREQ is not present for the GWAS sample.\n', sep='')
      sink()
    }
    
    #Plot the SD check (and if freq column is present, plot a 2nd panel with N imputed from freq)
    #To reduce overplotting, only draw 1/100 snps
    subsamp<- sample(x=c(TRUE,FALSE),size=nrow(GWAS),replace=TRUE,prob=c(0.01,0.99))
    
    bitmap(paste0(opt$output,'.SD_check_plot.png'), unit='px', res=300, width=pltwidth, height=1200)
    
    if("FREQ" %in% colnames(GWAS)){par(mfrow=c(1,2))} #Accommodate 2 panels
    
    #Plot from the main dataset
    plot(GWAS$sd_af[subsamp],GWAS$sd_ss[subsamp],col=ifelse(GWAS$is_bad[subsamp],"black","chartreuse2"),
         xlab="SD from Reference AF",
         ylab = "SD from sumstats",
         main = "Reported N"
         #xlim=c(0,1),ylim=c(0,1)
    )
    abline(a=0,b=1,lty=2)
    legend("topleft", inset=c(0,0), legend=c("TRUE","FALSE"),pch=1, col=c("black","chartreuse2"), title="Failed check?")
    
    if("FREQ" %in% colnames(GWAS)){ #plot 2nd panel
      #Plot from the dataset using an imputed N
      plot(GWAS$sd_af[subsamp],GWAS$sd_ss_impN[subsamp],col=ifelse(GWAS$is_bad_impN[subsamp],"black","chartreuse2"),
           xlab="SD from Reference AF",
           ylab = "SD from sumstats",
           main = "Imputed N"
           #xlim=c(0,1),ylim=c(0,1)
      )
      abline(a=0,b=1,lty=2)
      
      #Drop additional cols
      GWAS$sd_sumstataf <- GWAS$impN <- GWAS$sd_ss_impN <- GWAS$is_bad_impN <- NULL
    }
    dev.off()
    
    #Remove snps failing QC step
    GWAS <- GWAS[!GWAS$is_bad,]
    
    sink(file = paste(opt$output,'.log',sep=''), append = T)
    cat('After removal of SNPs failing SD check based on reported N and reference AF,',dim(GWAS)[1],' variants remain.\n')
    sink()
  
    #Drop columns, and remove ref freq only if not being retained
    GWAS$sd_af <- NULL
    GWAS$sd_ss <- NULL
    GWAS$is_bad <- NULL
  }
  
  #Lastly, keep or remove the reference allele frequency
  if(!keeprefmaf) GWAS$REF.FREQ <- NULL
  rm(keeprefmaf)
}


#####
# Remove SNPs with out-of-bounds p-values
#####

GWAS<-GWAS[GWAS$P <= 1 & GWAS$P > 0,]

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('After removal of SNPs with out-of-bound P values, ',dim(GWAS)[1],' variants remain.\n', sep='')
sink()

#####
# Remove SNPs with duplicated rs numbers
#####

dups<-GWAS$SNP[duplicated(GWAS$SNP)]
GWAS<-GWAS[!(GWAS$SNP %in% dups),]

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('After removal of SNPs with duplicate IDs, ',dim(GWAS)[1],' variants remain.\n', sep='')
sink()

#####
# Remove SNPs with N < 3SD from median N
#####

if(length(unique(GWAS$N)) > 1){
  N_sd<-sd(GWAS$N)
  GWAS<-GWAS[GWAS$N < median(GWAS$N)+(3*N_sd) & GWAS$N > median(GWAS$N)-(3*N_sd),]
  
  sink(file = paste(opt$output,'.log',sep=''), append = T)
  cat('After removal of SNPs with N > ',median(GWAS$N)+(3*N_sd),' or < ',median(GWAS$N)-(3*N_sd),', ',dim(GWAS)[1],' variants remain.\n', sep='')
  sink()
} else {
  sink(file = paste(opt$output,'.log',sep=''), append = T)
  cat('N column is not present or invariant.\n', sep='')
  sink()
}

#####
# Check for genomic control
#####

if(sum(names(GWAS) == 'SE') == 1){
  if(sum(names(GWAS) == 'OR') == 1){
    GWAS$BETA<-log(GWAS$OR)
    GWAS$Z<-GWAS$BETA/GWAS$SE
    GWAS$P_check<-2*pnorm(-abs(GWAS$Z))
    GWAS$Z<-NULL
    GWAS$BETA<-NULL
    
    if(abs(mean(GWAS$P[!is.na(GWAS$P_check)]) - mean(GWAS$P_check[!is.na(GWAS$P_check)])) > 0.01){
       GWAS$P<-GWAS$P_check
       GWAS$P_check<-NULL
      
       sink(file = paste(opt$output,'.log',sep=''), append = T)
       cat('Genomic control detected. P-value recomputed using OR and SE.\n', sep='')
       sink()
    } else {
       sink(file = paste(opt$output,'.log',sep=''), append = T)
       cat('Genomic control was not detected.\n', sep='')
       sink()
       GWAS$P_check<-NULL
    }
  }
  
  if(sum(names(GWAS) == 'BETA') == 1){
    GWAS$Z<-GWAS$BETA/GWAS$SE
    GWAS$P_check<-2*pnorm(-abs(GWAS$Z))
    GWAS$Z<-NULL

    if(abs(mean(GWAS$P[!is.na(GWAS$P_check)]) - mean(GWAS$P_check[!is.na(GWAS$P_check)])) > 0.01){
      GWAS$P<-GWAS$P_check
      GWAS$P_check<-NULL
      
      sink(file = paste(opt$output,'.log',sep=''), append = T)
      cat('Genomic control detected. P-value recomputed using BETA and SE.\n', sep='')
      sink()
    } else {
      sink(file = paste(opt$output,'.log',sep=''), append = T)
      cat('Genomic control was not detected.\n', sep='')
      sink()
      GWAS$P_check<-NULL
    }
    
  }
} else {
  sink(file = paste(opt$output,'.log',sep=''), append = T)
  cat('SE column is not present, genomic control cannot be detected.\n', sep='')
  sink()
}

#####
# Insert SE column
#####

if(sum(names(GWAS) == 'SE') == 0){
  if(sum(names(GWAS) == 'BETA') == 1){
    GWAS$Z<-abs(qnorm(GWAS$P/2))
    GWAS$SE<-abs(GWAS$BETA/GWAS$Z)
  } else {
    GWAS$Z<-abs(qnorm(GWAS$P/2))
    GWAS$SE<-abs(log(GWAS$OR)/GWAS$Z)
  }
  
  GWAS$Z<-NULL
  
  sink(file = paste(opt$output,'.log',sep=''), append = T)
  cat('SE column inserted based on BETA/OR and P.\n', sep='')
  sink()
}

#####
# Remove SNPs with SE == 0
#####

GWAS<-GWAS[GWAS$SE != 0,]

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('After removal of SNPs with SE == 0, ',dim(GWAS)[1],' variants remain.\n', sep='')
sink()

#####
# Write out results
#####

if(file.exists(paste0(opt$output,'.gz'))){
  system(paste0(paste0('rm ',opt$output,'.gz')))
}
if(file.exists(paste0(opt$output))){
  system(paste0(paste0('rm ',opt$output)))
}

fwrite(GWAS, opt$output, sep='\t')

if(opt$gz == T){
  # Compress the GWAS sumstats
  system(paste0('gzip ', opt$output))
}

end.time <- Sys.time()
time.taken <- end.time - start.time
sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Analysis finished at',as.character(end.time),'\n')
cat('Analysis duration was',as.character(round(time.taken,2)),attr(time.taken, 'units'),'\n')
sink()
