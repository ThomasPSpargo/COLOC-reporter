# This script will analyse the pariwise bivariate local genetic correlations between any two phenotypes of interest
# call from the command line as: Rscript lava_rg.R $REF_PREFIX $LOCUS_FILE $INPUT_INFO_FILE $SAMPLE_OVERLAP_FILE "$PHENO1;$PHENO2" $OUTPUT_FILENAME
# (e.g. Rscript lava_rg.R g1000_eur.maf01 blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile input.info.txt sample.overlap.txt 'mdd;scz' 'mdd.scz')
arg = commandArgs(T); refdat = arg[1]; locfile = arg[2]; infofile = arg[3]; overlapfile = arg[4]; phenos = unlist(strsplit(arg[5],";")); outname = arg[6]; ncompar = as.numeric(arg[7])
library(LAVA); library(parallel)

# univ threshold for bivar analysis
univ.thresh = .05 / ncompar

### READ IN DATA ###
loci = read.loci(locfile);  n.loc = nrow(loci)
input = process.input(infofile, overlapfile, refdat, phenos)

#### ANALYSE ####
print(paste("Starting LAVA analysis for",n.loc,"loci"))
progress = ceiling(quantile(1:n.loc, seq(.05,1,.05)))   # (used for printing the progress)
u = b = NULL

# Using mclapply to analyse the loci in parallel
# adjust the mc.cores argument according to the number of available cores (e.g. N avail cores - 1)
# if you dont want to parallelise you can use lapply or a for loop instead
out = mclapply(1:n.loc, mc.cores=15, function(i) {
	if (i %in% progress) print(paste("..",names(progress[which(progress==i)])))  # print progress

	# process locus #https://rdrr.io/github/josefin-werme/lava/man/process.locus.html
	locus = process.locus(loci[i,], input)

	##### Check for genome-wide significance of SNPs in region analysed
	
  ## Identify minimum p-value (imputed from LAVA-calculated Z-score) for SNPs in locus per-trait
  min_snppvals<- lapply(input$sum.stats,function(x,snplist){
    y <- x[x$SNP %in% snplist,]  #Get SNPs in region
    pval<- 2*pnorm(-abs(y$STAT)) #Compute pval from z-score 
    return(min(pval))},snplist=locus$snps)
  
  ## Check for any genome-wide significance
  pvals_aresig<- as.list(min_snppvals<5e-08)
  
  ## Assign informative names for the check
  names(pvals_aresig) <- paste0("min_pval_is_sig_5e.08_phen",1:2)
  names(min_snppvals) <- paste0("min_snp_pval_phen",1:2)
  
  
	# in some cases the locus object cannot be created due to e.g too few SNPs or negative variances in all analysed phenotypes, hence this check
	if (!is.null(locus)) {
		# extract general locus info for output
		loc.info = data.frame(locus = locus$id, chr = locus$chr, start = locus$start, stop = locus$stop, n.snps = locus$n.snps, n.pcs = locus$K)

		# run univ & bivar analysis functions & store results #https://rdrr.io/github/josefin-werme/lava/man/run.univ.bivar.html
		ub = run.univ.bivar(locus, univ.thresh=univ.thresh)
		u = cbind(loc.info, ub$univ)
		if (!is.null(ub$bivar)) {
			b = cbind(loc.info, ub$bivar,min_snppvals,pvals_aresig)
		}
	}
	return(list(univ=u, bivar=b))
})

### WRITE OUTPUT ###
write.table(do.call(rbind, lapply(out,"[[","univ")), paste0(outname,".univ"),sep="\t", row.names=F, quote=F)
write.table(do.call(rbind, lapply(out, "[[", "bivar")), paste0(outname,".bivar"),sep="\t", row.names=F, quote=F)


