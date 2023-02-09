## COLOC-reporter
___This repository is maintained by Thomas Spargo (thomas.spargo@kcl.ac.uk) - please reach out with any questions___
___

### Performing colocalisation analysis

When describing file structures in this README, "`.`" indicates the main `COLOC-reporter/` directory.

Analysis can be performed using the `./scripts/loopColocaliseRegion.sh` script, which by default runs the most comprehensive analysis protocol available within `./scripts/colocaliseRegion.R`, and applies `./scripts/gwasSummaryPlotter.R` to compare harmonised summary statistic results between pairs of traits prepared for colocalisation analysis.

Within this bash script, analysis is performed for each row of the `set.regions.txt` file, which indicates IDs for GWAS summary statistics and the genomic regions to examine. The `GWAS_samples.txt` file should indicate metadata linked to the IDs defined in `set.regions.txt`, including the configuration of GWAS summary statistic input files.

Plink binary files are also required to provide a reference panel for SNP alignment and determining linkage disequilibrium between SNPs. By default, these scripts are configured to expect the following structure: `./ld_reference/EUR_phase3_chr[1..22].(bed/bim/fam)`. 

The main `colocaliseRegion.R` script can accept several input and analysis configurations. Documentation for these operations are available by running `Rscript ./scripts/colocaliseRegion.R --help`.

Likewise, `gwasSummaryPlotter.R` can be used to compare measures provided for summary statistic pairs (e.g. for beta values, p-values, finemapping PIPs from the SuSiE software). Documentation for this script is available by calling: `Rscript ./scripts/gwasSummaryPlotter.R --help`

The structure of the `runColocaliseRegion.job` script provides an example of how the analysis can be configured when running `loopColocaliseRegion.sh`.

Basic setup instructions and details for preparing the input files are provided below. Some further information about valid analysis configurations is also provided.


### Setup

__Clone the repository__

```
git clone https://github.com/ThomasPSpargo/COLOC-reporter.git
cd COLOC-reporter
```

__Install dependencies__

The most straightforward way to ensure all the necessary dependencies are installed is to run analyses within the provided [conda](https://docs.conda.io/en/latest/miniconda.html) environment.

The steps for this are as follows:

```
#Create the environment from the .yml file provided
conda env create --file coloc_reporter.yml

#Activate the conda environment
conda activate coloc_reporter

```

Alternatively, COLOC reporter dependencies may be configured manually:

___R packages___

Various R packages are required. Running the following in `R` checks for necessary packages and the `pkg_missing` object will name any requiring installation.

```
#Name required packages
req_packages <- c("tidyverse","optparse","data.table","R.utils","coloc","susieR", "biomaRt","ggrepel", "patchwork")

#Check for presence of packages
pkg_missing <- req_packages[!req_packages %in% installed.packages()[,"Package"]]

#Return names of missing packages
pkg_missing
```

Most of the required packages can be installed from cran using `install.packages()`. However, `biomarRt` is installed using [bioconductor](https://bioconductor.org/packages/release/bioc/html/biomaRt.html) `BiocManager::install("biomaRt")`.

___PLINK executables___

The PLINK software (v1.90) must be installed before calling to `colocaliseRegion.R`, where it is used for calculating linkage disequilibrium betweeen SNPs.

PLINK2 is not a dependency of the main workflow but is required for the `./scripts/prep_1kg.sh` script, which is used to prepare the LD reference files with GRCh37 alignment for the 1000 genomes sample (see the relevant [section](https://github.com/ThomasPSpargo/COLOC-reporter#loopcolocaliseregionsh-input-files) below).


### `loopColocaliseRegion.sh` input files

__set.regions.txt__

Analyses to run should be identified in the `./scripts/set.regions.txt` file. The first column provides a comma-delimited list of IDs for traits with configurations declared in the `GWAS_samples.txt` file (see below). Column two indicates the genomic region to analyse. The genomic locus should be specified as a comma separated list of the format "chromosome,start_position,end_position" (e.g. 17,43460501,44865832).  Any subsequent columns within this file are ignored.


 
The file should contain a header row (which is ignored) and fields should be tab-delimited. An example of the minimal format for the file, including several possible trait configurations is visualised here:
 
 |traits|region|
 |---|---|
 |ALS,PD|17,43460501,44865832|
 |PD,SZ	|17,43460501,44865832|
 |SZ,ALS|17,43460501,44865832|
 |SZ,ALS,PD|17,43460501,44865832|
 |PD|17,43460501,44865832|

The first 3 rows of this example file indicate that each analysis (a row) should be performed using pairs of traits. Indicating two traits is compatible with the full analysis protocol from initial harmonisation of summary statistics with the reference panel to colocalisation analysis.

Rows 4 and 5 respectively indicate that a region should be analysed across 3 traits or a single trait. All steps prior to colocalisation analysis are performed when supplying either 1 or 3+ traits. When 1 trait is supplied, analysis stops after fine-mapping is complete. When 3+ traits are supplied, only the first two traits are currently taken forward for colocalisation analysis.

The `--runMode` option of `colocaliseRegion.R` allows further control over the analysis protocol (e.g. setting `finemapOnly` declares that only processing steps until and including finemapping should be performed)

Note: If analysis regions have been identified using the LAVA local genetic correlation [software](https://github.com/josefin-werme/lava), then we provide a script for readily identifying target regions from LAVA output files and generating the file `set.regions.txt` (See the __[LAVA integration](https://github.com/ThomasPSpargo/COLOC-reporter#lava-integration])__ section below).


__GWAS_samples.txt__

Information about the input configurations for traits to analyse should be detailed in the `./scripts/GWAS_samples.txt` file. An example this file is as follows:

ID|type|prop|traitSD|p_col|stat_col|N_col|chr_col|pos_col|se_col|snp_col|freq_col|traitLabel|FILEPATH
---|---|---|---|---|---|---|---|---|---|---|---|---|---
ALS|cc|0.5|NA|P|BETA|N|CHR|BP|SE|SNP|FREQ|Amyotrophic lateral sclerosis|/path/to/sumstats/file.extension
AZ|cc|0.5|NA|P|BETA|N|CHR|BP|SE|SNP|REF.FREQ|Alzheimer's Disease|/path/to/sumstats/file.extension
FTD|cc|0.5|NA|P|BETA|N|CHR|BP|SE|SNP|REF.FREQ|Frontotemporal dementia|/path/to/sumstats/file.extension
PD|cc|0.5|NA|P|BETA|N|CHR|BP|SE|SNP|FREQ|Parkinson's Disease|/path/to/sumstats/file.extension
SZ|cc|0.5|NA|P|BETA|N|CHR|BP|SE|SNP|REF.FREQ|Schizophrenia|/path/to/sumstats/file.extension

The column headers used in the template should be included. Each column indicates the following:
- ID: Character string indicating trait ID strings set in the `set.regions.txt` input file (see above)
- type: Character string, either 'quant' or 'cc' indicating respectively whether the trait is quantitative or binary (case control)
- prop: For trait where `type=cc`, Indicate proportion of data from cases; ignored if type is quant
- traitSD: For trait where `type=quant`, optionally indicate standard deviation of trait in the population; ignored if type is cc, and should be set to `NA` if unknown. This will be imputed by coloc for quantitative traits when unknown (see: R function `coloc::sdY.est`).
- p_col: Column name for p-values
- stat_col: Column name for test statistic, expects column referring to the beta coefficients, if odds ratios given, column must be called 'OR', and these will be converted to the betas. 
- N_col: Column name detailing sample size per snp. For `type=cc` trait, effective sample size is recommended, in which case `prop` option should be `0.5`
- chr_col: Column name for chromosome
- pos_col: Column name for genomic position
- se_col: Column name for standard error of beta coefficients
- snp_col: Column name for SNP ids, rsID is expected
- freq_col: Column name for allele frequency. Any SNPs with `freq>0.5` will be adjusted (`1-freq`) to obtain minor allele frequency (MAF).
- traitLabel: A long name for the trait, to be used for plotting
- FILEPATH: Path to summary statistics file (`.gz` compressed files can be used)

__LD reference panel__

A genomic reference dataset should be provided, and the path to per-chromosome plink binary files should be specified in the `loopColocaliseRegion.sh` script.

By default `loopColocaliseRegion.sh` expects a file structure of: `./ld_reference/EUR_phase3_chr[1..22].[bed/bim/fam]`. However, any per-chromosome LD reference panel can be used.

The default file structure can be generated by running the `./scripts/prep_1kg.sh` script to generate LD reference panels for chromosomes 1-22 based on the 1000 genomes phase 3 data release. The script accepts 3 trailing arguments, which define:
1. Directory in which to save the extracted files

1. Path to plink2 executable, defaults to 'plink2'

1. Desired population to extract, defaults to 'EUR'

To generate the default file structure (within a slurm scheduler), run 
```
sbatch -p <partition> ./scripts/prep_1kg.sh ./ld_reference
```

Note that running `prep_1kg.sh` produces a directory of final size `2.2GB` but additional storage is required for intermediate files. 

### Analysis protocol

Analysis can be include only finemapping, only colocalisation, or a combination of these steps.

Colocalisation analysis performed without prior finemapping uses the `coloc.abf` function. Colocalisation analysis including the SuSiE finemapping step is performed with the `coloc.susie` function. Finemapping analysis is performed using the `runsusie` function, which is a wrapper around `susieR::susie_rss` and is provided within the COLOC package. Note that colocalisation analysis with `coloc.susie` is only possible when at least one credible set is identified in both of two traits during the finemapping step.

Finemapping results will be returned for all traits defined for a given analysis. Any number of traits can be defined up to and including this step. However, note that only SNPs common to all supplied traits and the LD reference will be retained.

Colocalisation analysis can be performed only for pairs of traits. Therefore, this is currently performed using only the first two traits defined. If more than two traits are supplied for an analysis which is declared to include colocalisation analysis, then steps prior to colocalisation are completed for all traits before subsetting to the first two traits declared.

The main analysis protocol is controlled using the `--runMode` option of `colocaliseRegion.R`. This has 4 possible settings:
- `doBoth` (the default):  both `coloc.abf` and `coloc.susie` will run, assuming finemapping is successful for both traits supplied to `coloc.susie`
- `trySusie`: `coloc.susie` will be used if SuSiE finemapping identifies at least 1 credible set in each trait and `coloc.abf` is returned otherwise
- `skipSusie`: only `coloc.abf` will be applied, and processes necessary for `coloc.susie` are skipped (e.g. no need to call to plink and generate the linkage disequilibrium matrix)
- `finemapOnly`: finemapping is performed for each trait and any colocalisation analyses are skipped

### Quality control

Quality control checks are performed prior to finemapping analysis to test whether the GWAS summary statistic effect estimates are consistent with the LD matrix supplied. When an LD matrix is provided from an out-of-sample population, inconsistency could indicate a poor ancestry match between the reference and GWAS samples.

The main results from quality control checks performed are returned in the `./coloc/results/<prefix>_coloc/finemapQC` directory.

__Testing consistency between GWAS and LD matrix__

The SuSiE [`estimate_s_rss`](https://stephenslab.github.io/susieR/reference/estimate_s_rss.html) function returns a global measure of consistency between the summary statistic Z-scores (calculated from `beta` and `SE`) and LD matrix provided (see function documentation for details). This measure is a value between 0 and 1 where higher numbers indicate greater inconsistency. This can help identify whether to be cautious in interpreting finemapping results. However, there is not a clear threshold defining which values are reasonable. 

We perform this check in 2 stages.

First, consistency is checked across SNPs available (using all 3 estimation methods provided) after harmonising data across traits and to the reference population. If any SNPs are identified in the Z-score outliers step (see below) as having potentially reversed allele order, the check is repeated after flipping the valence of the effect estimate for these putatively reversed SNPs.

Second, consistency is checked using the default 'null-mle' method only  for the dataset analysed by susie, and this result is returned alongside the overall finemapping results in the column 'LD_Zscore_consistency'. Note that this dataset may correspond exactly to the one from the first check, but could contain either fewer SNPs or some SNPs with flipped effect estimates according to the results of the Z-score outliers check and analysis settings.

We suggest evaluating this measure in the context of your results. In particular, the reasonability of any credible sets identified and the SNP-wise PIP estimates should be checked against GWAS p-values, betas, or Z-scores.

__Identifying Z-score outliers__

The SuSiE [`kriging_rss`](https://stephenslab.github.io/susieR/reference/kriging_rss.html) function checks whether observed Z-scores correspond with those expected based on information from the LD matrix and other SNPs (see function documentation for details). As part of this check, this identifies any SNPs with test statistic effect estimates which might be inverted. These SNPs are marked in red on the Z-score plot returned by the function.

We use this function to:
- Obtain a visual indication of discordance between observed and expected Z-scores
- Identify SNPs which may have flipped effect estimate encoding (Note: allele order was earlier harmonised with the reference and therefore the valence of the effect estimates should be congruent with the LD matrix)
- Visualise positions of SNPs assigned to in credible sets relative to observed and expected Z-scores

If any SNPs are flagged as having a potentially flipped effect estimate, we:
- Always repeat the GWAS and LD matrix consistency check in a dataset with these values flipped (see above)
- Offer options (controlled by `--finemapQC_handleBetaFlips`) for proceding with finemapping and subsequent colocalisation analysis steps. These options are to:
 - continue with an unchanged datset
 - omit (across all traits) any SNPs identified as having potentially erroneous encoding
 - to flip the direction of the beta statistic for those SNPs flagged on a trait-by-trait basis
 
 By default, we omit SNPs with putatively flipped effect estimates from the analysis, since these can substantially affect the SuSiE finemapping result.
 
 The option to invert the effect estimates may improve convergence of SuSiE but should be used cautiously. This is because allele order has been harmonised between the the LD reference and the summary statistics, and therefore this flip is performed under the assumption of an allele order encoding error in the summary statistics.

__Checking convergence__

The number of iterations required for successful convergence of the SuSiE finemapping model is returned in the analysis log file. The [manuscript](https://doi.org/10.1371/journal.pgen.1010299) describing the approach notes that a high number of iterations can result from large inconsistencies between the LD matrix and the summary statistics. Considering the number of iterations required may help further identify if finemapping issues have occured.

### Outputs

__Output files per-analysis__

_Analysis summary_

The file `./coloc/results/<prefix>_coloc/colocalisation.log`, where `<prefix>` identifies the name given to a particular analysis, gives an overview of the results obtained from `colocaliseRegion.R`. The summary includes description of initial summary statistic processing and indicates where more detailed analysis outputs can be found.

The directory `./coloc/results/<prefix>_coloc/tables/` contains all tabular summaries generated within completed analyses. The most important outputs, which overview finemapping and colocalisation result, can be found in:
- `results_summary_coloc_abf.csv`
- `results_summary_coloc_susie.csv`
- `results_summary_finemapping.csv`

The directory `./coloc/results/<prefix>_coloc/plots/` contains various figures generated to facilitate interpretation of finemapping and colocalisation results. These should be interpreted in the context the main finemapping and colocalisation results.

The `colocalisation.log` file gives some details about each file returned.

_Files for further analysis_

The directory `./coloc/results/<prefix>_coloc/data/` stores primarily .Rds files generated across the main analysis steps that can be readily read into R (using the `readRDS` function) for additional analyses.

- `data/datasets/` contains:
	- coloc-formatted datasets for each trait analysed in the region after snps have been: harmonised across traits, to the reference provided, and then formatted for compatibility with the functions `coloc.susie`, `coloc.abf`, and `runsusie`.
	- a .csv file of harmonised summary statistics for the region across traits analysed, which can be supplied to the `gwasSummaryPlotter.R` script.

- `data/finemapping/` contains outputs of the `runsusie` function for each trait with at least 1 credible set identified. This can be supplied directly to `coloc.susie`, or further examined in accordance with your needs.

- `data/colocalisation/` contains outputs of completed analyses between trait pairs using `coloc.susie` and `coloc.abf`. These can be further examined in accordance with your neeeds.

Some files returned within `./coloc/results/<prefix>_coloc/tables/` may also be helpful for performing additional investigations not implemented within this workflow.

__Concatentating across multiple analyses__

Summary tables of the overall results of finemapping and colocalisation analyses (performed using the functions `coloc.abf`, `coloc.susie` and, `runsusie` from the R `coloc` [package](https://chr1swallace.github.io/coloc/articles/a01_intro.html)) are returned when running `colocaliseRegion.R`. 

When performing a series of analyses with `loopColocaliseRegion.sh`, this summary will be returned for each analysis individually. Running the following script will concatenate results summaries across these steps from multiple analyses into summary files:

```
#The trailing argument points to the parent directory to which all results have been returned
bash ./scripts/regionLoopSummary.sh ./coloc/results
```

The file `./results/summary_all_coloc_abf.csv` is returned if any results from `coloc.abf` are identified, and equivalently `./results/summary_all_coloc_susie.csv` is returned when identifying results from `coloc.susie`.

The results from univariate finemapping performed with SuSiE are returned in the file `./results/summary_all_finemapping.csv`.


### LAVA integration

[LAVA](https://github.com/josefin-werme/lava) is a software for performing local genetic correlation analysis. It is one approach for identifying genomic regions where variants may colocalise.

We provide facilities for extracting and visualising outputs from LAVA.

__Extracting regions with significant local genetic correlation from LAVA outputs__

The `set.regions.txt` input file can be readily generated from LAVA outputs by calling the `./scripts/extractLAVA.R` script, directing to a directory which contains the results of bivariate genetic correlations using LAVA. By default, the script expects files with the extension '.bivar'

An example use of this script is provided in `runextractLAVA.sh`.

Further details for options when using this R script are available by calling `Rscript ./scripts/extractLAVA.R --help`.

__Visualising univariate and bivariate analysis from LAVA__

The `./scripts/plotLAVAoutputs.R` script generates plots to visualise analyses from LAVA.

When applied to a bivariate correlation output file, two plots are returned: First, a plot which visualises genetic correlations at examined loci; Second, a manhattan-style plot indicating p-values for each locus.

When applied to univariate heritability analysis from LAVA, a single manhattan plot is returned.
This analysis is used as a filtering step in LAVA prior to bivariate analysis. Only loci passing the set p-value threshold of univariate analysis for both traits will have been carried forward to bivariate genetic correlation analysis.

This visualisation script can be applied recursively to all univariate and bivariate outputs contained within a given directory (and indicated respectively using `.univ` and `.bivar` file extensions), by applying the `./scripts/runplotLAVAoutputs.sh` script.

Further details for options when using the LAVA visualisation R script are available by calling `Rscript ./scripts/plotLAVAoutputs.R --help`.
