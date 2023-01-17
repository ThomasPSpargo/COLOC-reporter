## COLOC-reporter
___This repository is maintained by Thomas Spargo (thomas.spargo@kcl.ac.uk) - please reach out with any questions___
___

### Performing colocalisation analysis

When describing file structures in this README, "`.`" indicates the main `COLOC-reporter/` directory.

Analysis can be performed using the `./scripts/loopColocaliseRegion.sh` script, which by default runs the most comprehensive analysis protocol available within `./scripts/colocaliseRegion.R`.

Within this bash script, analysis is performed for each row of the `coloc.phenopairs.txt` file, which indicates IDs for pairings of GWAS summary statistics and the genomic regions to examine. The `GWAS_samples.txt` file should indicate metadata linked to the IDs defined in `coloc.phenopairs.txt`, including the configuration of GWAS summary statistic input files.

Plink binary files are also required to provide a reference panel for SNP alignment and determining linkage disequilibrium between SNPs. By default, these files are expected to have the following structure: `./ld_reference/EUR_phase3_chr[1..22].(bed/bim/fam)`. 

The main `colocaliseRegion.R` script can accept several input and analysis configurations. Documentation for these operations are available by running `Rscript ./scripts/colocaliseRegion.R --help`. Alternatively, refer to the `runColocaliseRegion.job` script to view the structure of jobs submitted when running `loopColocaliseRegion.sh`.

Basic setup instructions and details for preparing the input files are provided below.


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
req_packages <- c("tidyverse","optparse","data.table","R.utils","coloc","susieR", "biomaRt","ggrepel", "egg")

#Check for presence of packages
pkg_missing <- req_packages[!req_packages %in% installed.packages()[,"Package"]]

#Return names of missing packages
pkg_missing
```

Most of the required packages can be installed from cran using `install.packages()`. However, `biomarRt` is installed using [bioconductor](https://bioconductor.org/packages/release/bioc/html/biomaRt.html) `BiocManager::install("biomaRt")`.

___PLINK executables___

The PLINK software (v1.90) must be installed before calling to `colocaliseRegion.R`, where it is used for calculating linkage disequilibrium betweeen SNPs.

PLINK2 is not a dependency of the main workflow but is required for the `./scripts/prep_1kg.sh` script, which is used to prepare the LD reference files for the 1000 genomes sample (see the relevant [section](https://github.com/ThomasPSpargo/COLOC-reporter#loopcolocaliseregionsh-input-files) below).


### `loopColocaliseRegion.sh` input files

__coloc.phenopairs.txt__

Analyses to run should be identified in the `./scripts/coloc.phenopairs.txt` file. Columns 1 and 2 provide IDs for two traits whose configurations are detailed in the `GWAS_samples.txt` files (see below). Column 3 indicates the genomic locus to analyse. The genomic locus should be specified as a comma separated list of the format "chromosome,start_position,end_position" (e.g. 17,43460501,44865832).
 
The file should not contain a header, and the main fields must not be comma separated. An example (of the format) is visualised here:
 
 ||||
 |---|---|---|
 |ALS|PD|17,43460501,44865832|
 |PD|SZ	|17,43460501,44865832|
 |SZ|ALS|17,43460501,44865832|

Note: If analysis regions have been identified using the LAVA local genetic correlation [software](https://github.com/josefin-werme/lava), then we provide a script for readily identifying target regions from LAVA output files and generating the file `coloc.phenopairs.txt` (See the __[LAVA integration](https://github.com/ThomasPSpargo/COLOC-reporter#lava-integration])__ section below).


__GWAS_samples.txt__

Information about the input configurations for traits to analyse should be detailed in the `./scripts/GWAS_samples.txt` file. An example this file is as follows:

ID|type|prop|traitSD|pcolumn|statcol|Ncol|chromosome|positions|error|snpcol|freq|traitLabel|FILEPATH
---|---|---|---|---|---|---|---|---|---|---|---|---|---
ALS|cc|0.5|NA|P|BETA|N|CHR|BP|SE|SNP|FREQ|Amyotrophic lateral sclerosis|/path/to/sumstats/file.extension
AZ|cc|0.5|NA|P|BETA|N|CHR|BP|SE|SNP|REF.FREQ|Alzheimer's Disease|/path/to/sumstats/file.extension
FTD|cc|0.5|NA|P|BETA|N|CHR|BP|SE|SNP|REF.FREQ|Frontotemporal dementia|/path/to/sumstats/file.extension
PD|cc|0.5|NA|P|BETA|N|CHR|BP|SE|SNP|FREQ|Parkinson's Disease|/path/to/sumstats/file.extension
SZ|cc|0.5|NA|P|BETA|N|CHR|BP|SE|SNP|REF.FREQ|Schizophrenia|/path/to/sumstats/file.extension

The column headers used in the template should be included. Each column indicates the following:
- ID: Character string indicating identifying string set in the `coloc.phenopairs.txt` input file (see above)
- type: Character string, either 'quant' or 'cc' indicating respectively whether the trait is quantitative or binary (case control)
- prop: For trait where `type=cc`, Indicate proportion of data from cases; ignored if type is quant.
- traitSD: For trait where `type=quant`, optionally indicate standard deviation of trait in the population; ignored if type is cc, and should be set to `NA` if unknown. This will be imputed by coloc for quantitative traits when unknown (see: R function `coloc::sdY.est`).
- pcolumn: Column name for p-values
- statcol: Column name for test statistic, expects column referring to the beta coefficients, if odds ratios given, column must be called 'OR', and these will be converted to the betas. 
- Ncol: Column name detailing sample size per snp. For `type=cc` trait, effective sample size is recommended, in which case `prop` option should be `0.5`
- chromosome: Column name for chromosome
- positions: Column name for genomic position
- error: Column name for test statistic standard error
- snpcol: Column name for SNP ids, rsID is expected
- freq: Column name for allele frequency. Any SNPs with `freq>0.5` will be adjusted (`1-freq`) to obtain minor allele frequency (MAF).
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
sbatch ./scripts/prep_1kg.sh ./ld_reference
```

Note that running `prep_1kg.sh` produces a directory of final size `2.2GB` but additional storage is required for intermediate files. 


### Outputs (documentation incomplete)

__Output files per-analysis__

The file `./coloc/results/<prefix>_coloc/colocalisation.log`, where `<prefix>` identifies the name given to a particular analysis, gives an overview of the results obtained from `colocaliseRegion.R`. The summary includes description of initial summary statistic processing and will also indicate where more detailed outputs for the analyses can be found.

If colocalisation analysis without finemapping is requested, the results of analysis using `coloc.abf` will be returned.

If colocalisation analysis with SuSiE finemapping is requested, SuSiE finemapping results will be returned for each of traits analysed. The results summary of colocalisation analysis with `coloc.susie` will be returned if at least one credible set can be identified per trait. If the finemapping step fails for one or both traits, colocalisation analysis will default to `coloc.abf`, skipping finemapping.

The results of both `coloc.abf` and `coloc.susie` will be returned analysis with both approaches are requested (the default), and if finemapping is successful.

__Concatentating across multiple analyses__

Summary tables of the overall results of colocalisation analyses (performed using the functions `coloc.abf` and/or `coloc.susie` from the R `coloc` [package](https://chr1swallace.github.io/coloc/articles/a01_intro.html)) are returned when running `colocaliseRegion.R`. 

When performing a series of analyses using `loopColocaliseRegion.sh`, this summary will be returned for each analysis individually. Running the following script will concatenate the colocalisation results from multiple analyses into summary files:

```
#The trailing argument points to the parent directory to which all results have been returned
bash ./scripts/regionLoopSummary.sh ./coloc/results
```

The file `./results/summary_all_coloc_abf.csv` is returned if any results from coloc.abf are identified, and equivalently `./results/summary_all_coloc_susie.csv` is returned when identifying results from `coloc.susie`.


### LAVA integration

[LAVA](https://github.com/josefin-werme/lava) is a software for performing local genetic correlation analysis. It is one approach for identifying genomic regions where variants may colocalise.

We provide facilities for extracting and visualising outputs from LAVA.

__Extracting regions with significant local genetic correlation from LAVA outputs__

The `coloc.phenopairs.txt` input file can be readily generated from LAVA outputs by calling the `./scripts/extractLAVA.R` script, directing to a directory which contains the results of bivariate genetic correlations using LAVA. By default, the script expects files with the extension '.bivar'

An example use of this script is provided in `runextractLAVA.sh`.

Further details for options when using this R script are available by calling `Rscript ./scripts/extractLAVA.R --help`.

__Visualising univariate and bivariate analysis from LAVA__

The `./scripts/plotLAVAoutputs.R` script generates plots to visualise analyses from LAVA.

When applied to a bivariate correlation output file, two plots are returned: First, a plot which visualises genetic correlations at examined loci; Second, a manhattan-style plot indicating p-values for each locus.

When applied to univariate heritability analysis from LAVA, a single manhattan plot is returned.
This analysis is used as a filtering step in LAVA prior to bivariate analysis. Only loci passing the set p-value threshold of univariate analysis for both traits will have been carried forward to bivariate genetic correlation analysis.

This visualisation script can be applied recursively to all univariate and bivariate outputs contained within a given directory (and indicated respectively using `.univ` and `.bivar` file extensions), by applying the `./scripts/runplotLAVAoutputs.sh` script.

Further details for options when using the LAVA visualisation R script are available by calling `Rscript ./scripts/plotLAVAoutputs.R --help`.
