## COLOC-reporter
___This repository is maintained by Thomas Spargo (thomas.spargo@kcl.ac.uk) - please reach out with any questions___

_Last update 04/11/2022_

___

### Performing colocalisation analysis

When describing file structures in this README, "`.`" indicates the main `COLOC-reporter/` directory.

Analysis can be performed using the `./scripts/loopColocaliseRegion.sh` script, which by default runs the most comprehensive analysis protocol available within `./scripts/colocaliseRegion.R`.

Within this bash script, analysis is performed for each row of the `coloc.phenopairs.txt` file, which indicates IDs for pairings of GWAS summary statistics and the genomic regions to examine. The configuration of traits and GWAS summary statistic inputs should be specified in the `GWAS_samples.txt` file, with IDs linked to those specified within `coloc.phenopairs.txt`.

Plink binary files are also required to provide a reference panel for SNP alignment and determining linkage disequilibrium between SNPs. By default, these files are expected to have the following structure: `./ld_reference/EUR_phase3_chr[1..22].(bed/bim/fam)`. 

The main `colocaliseRegion.R` script can accept several input and analysis configurations. Documentation these operations are available by running `Rscript ./scripts/colocaliseRegion.R --help`. Alternatively, refer to the `runColocaliseRegion.job` script to view the structure of jobs submitted when running `loopColocaliseRegion.sh`

Basic setup instructions and details for preparing the input files are provided below.


### Setup

__Clone the repository__

```
git clone https://github.com/ThomasPSpargo/COLOC-reporter.git
cd COLOC-reporter
```

__Install R package dependencies__

Various R packages are required. Please ensure all the required packages are available to load.

Running the following in `R` checks for necessary packages and the `pkg_missing` object will name any packages that need to be installed.

```
#Name required packages
req_packages <- c("optparse", "dplyr","data.table","R.utils","dplyr","coloc","susieR", "biomaRt","ggrepel","ggplot2", "egg")

#Check for presence of packages
pkg_missing <- req_packages[!req_packages %in% installed.packages()[ , "Package"]]

#Return names of missing packages
pkg_missing
```

Most of the required packages can be installed from cran using `install.packages()`. However, `biomarRt` is installed using [bioconductor](https://bioconductor.org/packages/release/bioc/html/biomaRt.html) `BiocManager::install("biomaRt")`.

### `loopColocaliseRegion.sh` input files

__coloc.phenopairs.txt__

Analyses to run should be identified in the `./scripts/coloc.phenopairs.txt` file. Columns 1 and 2 provide IDs for two traits whose configurations are detailed in the `GWAS_samples.txt` files (see below). Column 3 indicates the genomic locus to analyse. The genomic locus should be specified as a comma separated list of the format "chromosome,start_position,end_position" (e.g. 17,43460501,44865832).
 
The file should not contain a header, and the main fields must not be comma separated. An example (of the format) is visualised here:
 
 ||||
 |---|---|---|
 |ALS|PD|17,43460501,44865832|
 |PD|SZ	|17,43460501,44865832|
 |SZ|ALS|17,43460501,44865832|

Note: If analysis regions have been identified using the LAVA local genetic correlation [software](https://github.com/josefin-werme/lava), then we provide a script for readily identifying target regions from LAVA output files and generating the file `coloc.phenopairs.txt` (See the __LAVA integration__ section below).


__GWAS_samples.txt__

Information about the input configurations for traits to analyse should be detailed in the `./scripts/GWAS_samples.txt` file. An example this file is as follows:

ID|type|prop|pcolumn|statcol|Ncol|chromosome|positions|error|snpcol|MAF|traitLabel|FILEPATH
---|---|---|---|---|---|---|---|---|---|---|---|---
ALS|cc|0.5|P|BETA|N|CHR|BP|SE|SNP|FREQ|Amyotrophic lateral sclerosis|/path/to/sumstats/file.extension
AZ|cc|0.5|P|BETA|N|CHR|BP|SE|SNP|REF.FREQ|Alzheimer's Disease|/path/to/sumstats/file.extension
FTD|cc|0.5|P|BETA|N|CHR|BP|SE|SNP|REF.FREQ|Frontotemporal dementia|/path/to/sumstats/file.extension
PD|cc|0.5|P|BETA|N|CHR|BP|SE|SNP|FREQ|Parkinson's Disease|/path/to/sumstats/file.extension
SZ|cc|0.5|P|BETA|N|CHR|BP|SE|SNP|REF.FREQ|Schizophrenia|/path/to/sumstats/file.extension

The column headers used in the template should be included. Each column indicates the following:
- ID: Character string indicating identifying string set in the `coloc.phenopairs.txt` input file (see above)
- type: Character string, either 'quant' or 'cc' indicating respectively whether the trait is quantitative or binary (case control)
- prop: For trait where `type=cc`, Indicate proportion of data from cases; ignored if type is quant
- pcolumn: Column name for p-values
- statcol: Column name for test statistic, expects column referring to the beta coefficients, if odds ratios given, column must be called 'OR', and these will be converted to the betas. 
- Ncol: Column name detailing sample size per snp. For type=cc trait, effective sample size is recommended, in which case `prop` option should be `0.5`
- chromosome: Column name for chromosome
- positions: Column name for genomic position
- error: Column name for test statistic standard error
- snpcol: Column name for SNP ids, rsID is expected
- MAF: Column name for allele frequency
- traitLabel: A long name for the trait, to be used for plotting
- FILEPATH: Path to summary statistics file (.gz compressed files are allowed)

__LD reference panel__

A genomic reference dataset should be provided, and the path to per-chromosome plink binary files should be specified in the `loopColocaliseRegion.sh` script.

By default `loopColocaliseRegion.sh` expects a file structure of: `./ld_reference/EUR_phase3_chr[1..22].[bed/bim/fam]`. However, any per-chromosome LD reference panel can be used.

The default file structure can be generated by running the `./scripts/prep_1kg.sh` script to generate LD reference panels for chromosomes 1-22 based on the 1000 genomes phase 3 data release. The script accepts 3 trailing arguments, which define:
1. Path in which to save the extracted files

1. Path to plink2 executable, defaults to 'plink2'

1. Desired population to extract, defaults to 'EUR'

To generate the default file structure (within a slurm scheduler), run 
```
sbatch ./scripts/prep_1kg.sh ./ld_reference
```

Note that running `prep_1kg.sh` produces a directory of final size `2.2GB` but additional storage is required for intermediate files.

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
