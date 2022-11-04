#!/bin/bash
#SBATCH --job-name=prep_1kg
#SBATCH --mem=20G
#SBATCH --time=0-10:00
#SBATCH --ntasks=5
#SBATCH --output=%x.log


#Syntax adapted from:
#https://cran.r-project.org/web/packages/snpsettest/vignettes/reference_1000Genomes.html

#Data available as plink2 resources, dropbox links may change
#https://www.cog-genomics.org/plink/2.0/resources#1kg_phase3

#Note to use 64bit version of plink2, else memory is constrained to 2047mb


#####
## INPUT VARIABLES
#####

#Make output directory if does not exist and set as wd
#cd /scratch/users/k1802739/coloc/ld_reference
mkdir -p ${1}
cd ${1}

#Define path to plink2 software; defaults to just 'plink2'
plink2path=${2:-plink2}

#Define population subset for which to extract LD reference; defaults to EUR
pop=${3:-EUR}


#####
## BEGIN
#####

# The links may change in the future
# -O to specify output file name, -q for quiet
if [ ! -f "all_phase3.psam" ]; then
	wget -q -O all_phase3.psam "https://www.dropbox.com/s/6ppo144ikdzery5/phase3_corrected.psam?dl=1"
fi
if [ ! -f "all_phase3.pgen.zst" ]; then
	wget -q -O all_phase3.pgen.zst "https://www.dropbox.com/s/y6ytfoybz48dc0u/all_phase3.pgen.zst?dl=1"
fi
if [ ! -f "all_phase3.pvar.zst" ]; then
	wget -q -O all_phase3.pvar.zst "https://www.dropbox.com/s/odlexvo8fummcvt/all_phase3.pvar.zst?dl=1"
fi

if [ ! -f "all_phase3.pgen" ]; then
	# Decompress pgen.zst to pgen
	plink2 --zst-decompress all_phase3.pgen.zst > all_phase3.pgen
fi

if [ ! -f "${pop}_1kg_samples.txt" ]; then
# Prepare sub-population filter file
	awk -v pop=${pop} 'NR == 1 || $5 == pop {print $1}' all_phase3.psam > ${pop}_1kg_samples.txt
fi

# "vzs" modifier to directly operate with pvar.zst
# "--chr 1-22" excludes all variants not on the listed chromosomes
# "--output-chr 26" uses numeric chromosome codes
# "--max-alleles 2": PLINK 1 binary does not allow multi-allelic variants
# "--rm-dup" removes duplicate-ID variants
# "--set-missing-var-id" replaces missing IDs with a pattern

# Run recursively by chromosome to reduce burden of intermediate file storage
for i in {1..22}; do
	$plink2path --pfile all_phase3 vzs \
				--chr $i \
				--output-chr 26 \
				--max-alleles 2 \
				--rm-dup exclude-mismatch \
				--set-missing-var-ids '@_#_$1_$2' \
				--make-pgen \
				--keep ${pop}_1kg_samples.txt \
				--out ${pop}_phase3_chr${i}

	# pgen to bed
	# "--maf 0.005" remove most monomorphic SNPs 
	# (still may have some when all samples are heterozyguous -> maf=0.5)
	$plink2path --pfile ${pop}_phase3_chr${i} \
				--maf 0.001 \
				--make-bed \
				--out ${pop}_phase3_chr${i}
	
	#Clear per-chromosome intermediate files
	rm ${pop}_phase3_chr${i}.p* 

# Split bed/bim/fam by chromosome
#for i in {1..22}
#	do $plink2path --bfile ${pop}_phase3_autosomes --chr $i --make-bed --out ${pop}_phase3_chr$i
done

#Clear initial files
rm all_phase3.*
