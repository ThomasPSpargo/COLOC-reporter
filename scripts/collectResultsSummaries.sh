#!/bin/bash
#####
# Author: Thomas Spargo (thomas.spargo@kcl.ac.uk)
# Obtained from GitHub repository: https://github.com/ThomasPSpargo/COLOC-reporter
#
# This script concatenates finemapping and colocalisation results produced while running COLOC-reporter across a series of inputs.
# This script should be called from the command line with a trailing argument indicating the directory to which all different analyses were returned
#####

#Directory within which to search
dir=${1:-.}

######
### Concatenate across all coloc.abf results files
######
echo "" #Echo blank line here and below for readability of command line output
declare -a abfRes=($(echo ${dir}/**/tables/results_summary_coloc_abf.csv))

#Logic check testing whether files have been identified
if [[ $abfRes == *"**"* ]]; then
	echo "No coloc.abf results identified"
else
	echo "Concatenating coloc.abf results:"
	#Extract header from the first file
	cat ${abfRes[0]} | head -n 1 > ${dir}/summary_all_coloc_abf.csv

	#Across all files minus header
	for i in ${abfRes[@]}; do
		echo ${i}
		cat ${i} | sed -n '1!p' >> ${dir}/summary_all_coloc_abf.csv
	done
fi


######
### Concatenate across all coloc.susie results files
######
echo ""
declare -a clcSusieRes=($(echo ${dir}/**/tables/results_summary_coloc_susie.csv))

#Logic check testing whether files have been identified
if [[ $clcSusieRes == *"**"* ]]; then
	echo "No coloc.susie results identified"
else
	echo "Concatenating coloc.susie results:"
	
	#Extract header from the first file
	cat ${clcSusieRes[0]} | head -n 1 > ${dir}/summary_all_coloc_susie.csv

	#concatenate across all files minus header
	for j in ${clcSusieRes[@]}; do
		echo ${j}
		cat ${j} | sed -n '1!p' >> ${dir}/summary_all_coloc_susie.csv
	done
fi

######
### Concatenate across all susie finemapping results summaries
######
echo ""
declare -a finemapRes=($(echo ${dir}/**/tables/results_summary_finemapping.csv))

#Logic check testing whether files have been identified
if [[ $finemapRes == *"**"* ]]; then
	echo "No finemapping results from SuSiE identified"
else
	echo "Concatenating SuSiE finemapping results:"
	#Extract header from the first file
	cat ${finemapRes[0]} | head -n 1 > ${dir}/summary_all_finemapping.csv

	#Across all files minus header
	for i in ${finemapRes[@]}; do
		echo ${i}
		cat ${i} | sed -n '1!p' >> ${dir}/summary_all_finemapping.csv
	done
fi
echo ""