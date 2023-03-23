# this 'settings' file contains relevant directories and file names used for the analyses
# you need to at least adapt the MAIN path 

# key directories
MAIN=/scratch/users/k1802739/LAVA 	# base directory ( ADAPT !!! )
SCRIPTS=$MAIN/scripts/02_LAVA_analyses		# scripts dir
OUT=$MAIN/results		# results dir 
DATA=$MAIN/inputs			# data dir

# names of key input files
LOCFILE=blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile
INFOFILE=input.info.txt
OVERLAP=sample.overlap.txt

# reference data dir & prefix. Since 
REFDIR=$DATA # (usually this is not in the main data dir since I only keep a single copy of these files)
REFDAT=EUR_autosomes