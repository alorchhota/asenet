#!/bin/sh
homedir='/home/asaha6/github/asenet'
outdir='results/mono-allelic-filter'
cd $homedir

# create output directory structure
mkdir $outdir

# creating validity data with mono-allelic expression
asedest=$outdir/ase_mono.txt
validdest_old=$outdir/ase_validity_old.txt
minreads=20
pseudo=50
python calculate_ase.py -asedest $asedest -validdest $validdest_old -minreads $minreads -pseudo $pseudo

# creating validity data without mono-allelic expression
asedest=$outdir/ase_wo_mono.txt
validdest_new=$outdir/ase_validity_new.txt
minreads=20
pseudo=50
python calculate_ase.py -asedest $asedest -validdest $validdest_new -minreads $minreads -pseudo $pseudo -monofilter

# analyzing filters
fildest=$outdir/monoallelic_filtered_loci.txt
python analyze_filter.py -newfdata $validdest_new -oldfdata $validdest_old -fildest $fildest