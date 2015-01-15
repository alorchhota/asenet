#!/bin/sh
homedir='/home/asaha6/github/asenet'
outdir='results'

cd $homedir

# create output directory structure
mkdir $outdir

# calculate isoform ratio correlations
python map_ase_locus_to_gene.py -home $homedir

