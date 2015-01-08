homedir='/home/asaha6/github/asenet'
outdir='results/TF-ASE-Correlation-2015-01-07'
cd $homedir

# calculating ase
asedest='results/ase.txt'
validdest='results/ase_validity.txt'
minreads=20
pseudo=50
python calculate_ase.py -home $homedir -asedest $asedest -validdest $validdest -minreads $minreads -pseudo $pseudo

# moving ase into the data folder
asedata='data/ase.txt'
validdata='data/ase_validity.txt'
mv $asedest $asedata 
mv $validdest $validdata

# calculating tf-ase correlations
minsamples=30
statdest=$outdir/test_statistics.txt
python calculate_tf_ase_correlations.py -home $homedir -asedata $asedata -validdata $validdata -minsamples $minsamples -statdest $statdest

# create output directory structure
mkdir $outdir
mkdir $outdir/fig

# fdr correction
statdata=$statdest
method='fdr_bh'
sigdest='results/test_statistics_fdr_bh.txt'
python multipletests_correction.py -home $homedir -statdata $statdata -method $method -sigdest $sigdest

# annotate ase loci
sigdata=$sigdest
annotdest=$outdir/significant_tf_ase_fdr_bh.txt
python extract_gene_names_from_ase_locus.py -home $homedir -sigdata $sigdata -annotdest $annotdest

# bonferroni correction
method='bonferroni'
sigdest='results/test_statistics_bonferroni.txt'
python multipletests_correction.py -home $homedir -statdata $statdata -method $method -sigdest $sigdest

# annotate ase loci
sigdata=$sigdest
annotdest=$outdir/significant_tf_ase_bonferroni.txt
python extract_gene_names_from_ase_locus.py -home $homedir -sigdata $sigdata -annotdest $annotdest

# create plots with bonferroni corrected tf-ase-pairs
sigdata=$outdir/significant_tf_ase_bonferroni.txt
plotdir=$outdir/fig
python generate_tf_ase_plot.py -home $homedir -sigdata $sigdata -plotdir $plotdir

# clean unnecessary files
rm results/test_statistics_fdr_bh.txt
rm results/test_statistics_bonferroni.txt

