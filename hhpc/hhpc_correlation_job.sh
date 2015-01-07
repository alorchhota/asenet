homedir='/home/asaha6/github/asenet'
cd $homedir

# calculating ase
asedest='results/ase.txt'
validdest='results/ase_validity.txt'
minreads=20
pseudo=20
#python calculate_ase.py -home $homedir -asedest $asedest -validdest $validdest -minreads $minreads -pseudo $pseudo

# moving ase into the data folder
asedata='data/ase.txt'
validdata='data/ase_validity.txt'
mv $asedest $asedata 
mv $validdest $validdata

# calculating tf-ase correlations
minsamples=30
statdest='results/test_statistics.txt'
#python calculate_tf_ase_correlations.py -home $homedir -asedata $asedata -validdata $validdata -minsamples $minsamples -statdest $statdest

# fdr correction
statdata=$statdest
method='fdr_bh'
sigdest='results/test_statistics_fdr_bh.txt'
python multipletests_correction.py -home $homedir -statdata $statdata -method $method -sigdest $sigdest

# bonferroni correction
method='bonferroni'
sigdest='results/test_statistics_bonferroni.txt'
python multipletests_correction.py -home $homedir -statdata $statdata -method $method -sigdest $sigdest

#outFile='results/out_'$(date +%Y-%m-%d_%H-%M-%S)'.txt'
#/home/asaha6/programs/anaconda/bin/python  calculate_tf_ase_correlations.py > $outFile


