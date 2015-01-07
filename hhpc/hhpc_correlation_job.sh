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
python calculate_tf_ase_correlations.py -home $homedir -asedata $asedata -validdata $validdata -minsamples $minsamples

#outFile='results/out_'$(date +%Y-%m-%d_%H-%M-%S)'.txt'
#/home/asaha6/programs/anaconda/bin/python  calculate_tf_ase_correlations.py > $outFile


