homedir='/home/asaha6/github/asenet'
cd $homedir

asedest='results/ase.txt'
validdest='results/ase_validity.txt'
minreads=20
pseudo=20
python calculate_ase.py -home $homedir -asedest $asedest -validdest $validdest -minreads $minreads -pseudo $pseudo

asedata='data/ase.txt'
validdata='data/ase_validity.txt'
mv $asedest $asedata 
mv $validdest $validdata

#outFile='results/out_'$(date +%Y-%m-%d_%H-%M-%S)'.txt'
#/home/asaha6/programs/anaconda/bin/python  calculate_tf_ase_correlations.py > $outFile

