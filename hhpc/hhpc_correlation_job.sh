cd /home/asaha6/github/asenet
outFile='results/out_'$(date +%Y-%m-%d_%H-%M-%S)'.txt'
/home/asaha6/programs/anaconda/bin/python  calculate_tf_ase_correlations.py > $outFile
