#!/bin/bash


brutus=`grep brutus start_eigenvalues.m  | awk 'END {print NR}'`
euler=`grep euler start_eigenvalues.m  | awk 'END {print NR}'`
monch=`grep monch start_eigenvalues.m  | awk 'END {print NR}'`

if [ $brutus -gt 1 ] || [ $euler -gt 1 ]; then
    module load matlab/8.5
fi

if [ $monch -gt 1 ]; then
    module load matlab/r2015a
fi



while true; do
    read -p "Do conversion to mex-functions? " yn
    case $yn in
        [Yy]* ) cd ../code/mex_functions; rm -f run_*; matlab -nodisplay < compile_all.m; cd ../../inversion; break;;
        [Nn]* ) break;;
        * ) echo "Please answer yes or no.";;
    esac
done



# START JOB ON EULER/BRUTUS
if [ $brutus -gt 1 ] || [ $euler -gt 1 ]; then
	bsub -W "120:00" -R "rusage[mem=10240]" -o "logs/matlab_%J.out" -e "logs/matlab_%J.err" -n 1 matlab -nodisplay -singleCompThread -r start_eigenvalues
fi



# START JOB ON MONCH
if [ $monch -gt 1 ]; then

cat <<EOF > eigenvalues.sh
#!/bin/bash -l								

#SBATCH --partition=fichtner_compute_wk
#SBATCH --job-name=inversion
#SBATCH --output=logs/matlab_%j.out
#SBATCH --error=logs/matlab_%j.err
#SBATCH --time=07-00:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=32768

######################
# Begin work section #
######################

module load matlab/r2015a
matlab -nodisplay -r start_eigenvalues

EOF

sbatch eigenvalues.sh
rm eigenvalues.sh

fi

