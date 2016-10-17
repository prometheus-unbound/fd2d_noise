#!/bin/bash


brutus=`grep brutus calculate_hessian.m  | awk 'END {print NR}'`
euler=`grep euler calculate_hessian.m  | awk 'END {print NR}'`
monch=`grep monch calculate_hessian.m  | awk 'END {print NR}'`

if [ $brutus -gt 1 ] || [ $euler -gt 1 ]; then
    module load matlab/8.5
fi

if [ $monch -gt 1 ]; then
    module load matlab/r2015a
fi


while true; do
    read -p "Do conversion to mex-functions? " yn
    case $yn in
        [Yy]* ) cd mex_functions; rm -f run_*; matlab -nodisplay < compile_all.m; cd ..; break;;
        [Nn]* ) break;;
        * ) echo "Please answer yes or no.";;
    esac
done


# START JOB ON EULER/BRUTUS
if [ $brutus -gt 1 ] || [ $euler -gt 1 ]; then
	bsub -W "12:00" -R "rusage[mem=3072]" -o "logs/matlab_%J.out" -e "logs/matlab_%J.err" -n 1 matlab -nodisplay -singleCompThread -r calculate_hessian
fi



# START JOB ON MONCH
if [ $monch -gt 1 ]; then

cat <<EOF > hessian.sh
#!/bin/bash -l								

#SBATCH --partition=other_hugemem
#SBATCH --job-name=hessian
#SBATCH --output=logs/matlab_%j.out
#SBATCH --error=logs/matlab_%j.err
#SBATCH --time=01-00:00:00
#SBATCH --ntasks=1
#SBATCH --mem=65536


######################
# Begin work section #
######################

module load matlab/r2015a
matlab -nodisplay -r calculate_hessian

EOF

sbatch hessian.sh
rm hessian.sh

fi

