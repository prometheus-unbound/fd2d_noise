#!/bin/bash

module load matlab/r2015a

while true; do
    read -p "Do conversion to mex-functions? " yn
    case $yn in
        [Yy]* ) cd ../code/mex_functions; rm -f run_*; matlab -nodisplay < compile_all.m; cd ../../inversion; break;;
        [Nn]* ) break;;
        * ) echo "Please answer yes or no.";;
    esac
done


nf=`ls -l ../output/interferometry/G_2_* 2>/dev/null | grep -v ^l | wc -l`

if [ $nf -gt 0 ]; then
	while true; do
	    read -p "Move Green functions in output/interferometry? " yn
	    case $yn in
	        [Yy]* ) var=$(date '+%Y_%m_%d_%H_%M_%S'); mkdir ../../backup/$var; mv ../output/interferometry/G_2_*.mat ../../backup/$var/; break;;
	        [Nn]* ) break;;
	        * ) echo "Please answer yes or no.";;
	    esac
	done
fi


# START JOB ON EULER/BRUTUS
# bsub -W "12:00" -R "rusage[mem=3072]" -o "logs/matlab_%J.out" -e "logs/matlab_%J.err" -n 1 matlab -nodisplay -singleCompThread -r start_inversion


# START JOB ON MONCH
cat <<EOF > inversion.sh
#!/bin/bash -l								

#SBATCH --partition=fichtner_compute
#SBATCH --job-name=inversion
#SBATCH --output=logs/matlab_%j.out
#SBATCH --error=logs/matlab_%j.err
#SBATCH --time=23:59:00
#SBATCH --ntasks=1
#SBATCH --mem=8192


######################
# Begin work section #
######################

module load matlab/r2015a
matlab -nodisplay -singleCompThread -r start_inversion

EOF

sbatch inversion.sh
rm inversion.sh

