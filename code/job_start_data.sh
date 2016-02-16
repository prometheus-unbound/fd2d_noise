#!/bin/bash


brutus=`grep brutus start_inversion.m  | awk 'END {print NR}'`
euler=`grep euler start_inversion.m  | awk 'END {print NR}'`
monch=`grep monch start_inversion.m  | awk 'END {print NR}'`

if [ $brutus -gt 1 ] || [ $euler -gt 1 ]; then
    module load matlab/8.5
fi

if [ $monch -gt 1 ]; then
    module load matlab/r2015a
fi


while true; do
    read -p "Do conversion to mex-functions? " yn
    case $yn in
        [Yy]* ) cd mex_functions; rm -f run_*; matlab -nodisplay < compile_data.m; cd ..; break;;
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
if [ $brutus -gt 1 ] || [ $euler -gt 1 ]; then
	bsub -W "12:00" -R "rusage[mem=3072]" -o "logs/matlab_%J.out" -e "logs/matlab_%J.err" -n 1 matlab -nodisplay -singleCompThread -r start_inversion
fi



# START JOB ON MONCH
if [ $monch -gt 1 ]; then

cat <<EOF > data.sh
#!/bin/bash -l								

#SBATCH --partition=fichtner_compute
#SBATCH --job-name=data
#SBATCH --output=logs/matlab_%j.out
#SBATCH --error=logs/matlab_%j.err
#SBATCH --time=00-05:00:00
#SBATCH --ntasks=1
#SBATCH --mem=3072


######################
# Begin work section #
######################

module load matlab/r2015a
matlab -nodisplay -singleCompThread -r calculate_data

EOF

sbatch data.sh
rm data.sh

fi

