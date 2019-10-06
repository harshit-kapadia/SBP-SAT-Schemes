#!/usr/local_rwth/bin/zsh
 
### Job name
#SBATCH --job-name=MATLAB_SERIAL
 
### File / path where STDOUT will be written, the %J is the job id
#SBATCH --output=log_files_new/lid_driven_cavity_n20_M12
 
### Request the time you need for execution. The full format is D-HH:MM:SS
### You must at least specify minutes or days and hours and may add or
### leave out any other parameters
#SBATCH --time=1000
 
### Request the memory you need for your job. You can specify this
### in either MB (1024M) or GB (4G).
#SBATCH --mem-per-cpu=5G
  
### Change to the work directory, if not in submit directory
#cd /home/xx505837/git-repo/SBP-SAT_Schemes/sbp_sat
 
### Load the required modules
module load MISC
module load matlab
 
 
# start non-interactive batch job
matlab -singleCompThread -nodisplay -nodesktop -nosplash -logfile log_files_new/lid_driven_cavity_n20_M12.txt <<EOF
run ex_lid_driven_cavity(12);
quit();
EOF
