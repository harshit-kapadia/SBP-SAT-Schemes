#!/usr/bin/env zsh
 
### Job name
#BSUB -J "MATLAB_ARRAY[15]"
 
### File / path where STDOUT will be written, the %J is the job id
#BSUB -o log_files/heated_cavity_DVM_%I
 
### Request the time you need for execution in minutes
### The format for the parameter is: [hour:]minute,
### that means for 80 minutes you could also use this: 1:20
#BSUB -W 50:00

### Request memory you need for your job in MB
#BSUB -M 10000
 
### Change to the work directory
cd /home/xx505837/SBP-SAT_Schemes/sbp_sat/DVM
 
### load modules and execute
module load MISC
module load matlab
  
# start non-interactive batch job
matlab -singleCompThread -nodisplay -nodesktop -nosplash -logfile log_files/heated_cavity_DVM_$LSB_JOBINDEX.log <<EOF
run ex_heated_cavity($LSB_JOBINDEX);
quit();
EOF
