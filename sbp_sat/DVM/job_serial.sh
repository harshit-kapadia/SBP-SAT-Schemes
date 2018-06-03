#!/usr/bin/env zsh
 
### Job name
#BSUB -J "MATLAB_ARRAY[15]"
 
### File / path where STDOUT will be written, the %J is the job id
#BSUB -o log_files/solving_Inflow2D_DVM%I
 
### Request the time you need for execution in minutes
### The format for the parameter is: [hour:]minute,
### that means for 80 minutes you could also use this: 1:20
#BSUB -W 1:20
 
### Request memory you need for your job in MB
#BSUB -M 3000
 
### Change to the work directory
cd /home/ns179556/SBP_1D/
 
### load modules and execute
module load MISC
module load matlab
 
 
# start non-interactive batch job
matlab -singleCompThread -nodisplay -nodesktop -nosplash -logfile log_files/solver_HC2D_$LSB_JOBINDEX.log <<EOF
run solve_Inflow2D_DVM($LSB_JOBINDEX);
quit();
EOF
