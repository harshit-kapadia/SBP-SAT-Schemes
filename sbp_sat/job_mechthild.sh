#!/bin/bash
# Your job name.
#SBATCH -J lid_driven_cavity
#
# Output files. .out. is for standard output .err. is for the error output. 
#SBATCH -o %x-%j.out-%N # name of the stdout, using the job name (%x), the job number (%j) and the first node (%N)
#
# Maximum expected runtime.  ( 00 Days, 1 hour, 00 minutes, 00 seconds) 
#SBATCH --time=00-10:00:00   
#
# Allocate one node with all 16 CPU cores 
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16 
#
# Choose Partition (Queue) 
#SBATCH --partition long 
#
# Choose Constraints (node type, here restricting to standard nodes)
#SBATCH --constraint=RAM192
#
# Mail Options 
#SBATCH --mail-type=FAIL,BEGIN,END    # An email is sent on begin, end, and failure of the job 
#SBATCH --mail-user=sarnam@mpi-magdeburg.mpg.de   # E-Mail for notification
#
### END OF THE SLURM SPECIFIC PART ###

# Setup OpenMP 
if [ -n "$SLURM_CPUS_PER_TASK" ]; then
  omp_threads=$SLURM_CPUS_PER_TASK
else
  omp_threads=1
fi
export OMP_NUM_THREADS=$omp_threads

# set up the environment (choose whats needed)
# load the modules system 
source /etc/profile.d/modules.sh


# Load MATLAB

module load apps/matlab/2018a

# Or one of the older versions 
#module load apps/matlab/2012b
#module load apps/matlab/2016b


# SETUP Your MATLAB environment 
# The default search path for m-files.
export MATLABPATH=.


#run your script
matlab -nojvm -nosplash -r "test4_grid_refinement(12,95)" < /dev/null
