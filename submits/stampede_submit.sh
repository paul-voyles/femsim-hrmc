#!/bin/sh

#SBATCH --job-name=hrmc                  # job name
#SBATCH --partition=normal                # default "univ" if not specified
#SBATCH --error=job.%J.err              # error file
#SBATCH --output=job.%J.out             # output file

#SBATCH --time=2-00:00:00               # run time in days-hh:mm:ss

#SBATCH --nodes=7                      # number of nodes requested (n)
#SBATCH --ntasks=112                    # required number of CPUs (n)
#SBATCH --ntasks-per-node=16             # default 16 (Set to 1 for OMP)
#SBATCH --cpus-per-task=1              # default 1 (Set to 16 for OMP)
##SBATCH --mem=16384                    # total RAM in MB, max 64GB  per node
##SBATCH --mem-per-cpu=4000              # RAM in MB (default 4GB, max 8GB)

##SBATCH --export=ALL

#SBATCH --mail-user=jjmaldonis@gmail.com
#SBATCH --mail-type=begin
#SBATCH --mail-type=end


echo "Date:"
date
echo "Github has:"
git rev-parse --verify HEAD
echo "Using Stampede."
echo "JobID = $SLURM_JOB_ID"
echo "Using $SLURM_NNODES nodes"
echo "Using $SLURM_NODELIST nodes."
echo "Number of cores per node: $SLURM_TASKS_PER_NODE"
echo "Submit directory: $SLURM_SUBMIT_DIR"
echo ""
cat $@

# Executable
ibrun ./hrmc $SLURM_JOB_ID $@

echo "Finished on:"
date
