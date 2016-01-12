#PBS -S /bin/bash
#PBS -A TG-DMR140035
#PBS -l size=224,walltime=24:00:00
#PBS -q batch

date
echo $PBS_O_WORKDIR
cd $PBS_O_WORKDIR
pwd
ls -l
aprun -n $PBS_NNODES ./hrmc $PBS_JOBID /lustre/medusa/maldonis/Al92Sm8/hrmc/param_file.2000000.in
