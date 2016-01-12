#!/bin/bash

# Submit jobs. Assume 16 cores per proc. Make sure line numbers are the same.
for (( i = 1; i<= 16 ; i++ ))
do
    sed -i "3s/.*/#PBS -l size=$(($i*16)),walltime=24:00:00/" submits/pbs.sh
    qsub submits/pbs.sh
    #cat submits/pbs.sh && echo ''
done

