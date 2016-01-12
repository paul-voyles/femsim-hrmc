#!/bin/bash

# This file will create directories for running models multiple times.
# The seeds are changed automatically but you may want to do them by hand.
# Copy this file above the hrmc directory then chmod and run it.
# You are expected to have the setup for your specific instance set up in hrmc/.
# You will have to change file names in this file to copy the correct ones.

# Make 10 directories, from 0 to 9.
for i in `seq 0 9`;
do
    echo r$i
    # If the run directory doesn't exist, make one.
    if [ ! -d "r$i" ]; then
        mkdir r$i
    fi
    # Copy in the paramfile and change the seed.
    cp hrmc/param_file.in r$i
    sed -i "11s/.*/$((10470+$i))           # seed/" r$i/param_file.in
    # Copy the eam file, fem file, and starting modelfile
    cp hrmc/PdSi.lammps.eam r$i
    cp hrmc/PdSi_1517atom_start.real.xyz r$i
    cp hrmc/fem_Pd82Si18_360C_1min.txt r$i
    # If the submit directory doesn't exist, create it
    if [ ! -d "r$i/submits" ]; then
        mkdir r$i/submits
    fi
    # Copy the submit files we need
    cp hrmc/submits/stampede_submit.py r$i/submits
    cp hrmc/submits/stampede_submit.sh r$i/submits
    # Symbolic link to hrmc/hrmc -- don't forget to make it!
    if [ ! -h "r$i/hrmc" ]; then
        ln -s /work/02916/maldonis/PdSi/t3/hrmc/hrmc r$i/hrmc
    fi
done

