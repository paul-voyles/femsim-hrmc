#!/bin/bash

# After running to_copy.sh, copy this file above hrmc/ as well
# and then run it. It will python submit (slurm) all the directories.
# You are expected to have the setup for your specific instance set
# up in hrmc/.

for i in `seq 0 9`;
do
    echo r$i
    cd r$i
    python submits/stampede_submit.py
    cd ..
done

