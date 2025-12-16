#! /bin/bash

source ~/miniconda3/etc/profile.d/conda.sh
conda activate trash
cd ~/TRASH
mkdir -p results
./TRASH_run.sh ~/TRASH/example_run/CP068268_39050443_39150442.fa --o ~/TRASH/results
