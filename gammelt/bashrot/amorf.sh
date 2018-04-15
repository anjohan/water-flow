#!/bin/bash

set -e -x
./lag_amorf.sh

cd data
echo "hei fra $(pwd)"
sbatch --ntasks-per-node=16 --nodes=8 --time="02:30:00" --job-name ajamorf job.sh amorf.in 128
