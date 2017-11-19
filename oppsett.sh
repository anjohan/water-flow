#!/bin/bash

set -e -x
./lag_oppsett.sh

cd data
echo "hei fra $(pwd)"
sbatch --ntasks-per-node=16 --nodes=8 --time="00:30:00" --job-name ajoppsett job.sh oppsett.in 128
