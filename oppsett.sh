#!/bin/bash

set -e -x

cd data
echo "hei fra $(pwd)"
sbatch --ntasks-per-node=16 --nodes=20 --time="10:00:00" --job-name ajoppsett job.sh oppsett.in 320
