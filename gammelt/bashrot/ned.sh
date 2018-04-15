#!/bin/bash

set -e -x

cd data
for i in $(cat ../i.txt)
do
    cd ilik$i
    echo "hei fra $(pwd)"
    sbatch --ntasks-per-node=16 --nodes=20 --time="100:00:00" --job-name aji$i job.sh ned.in 320
    cd ..
done
