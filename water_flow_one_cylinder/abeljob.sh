#!/bin/sh
#SBATCH --time='1-00:00:00'
#SBATCH --account=nn9272k
#SBATCH --nodes=32
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=anjohan@uio.no
#SBATCH --mem-per-cpu=3600M

source /cluster/bin/jobsetup
module load intel/2018.1
module load intelmpi.intel

$HOME/miniconda3/bin/make ITERATIONS=20
