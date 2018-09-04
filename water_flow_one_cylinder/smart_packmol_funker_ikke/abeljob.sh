#!/bin/sh
# Number of tasks (MPI ranks):
#SBATCH --time='01:00:00'
#SBATCH --account=nn9272k
#SBATCH --nodes=4
# #SBATCH --account=trocks
# #SBATCH --nodes=5
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=anjohan@uio.no
#SBATCH --mem-per-cpu=3600M
source /cluster/bin/jobsetup
module load intelmpi.intel
module load intel
module load python3

make data/add_water.inp