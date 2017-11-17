#!/bin/sh
# #SBATCH --job-name=cyl_long
#SBATCH --account=trocks
# #SBATCH --account=nn9272k
# GI SOM ARGUMENT
##SBATCH --time='00:10:00'
#SBATCH --mem-per-cpu=3600M
# Number of tasks (MPI ranks): GI SOM ARGUMENT
##SBATCH --nodes=4
##SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=anjohan@uio.no
source /cluster/bin/jobsetup
module load intelmpi.intel
module load intel
mpirun -n $2 lmp_intel_cpu_intelmpi -in $1
