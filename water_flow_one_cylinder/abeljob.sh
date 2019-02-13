#!/bin/sh
# Number of tasks (MPI ranks):
#SBATCH --time='7-00:00:00'
# #SBATCH --account=nn9272k
#SBATCH --nodes=20
#SBATCH --account=trocks
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=anjohan@uio.no
#SBATCH --mem-per-cpu=3600M

source /cluster/bin/jobsetup
module load intel/2019.2
module load intelmpi.intel

mpirun -n 1 /usit/abel/u1/anjohan/cmake_lammps/lammps/build/lmp -in make_amorphous.in
