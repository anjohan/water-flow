#!/bin/bash
#SBATCH --partition=normal
#SBATCH --ntasks=4
#SBATCH --gres=gpu:4
#SBATCH --job-name=sphere
echo $CUDA_VISIBLE_DEVICES
mpirun -n 4 /lammps/lammps_kokkos2/src/lmp_kokkos_cuda_mpi -k on g 4 -sf kk -pk kokkos newton on neigh half binsize 7.5 -in make_amorphous.in

