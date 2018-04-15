#!/bin/bash
## Project:
#SBATCH --account=nn9272k
## Job name:
#SBATCH --job-name=pillar
## Wall time limit:
#SBATCH --time=48:0:0
## Number of nodes:
#SBATCH --nodes=8
## Number of tasks to start on each node:
#SBATCH --ntasks-per-node=32
## Set OMP_NUM_THREADS
# #SBATCH --cpus-per-task=16

## Recommended safety settings:
set -o errexit # Make bash exit on any error
set -o nounset # Treat unset variables as errors

## Software modules
module restore system   # Restore loaded modules to the default
#module load LAMMPS/11Aug17-foss-2017a

## Prepare input files
#cp * $SCRATCH
#cd $SCRATCH

## Make sure output is copied back after job finishes
#copyback log.lammps

## Run the application
module load intel/2017b
module load tbb/2017.6.196

mpirun /nird/home/anderhaf/lammps_builds/lammps_intel/src/lmp_intel_cpu_intelmpi -in vann.in
