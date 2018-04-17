# Number of tasks (MPI ranks):
#SBATCH --time='02:00:00'
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=anjohan@uio.no
#SBATCH --mem-per-cpu=3600M
source /cluster/bin/jobsetup
module load intelmpi.intel
module load intel

make data/amorphous.data
