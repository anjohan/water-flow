./lag_passivering.sh
cd data
sbatch --ntasks-per-node=1 --nodes=1 --time="00:20:00" --job-name aji0pass job.sh passivering.in 1
cd ..
