#qlogin --account=nn9272k --nodes=$1 --ntasks-per-node=$2
#qlogin --account=trocks --nodes=$1 --ntasks-per-node=$2 --time=$3:00
#qlogin --account=trocks --nodes=$1 --ntasks-per-node=$2 --cpus-per-task=$3 --time=4:0:0
qlogin --account=trocks --nodes=16 --ntasks-per-node=16 --time=4:0:0
