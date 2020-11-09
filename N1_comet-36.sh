#!/bin/bash
#### this a bash script which we write out and then submit to the batch queue
# This script is intepreted by the Bourne Shell, sh
#
#SBATCH --account=csd562
#SBATCH --job-name=apf-strong-scale-36
#SBATCH --output="apf-strong-scale-36.out"
#SBATCH --partition="compute"
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=24
#SBATCH --export=ALL
####   #SBATCH --mail-type=BEGIN,END,FAIL
####   #SBATCH --mail-user=your_email@ucsd.edu
#SBATCH -t 00:03:00
# Commands go here

ibrun -np 36 ./apf -n 1800 -i 2000 -x 2 -y 18
ibrun -np 36 ./apf -n 1800 -i 2000 -x 3 -y 12
ibrun -np 36 ./apf -n 1800 -i 2000 -x 4 -y 9
ibrun -np 36 ./apf -n 1800 -i 2000 -x 6 -y 6
ibrun -np 36 ./apf -n 1800 -i 2000 -x 9 -y 4
ibrun -np 36 ./apf -n 1800 -i 2000 -x 12 -y 3
ibrun -np 36 ./apf -n 1800 -i 2000 -x 18 -y 2
ibrun -np 36 ./apf -n 1800 -i 2000 -x 2 -y 18 -k
ibrun -np 36 ./apf -n 1800 -i 2000 -x 3 -y 12 -k
ibrun -np 36 ./apf -n 1800 -i 2000 -x 4 -y 9 -k
ibrun -np 36 ./apf -n 1800 -i 2000 -x 6 -y 6 -k
ibrun -np 36 ./apf -n 1800 -i 2000 -x 9 -y 4 -k
ibrun -np 36 ./apf -n 1800 -i 2000 -x 12 -y 3 -k
ibrun -np 36 ./apf -n 1800 -i 2000 -x 18 -y 2 -k