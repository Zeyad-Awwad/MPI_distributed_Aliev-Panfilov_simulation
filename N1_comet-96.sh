#!/bin/bash
#### this a bash script which we write out and then submit to the batch queue
# This script is intepreted by the Bourne Shell, sh
#
#SBATCH --account=csd562
#SBATCH --job-name=apf-strong-scale-96
#SBATCH --output="apf-strong-scale-96.out"
#SBATCH --partition="compute"
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=24
#SBATCH --export=ALL
####   #SBATCH --mail-type=BEGIN,END,FAIL
####   #SBATCH --mail-user=your_email@ucsd.edu
#SBATCH -t 00:03:00
# Commands go here

ibrun -np 96 ./apf -n 1800 -i 2000 -x 2 -y 48
ibrun -np 96 ./apf -n 1800 -i 2000 -x 4 -y 24
ibrun -np 96 ./apf -n 1800 -i 2000 -x 6 -y 16
ibrun -np 96 ./apf -n 1800 -i 2000 -x 8 -y 12
ibrun -np 96 ./apf -n 1800 -i 2000 -x 12 -y 8
ibrun -np 96 ./apf -n 1800 -i 2000 -x 16 -y 6
ibrun -np 96 ./apf -n 1800 -i 2000 -x 24 -y 4
ibrun -np 96 ./apf -n 1800 -i 2000 -x 48 -y 2 
ibrun -np 96 ./apf -n 1800 -i 2000 -x 2 -y 48 -k
ibrun -np 96 ./apf -n 1800 -i 2000 -x 4 -y 24 -k
ibrun -np 96 ./apf -n 1800 -i 2000 -x 6 -y 16 -k
ibrun -np 96 ./apf -n 1800 -i 2000 -x 8 -y 12 -k
ibrun -np 96 ./apf -n 1800 -i 2000 -x 12 -y 8 -k
ibrun -np 96 ./apf -n 1800 -i 2000 -x 16 -y 6 -k
ibrun -np 96 ./apf -n 1800 -i 2000 -x 24 -y 4 -k
ibrun -np 96 ./apf -n 1800 -i 2000 -x 48 -y 2 -k
