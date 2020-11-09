#!/bin/bash
#### this a bash script which we write out and then submit to the batch queue
# This script is intepreted by the Bourne Shell, sh
#
#SBATCH --account=csd562
#SBATCH --job-name=apf-strong-scale-60
#SBATCH --output="apf-strong-scale-60.out"
#SBATCH --partition="compute"
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=24
#SBATCH --export=ALL
####   #SBATCH --mail-type=BEGIN,END,FAIL
####   #SBATCH --mail-user=your_email@ucsd.edu
#SBATCH -t 00:03:00
# Commands go here

ibrun -np 60 ./apf -n 1800 -i 2000 -x 2 -y 30
ibrun -np 60 ./apf -n 1800 -i 2000 -x 3 -y 20
ibrun -np 60 ./apf -n 1800 -i 2000 -x 4 -y 15
ibrun -np 60 ./apf -n 1800 -i 2000 -x 5 -y 12
ibrun -np 60 ./apf -n 1800 -i 2000 -x 6 -y 10
ibrun -np 60 ./apf -n 1800 -i 2000 -x 10 -y 6
ibrun -np 60 ./apf -n 1800 -i 2000 -x 12 -y 5
ibrun -np 60 ./apf -n 1800 -i 2000 -x 15 -y 4
ibrun -np 60 ./apf -n 1800 -i 2000 -x 20 -y 3
ibrun -np 60 ./apf -n 1800 -i 2000 -x 30 -y 2



ibrun -np 60 ./apf -n 1800 -i 2000 -x 2 -y 30 -k
ibrun -np 60 ./apf -n 1800 -i 2000 -x 3 -y 20 -k
ibrun -np 60 ./apf -n 1800 -i 2000 -x 4 -y 15 -k
ibrun -np 60 ./apf -n 1800 -i 2000 -x 5 -y 12 -k
ibrun -np 60 ./apf -n 1800 -i 2000 -x 6 -y 10 -k
ibrun -np 60 ./apf -n 1800 -i 2000 -x 10 -y 6 -k
ibrun -np 60 ./apf -n 1800 -i 2000 -x 12 -y 5 -k
ibrun -np 60 ./apf -n 1800 -i 2000 -x 15 -y 4 -k
ibrun -np 60 ./apf -n 1800 -i 2000 -x 20 -y 3 -k
ibrun -np 60 ./apf -n 1800 -i 2000 -x 30 -y 2 -k