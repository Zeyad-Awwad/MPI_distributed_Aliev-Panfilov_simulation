#!/bin/bash
#### this a bash script which we write out and then submit to the batch queue
# This script is intepreted by the Bourne Shell, sh
#
#SBATCH --account=csd562
#SBATCH --job-name=apf-strong-scale-96
#SBATCH --output="apf-strong-scale-96-8000.out"
#SBATCH --partition="compute"
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=24
#SBATCH --export=ALL
####   #SBATCH --mail-type=BEGIN,END,FAIL
####   #SBATCH --mail-user=your_email@ucsd.edu
#SBATCH -t 00:03:00
# Commands go here
ibrun -np 96 ./apf -n 8000 -i 2000 -x 16 -y 6
ibrun -np 96 ./apf -n 8000 -i 2000 -x 24 -y 4
ibrun -np 96 ./apf -n 8000 -i 2000 -x 48 -y 2 
ibrun -np 96 ./apf -n 8000 -i 2000 -x 16 -y 6 -k
ibrun -np 96 ./apf -n 8000 -i 2000 -x 24 -y 4 -k
ibrun -np 96 ./apf -n 8000 -i 2000 -x 48 -y 2 -k
