#!/bin/bash
#### this a bash script which we write out and then submit to the batch queue
# This script is intepreted by the Bourne Shell, sh
#
#SBATCH --account=csd562
#SBATCH --job-name=apf-strong-scale-384
#SBATCH --output="apf-strong-scale-384.out"
#SBATCH --partition="compute"
#SBATCH --nodes=16
#SBATCH --ntasks-per-node=24
#SBATCH --export=ALL
####   #SBATCH --mail-type=BEGIN,END,FAIL
####   #SBATCH --mail-user=your_email@ucsd.edu
#SBATCH -t 00:03:00
# Commands go here

ibrun -np 384 ./apf -n 8000 -i 2000 -x 16 -y 24
ibrun -np 384 ./apf -n 8000 -i 2000 -x 24 -y 16
ibrun -np 384 ./apf -n 8000 -i 2000 -x 32 -y 12
ibrun -np 384 ./apf -n 8000 -i 2000 -x 16 -y 24 -k
ibrun -np 384 ./apf -n 8000 -i 2000 -x 24 -y 16 -k
ibrun -np 384 ./apf -n 8000 -i 2000 -x 32 -y 12 -k