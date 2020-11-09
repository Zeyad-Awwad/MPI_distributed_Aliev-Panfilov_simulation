#!/bin/bash
#### this a bash script which we write out and then submit to the batch queue
# This script is intepreted by the Bourne Shell, sh
#
#SBATCH --account=csd562
#SBATCH --job-name=apf-strong-scale-240
#SBATCH --output="apf-strong-scale-240.out"
#SBATCH --partition="compute"
#SBATCH --nodes=10
#SBATCH --ntasks-per-node=24
#SBATCH --export=ALL
####   #SBATCH --mail-type=BEGIN,END,FAIL
####   #SBATCH --mail-user=your_email@ucsd.edu
#SBATCH -t 00:03:00
# Commands go here
ibrun -np 240 ./apf -n 8000 -i 2000 -x 30 -y 8
ibrun -np 240 ./apf -n 8000 -i 2000 -x 40 -y 6
ibrun -np 240 ./apf -n 8000 -i 2000 -x 48 -y 5
ibrun -np 240 ./apf -n 8000 -i 2000 -x 30 -y 8 -k
ibrun -np 240 ./apf -n 8000 -i 2000 -x 40 -y 6 -k
ibrun -np 240 ./apf -n 8000 -i 2000 -x 48 -y 5 -k
