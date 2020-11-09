#!/bin/bash
#### this a bash script which we write out and then submit to the batch queue
# This script is intepreted by the Bourne Shell, sh
#
#SBATCH --account=csd562
#SBATCH --job-name=apf-strong-scale-480
#SBATCH --output="apf-strong-scale-480.out"
#SBATCH --partition="compute"
#SBATCH --nodes=20
#SBATCH --ntasks-per-node=24
#SBATCH --export=ALL
####   #SBATCH --mail-type=BEGIN,END,FAIL
####   #SBATCH --mail-user=your_email@ucsd.edu
#SBATCH -t 00:03:00
# Commands go here
ibrun -np 480 ./apf -n 8000 -i 2000 -x 120 -y 4
ibrun -np 480 ./apf -n 8000 -i 2000 -x 96 -y 5
ibrun -np 480 ./apf -n 8000 -i 2000 -x 80 -y 6
ibrun -np 480 ./apf -n 8000 -i 2000 -x 60 -y 8
ibrun -np 480 ./apf -n 8000 -i 2000 -x 48 -y 10
ibrun -np 480 ./apf -n 8000 -i 2000 -x 120 -y 4 -k
ibrun -np 480 ./apf -n 8000 -i 2000 -x 96 -y 5 -k
ibrun -np 480 ./apf -n 8000 -i 2000 -x 80 -y 6 -k
ibrun -np 480 ./apf -n 8000 -i 2000 -x 60 -y 8 -k
ibrun -np 480 ./apf -n 8000 -i 2000 -x 48 -y 10 -k