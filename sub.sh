#!/bin/bash

#SBATCH --job-name=aka_stream_test0
#SBATCH --time=23:59:59
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40

#SBATCH --error=%x-%j.err
#SBATCH --output=%x-%j.out

#SBATCH --mail-user=ywppku@gmail.com
#SBATCH --mail-type=BEGIN,FAIL,END

#PROGRAM=./smilei_p3.sh
PROGRAM=/home/a/anticipa/weipeng/CODES/AKA52/aka.exe
SOURCE=/home/a/anticipa/weipeng/CODES/AKA52/compile_aka52_niagara.sh
NPROC=32 # Note: nodes x ntasks-per-node(=40)
NAMELIST=/home/a/anticipa/weipeng/CODES/AKA52/input/
source $SOURCE
mpirun -n $NPROC $PROGRAM  $NAMELIST
