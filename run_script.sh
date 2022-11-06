#!/bin/bash
#PBS -N hpc_data_compression
##PBS -l select=1:ncpus=4:mem=30gb:ngpus=1:gpu_model=v100:phase=21
##PBS -l select=1:ncpus=16:mem=60gb:ngpus=1:gpu_model=v100
#PBS -l select=16:ncpus=16:mem=60gb:mpiprocs=16:interconnect=fdr
#PBS -l walltime=01:00:00
#PBS -j oe
module purge
#module load courses/cpsc8200
module load openmpi/4.1.3-gcc/9.5.0-ucx
cd $PBS_O_WORKDIR
make clean;make
mpirun -n 256 ./SpMV -s 20 >> log
mpirun -n 256 ./SpMV -s 21 >> log
mpirun -n 256 ./SpMV -s 22 >> log
mpirun -n 256 ./SpMV -s 23 >> log
mpirun -n 256 ./SpMV -s 24 >> log
mpirun -n 256 ./SpMV -s 25 >> log
mpirun -n 256 ./SpMV -s 26 >> log
mpirun -n 256 ./SpMV -s 27 >> log
mpirun -n 256 ./SpMV -s 28 >> log
mpirun -n 256 ./SpMV -s 29 >> log
mpirun -n 256 ./SpMV -s 30 >> log