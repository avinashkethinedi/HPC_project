#!/bin/bash
#PBS -N hpc_data_compression
##PBS -l select=1:ncpus=4:mem=30gb:ngpus=1:gpu_model=v100:phase=21
##PBS -l select=1:ncpus=16:mem=60gb:ngpus=1:gpu_model=v100
#PBS -l select=4:ncpus=16:mem=60gb:mpiprocs=1:interconnect=fdr
#PBS -l walltime=01:00:00
#PBS -j oe
module purge
#module load courses/cpsc8200
module load openmpi/4.1.3-gcc/9.5.0-ucx
cd $PBS_O_WORKDIR
make clean;make
mpirun -n 4 ./SpMV -s 20