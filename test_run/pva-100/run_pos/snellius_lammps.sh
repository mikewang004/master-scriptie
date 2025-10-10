#!/bin/bash
#SBATCH -p rome
#SBATCH --nodes=2
#SBATCH --ntasks=192
#SBATCH --cpus-per-task=1
#SBATCH --mem=20G
#SBATCH --time=24:00:00

module load 2022
module load OpenMPI/4.1.4-NVHPC-22.7-CUDA-11.7.0
module load LAMMPS/23Jun2022-foss-2022a-kokkos
mpirun lmp -in in.cooling




