#!/bin/bash -x
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=64
#SBATCH --time=00:30:00
#SBATCH --job_name=JG_LATTICE_BOLTZMANN_PERFORMANCE_TEST
#SBATCH --mem=249gb
#SBATCH --export=ALL
#SBATCH --partition=dev_cpu_il
#SBATCH --output=test_results.%j.out
#SBATCH --error=test_results.%j.err
module load compiler/intel
module load mpi/impi

export FOR_COARRAY_CONFIG_FILE=/home/fr/fr_fr/fr_greiandr/FORTRAN/hello.conf
export EXE=/home/fr/fr_fr/fr_greiandr/FORTRAN/hello

echo "-genvall -genv I_MPI_FABRICS=shm:ofi -n "$SLURM_NTASKS" "$EXE > $FOR_COARRAY_CONFIG_FILE


echo "Running on ${SLURM_JOB_NUM_NODES} nodes with ${SLURM_JOB_CPUS_PER_NODE} cores each."
echo "Each node has ${SLURM_MEM_PER_NODE} of memory allocated to this job."
time $EXE > ./hello.txt 