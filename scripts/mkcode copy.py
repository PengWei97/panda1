import numpy as np

def writesh(outputname, threads, file,outputshfilename):
    shfile = """#!/bin/bash
# SBATCH --job-name=%s
# SBATCH --nodes=1
# SBATCH --ntasks-per-node=1
# SBATCH --cpus-per-task=8

source /PARA/app/scripts/cn-module.sh
module load intel-compilers/2018
module load cmake/3.13.4
module load gcc/5.4.0
module load fftw/3.3.4-gcc

# Set OMP_NUM_THREADS to the same value as -c
# with a fallback in case it isn't set.
# SLURM_CPUS_PER_TASK is set to the value of -c, but only if -c is explicitly set
if [ -n "$SLURM_CPUS_PER_TASK" ]; then
omp_threads=$SLURM_CPUS_PER_TASK
else
omp_threads=1
fi
export OMP_NUM_THREADS=$omp_threads

srun ./NiAl_omp_all 31 1000
    """ % (outputname)
    with open(outputshfilename, "w") as f:
        f.write(shfile)
    f.close()
    # outputname = "1.txt"

writesh( "1.txt", 1, 1,  "1.sh")

# print(shfile)