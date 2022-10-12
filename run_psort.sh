#!/bin/bash
#SBATCH --job-name=your_job_name # Job name
#SBATCH --nodes=1                    # Run all processes on a single node	
#SBATCH --ntasks=20                   # number of processes = 20
#SBATCH --cpus-per-task=1            # Number of CPU cores allocated to each process
#SBATCH --partition=Project            # Partition name: Project or Debug (Debug is default)

mpirun -np 4 ./psort 20 ./test_data/20a.in
./check 20 ./test_data/20a.in.parallel.out

# mpirun -np 8 ./psort 80 ./test_data/80a.in
# ./check 80 ./test_data/80a.in.parallel.out
# mpirun -np 8 ./psort 10000 ./test_data/10000a.in
# ./check 10000 ./test_data/10000a.in.parallel.out