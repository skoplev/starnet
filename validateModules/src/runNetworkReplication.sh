#!/bin/bash
# Submit job to server

module load rstudio

mkdir logs

# Worst time calculation for runtime at M=1,000 permutations based on test run for block 1 with M=50:
# ((9252 / 60 / 60) / 50) * 1000 = 51.4h
# 
# memory is per core, MB unit
# n is number of cores
# Walltime in hours:minutes
#

for i in {1..5}
do
	bsub -J NetRep \
		-P acc_STARNET \
		-q premium \
		-W 100:00 \
		-n 2 \
		-R rusage[mem=128000] \
		-e logs/error.%J \
		-o logs/output.%J \
		./src/networkReplication.R $i
done
