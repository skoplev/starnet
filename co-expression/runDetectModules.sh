#!/bin/bash
# Submit job to server

module load rstudio

mkdir logs

bsub -J findModules \
	-P acc_STARNET \
	-q alloc \
	-W 48:00 \
	-n 1 \
	-R rusage[mem=100000] \
	-e logs/error.%J \
	-o logs/output.%J \
	./detectModules.R between_within

# 	./detectModules.R complete
# 	./detectModules.R single
