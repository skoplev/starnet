#!/bin/bash
# Submit job to server

module load rstudio

mkdir logs

bsub -J determinePower \
	-P acc_STARNET \
	-q alloc \
	-W 1:00 \
	-n 1 \
	-R rusage[mem=8000] \
	-e logs/error.%J \
	-o logs/output.%J \
	./detectModules.R
