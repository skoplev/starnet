#!/bin/bash
# Submit 

module load rstudio

mkdir logs

bsub -J determinePower \
	-P acc_STARNET \
	-q alloc \
	-W 48:00 \
	-n 4 \
	-R rusage[mem=4000] \
	-e logs/error.%J \
	-o logs/output.%J \
	./determinePower.r
