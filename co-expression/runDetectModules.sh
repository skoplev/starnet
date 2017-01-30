#!/bin/bash
# Submit job to server

module load rstudio

mkdir logs

bsub -J determinePower \
	-P acc_STARNET \
	-q alloc \
	-W 80:00 \
	-n 1 \
	-R rusage[mem=72000] \
	-e logs/error.%J \
	-o logs/output.%J \
	./detectModules.R /hpc/users/koples01/links/STARNET/koples01/data/cross-tissue/gene_exp_norm_batch_imp/all.RData
