#!/bin/bash
# Submit job to server

module load rstudio

mkdir logs

bsub -J determinePower \
	-P acc_STARNET \
	-q premium \
	-W 80:00 \
	-n 1 \
	-R rusage[mem=40000] \
	-e logs/error.%J \
	-o logs/output.%J \
	./determinePower.R /hpc/users/koples01/links/STARNET/koples01/data/cross-tissue/gene_exp_norm_batch_imp/all.RData
