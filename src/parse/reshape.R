# Reshapes normalized expression matrices from multiple tissues into
# rows with genes and tissues, columns with samples (STARNET ID).
# Writes an RData object.

rm(list=ls())

library(reshape2)

setwd("/Users/sk/Google Drive/projects/cross-tissue")
source("src/parse.r")

# User input
data_dir = "/Users/sk/DataProjects/cross-tissue/STARNET/gene_exp_norm_batch"
# data_dir = "/sc/orga/projects/STARNET/koples01/expr-mat"

min_sd = 0.5  # sd threshold per tissue for including transcripts.

# Load data in folder, preprocess
expr_recast = loadNormData(data_dir,
	min_sd,
	exclude_files=c(
		"exp.GRCh38.GENCODE_r24.COR.mat",
		"exp.GRCh38.GENCODE_r24.FOC.mat",
		"exp.GRCh38.GENCODE_r24.MAC.mat"
	)
)

save(expr_recast, file="/Users/sk/DataProjects/cross-tissue/STARNET/gene_exp_norm_reshape/expr_recast.RData")
