#!/usr/bin/env Rscript

# Determines power of cross-tissue modules

rm(list=ls())

library(data.table)
library(reshape2)
library(plyr)
library(parallel)
library(RColorBrewer)
library(Matrix)

library(compiler)
enableJIT(3)

# Load and setup WGCNA
library(WGCNA)
options(stringsAsFactors = FALSE)

enableWGCNAThreads(nThreads=8)  # assuming 2 threads per core


# # Filter out low connections
# # k = cmat[cmat > 0.05]
# # scaleFreeFitIndex(k, nBreaks=100)
# # k = c(k1, k2)
# fastScaleFreeFitIndex = function(k, nBreaks=50, epsilon=1e-09) {
# 	# k = as.vector(k)
# 	# Cut connections into discrete groups
# 	# groups = cut(k, nBreaks)
# 	# group_mean = tapply(k, groups, mean)  # mean connectivity in each group

# 	# hist is faster than cut
# 	groups = hist(k, breaks=nBreaks, plot=FALSE)

# 	# groups = tabulate(.bincode(k, breaks=nBreaks))

# 	groups$pk = groups$counts / length(k)
# 	groups$log_pk = log(groups$pk + epsilon)
# 	groups$log_mids = log(groups$mids)

# 	# Count number of connections in each group
# 	# group_pk = table(groups) / length(k)

# 	# fit = fastLm(log(group_pk)~log(group_mean))
# 	# fit = lm(log(groups$pk) ~ log(groups$mids))

# 	fit = lm(groups$log_pk ~ groups$log_mids)

# 	out = data.frame(Rsquared.SFT=summary(fit)$r.squared)
# 	return(out)

# 	plot(groups$log_mids, groups$log_pk)
# }


# setwd("/Users/sk/Google Drive/projects/cross-tissue")

data_dir = "/Users/sk/DataProjects/cross-tissue"

setwd("/Users/sk/Google Drive/projects/cross-tissue/co-expression/determinePower")
# setwd("/sc/orga/projects/STARNET/koples01/cross-tissue/co-expression/determinePower")

source("../../src/base.R")


# # STARNET phenotype data
# pheno = fread(file.path(
# 	"/Volumes/SANDY/phenotype_data",
# 	"STARNET_main_phenotype_table.cases.Feb_29_2016.tbl"
# ))

# Load imputed recast gene expression matrix
load(file.path(data_dir, "STARNET/gene_exp_norm_batch_imp/all.RData"), verbose=TRUE)

if (any(is.na(mat))) {
	stop("Gene expression matrix contains missing data.")
}

# Standardize expression data
# mat_scaled = scale(t(mat))

# # Match phenotype data to selected gene expression matrix
# pheno_matched = pheno[match(colnames(mat), pheno$starnet.ID), ]

# Picking a soft power threshold to get scale-free correlation networks from each tissue alone
# powers = seq(1, 10, length.out=20)

powers = seq(0.1, 10, length.out=51)

# powers = seq(1, 5, length.out=50)
# powers = seq(0.1, 10, length.out=50)
# block_size = 4000

n_breaks = 50  # binning of connectivity

# abs_cor_min = 0.05  # minimum correlation to be considered
# abs_cor_min = 0.1  # minimum correlation to be considered
# abs_cor_min = 0.15  # minimum correlation to be considered
# abs_cor_min = 0.2  # minimum correlation to be considered
# abs_cor_min = 0.3  # minimum correlation to be considered
cor_quant = 0.95


# Evaluate sequence of soft network cutoffs
con_eval = list()
for (i in 1:length(unique(row_meta$tissue))) {
	tissue = unique(row_meta$tissue)[i]
	message("Estimating power for: ", tissue)

	# Symmetric correlation matrix within-tissue
	message("estimating correlation coefficients")
	cmat = corFast(t(mat[row_meta$tissue == tissue, ]),
		nThreads=8  # parallel computation not working
	)
	gc()
	message("dim: ", nrow(cmat), ", ", ncol(cmat))

	# Convert to adjacency matrix
	message("making adjacency matrix")
	cmat = abs(cmat)
	gc()

	abs_cor_min = quantile(cmat, cor_quant)
	gc()

	# Drop entries below threshold
	message("converting to sparse matrix, abs_cor_min: ", abs_cor_min)
	cmat[cmat < abs_cor_min] = 0  # delete entries
	gc()

	# Convert to sparse, symmetric matrix, for increased efficiency
	cmat = forceSymmetric(cmat)
	gc()

	cmat = Matrix(cmat, sparse=TRUE)
	gc()

	message("non-zero coefficients: ", nnzero(cmat))
	gc()

	# Returns scale-free distribution fits and connectivity vectors for each tissue
	# at the power series.
	fits = lapply(powers, function(pow) {
		message("power: ", pow)
		# message("\tadjusting weights...")
		cmat_power = cmat^pow
		gc()
		# message("\tcalculating connectivity")
		# Calculate connectivity for nodes in both tissues
		k = apply(cmat_power, 1, sum) - 1  # subtracting diagonal
		gc()

		# message("\tfitting scale-free connectivity model")
		fit = scaleFreeFitIndex(k, nBreaks=n_breaks, removeFirst=TRUE)
		gc()
		return(list(fit=fit, k=k))
	})

	con_eval[[i]] = fits
}


# Test optimal beta for combinations of tissues. Uses cross-correlations only.
con_eval_pairs = list()
paired_tissue = combn(unique(row_meta$tissue), 2)
for (i in 1:ncol(paired_tissue)) {
	# Cross-correlation only
	message("calculating cross-tissue correlation: ", paired_tissue[1, i], ", ", paired_tissue[2, i])

	message("estimating correlation coefficients")
	cmat = corFast(
		t(mat[row_meta$tissue == paired_tissue[1, i], ]),
		t(mat[row_meta$tissue == paired_tissue[2, i], ]),
		nThreads=8  # might not work?
	)
	gc()  # garbage collection

	# Convert to adjacency matrix
	message("making adjacency matrix")
	cmat = abs(cmat)
	gc()

	abs_cor_min = quantile(cmat, cor_quant)
	gc()

	# Drop entries below threshold
	message("converting to sparse matrix, abs_cor_min: ", abs_cor_min)
	cmat[cmat < abs_cor_min] = 0
	gc()

	# Convert to sparse matrix, for increased efficiency
	cmat = Matrix(cmat, sparse=TRUE)
	message("non-zero coefficients: ", nnzero(cmat))
	gc()

	# Returns scale-free distribution fits and connectivity vectors for each tissue
	# at the power series.
	fits = lapply(powers, function(pow) {
		message("power: ", pow)
		# message("\tadjusting weights...")
		cmat_power = cmat^pow
		gc()
		# message("\tcalculating connectivity")
		# Calculate connectivity for nodes in both tissues
		k1 = apply(cmat_power, 1, sum)
		k2 = apply(cmat_power, 2, sum)

		# message("\tfitting scale-free connectivity model")
		fit = scaleFreeFitIndex(c(k1, k2), nBreaks=n_breaks, removeFirst=TRUE)
		gc()
		return(list(fit=fit, k1=k1, k2=k2))
	})

	con_eval_pairs[[i]] = fits
	gc()
}
names(con_eval_pairs) = apply(paired_tissue, 2, paste, collapse="_")

dir.create("output")
save(con_eval, file="output/thresh_eval.RData")
save(con_eval_pairs, file="output/thresh_eval_pairs.RData")

