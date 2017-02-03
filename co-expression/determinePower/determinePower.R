#!/usr/bin/env Rscript

# Determines power of cross-tissue modules from reshaped gene expression matrix (tissue_gene x sample).
# 
# Input R object must contain reshaped gene expression matrix [mat] containing multiple
# tissues in rows and aligned column-wise by individual. Must also contain associated row_meta table [row_meta].
# 
# Makes folders: output in working directory containing power estimate statistics.
#
# example of usage:
# determinePower.R ~/DataProjects/cross-tissue/STARNET/gene_exp_norm_batch_imp/all.RData

# Set options
# -----------------------------------------------------------

opts = list()
# Picking a soft power threshold to get scale-free correlation networks from each tissue alone
# opts$powers = seq(1, 10, length.out=20)
# opts$powers = seq(0.1, 10, length.out=51)
opts$powers = seq(0.5, 10, length.out=20)
# opts$powers = seq(1, 5, length.out=50)
# opts$powers = seq(0.1, 10, length.out=50)
# block_size = 4000

opts$n_breaks = 50  # discrete bins of connectivity distribution

# opts$abs_cor_min = 0.2  # minimum correlation to be considered
# cor_quant = 0.95


# Parse and check user input
# -------------------------------------------------------------
args = commandArgs(trailingOnly=TRUE)  # R-style positional arguments
if (length(args) != 1) {
	stop("Invalid arguments. USAGE: determinePower.R <path to cross-tissue expression object>")
}
opts$emat_file = args[1]

# For running in interactive mode
if (!exists("args")) {
	setwd("/Users/sk/Google Drive/projects/cross-tissue/co-expression/determinePower")
	# setwd("/sc/orga/projects/STARNET/koples01/cross-tissue/co-expression/determinePower")

	data_dir = "/Users/sk/DataProjects/cross-tissue"

	# Load imputed recast gene expression matrix
	opts$emat_file = file.path(data_dir, "STARNET/gene_exp_norm_batch_imp/all.RData")
}


# Load reshaped gene expression matrix
load(opts$emat_file, verbose=TRUE)

# Convert table from reshape2
if (exists("expr_recast")) {
	mat = expr_recast[, 3:ncol(expr_recast)]
	mat = data.matrix(mat)
	row_meta = expr_recast[, 1:2]
	row_meta = as.data.frame(row_meta)
}

# Check input
if (!exists("mat")) {
	stop("R object does not contain mat variable.")
}

if (!exists("row_meta")) {
	stop("Input R object does not contain meta_row variable")
}

if (any(is.na("mat"))) {
	warning("Gene expression matrix contains missing data.")
}

# Set up R environment
# -------------------------------------------------------------
library(data.table)
# library(reshape2)
# library(plyr)
# library(parallel)
# library(RColorBrewer)
library(Matrix)

library(compiler)
enableJIT(3)

# Load and setup WGCNA
library(WGCNA)
options(stringsAsFactors = FALSE)

enableWGCNAThreads(nThreads=4)  # assuming 2 threads per core


# Evaluate sequence of soft network cutoffs within tissues
# ------------------------------------------------------------
con_eval = list()
for (i in 1:length(unique(row_meta$tissue))) {
	tissue = unique(row_meta$tissue)[i]
	message("Estimating power for: ", tissue)

	# Symmetric correlation matrix within-tissue
	message("estimating correlation coefficients")
	cmat = corFast(t(mat[row_meta$tissue == tissue, ]),
		nThreads=4  # parallel computation not working on OSX
	)
	gc()
	message("dim: ", nrow(cmat), ", ", ncol(cmat))

	# Convert to adjacency matrix
	message("making adjacency matrix")
	cmat = abs(cmat)
	gc()

	# opts$abs_cor_min = quantile(cmat, cor_quant)
	# gc()

	# # Drop entries below threshold
	# message("converting to sparse matrix, abs_cor_min: ", opts$abs_cor_min)
	# cmat[cmat < opts$abs_cor_min] = 0  # delete entries
	# gc()

	# # Convert to sparse, symmetric matrix, for increased efficiency
	# cmat = forceSymmetric(cmat)
	# gc()

	# cmat = Matrix(cmat, sparse=TRUE)
	# gc()

	# message("non-zero coefficients: ", nnzero(cmat))
	# gc()

	# Returns scale-free distribution fits and connectivity vectors for each tissue
	# at the power series.
	fits = sapply(opts$powers, function(pow) {
		message("power: ", pow)
		# message("\tadjusting weights...")
		cmat_power = cmat^pow
		gc()
		# message("\tcalculating connectivity")
		# Calculate connectivity for nodes in both tissues
		k = apply(cmat_power, 1, sum) - 1  # subtracting diagonal
		gc()

		# message("\tfitting scale-free connectivity model")
		fit = scaleFreeFitIndex(k, nBreaks=opts$n_breaks, removeFirst=TRUE)
		gc()
		return(fit)
		# return(list(fit=fit, k=k))
	})

	con_eval[[i]] = fits
	gc()
}


# Test optimal beta for combinations of tissues. Uses cross-correlations only.
# Between-tissues
# ----------------------------------------------------------
con_eval_pairs = list()
paired_tissue = combn(unique(row_meta$tissue), 2)
for (i in 1:ncol(paired_tissue)) {
	# Cross-correlation only
	message("calculating cross-tissue correlation: ", paired_tissue[1, i], ", ", paired_tissue[2, i])

	message("estimating correlation coefficients")
	cmat = corFast(
		t(mat[row_meta$tissue == paired_tissue[1, i], ]),
		t(mat[row_meta$tissue == paired_tissue[2, i], ]),
		nThreads=4
	)
	gc()  # garbage collection

	# Convert to adjacency matrix
	message("making adjacency matrix")
	cmat = abs(cmat)
	gc()

	# opts$abs_cor_min = quantile(cmat, cor_quant)
	# gc()

	# # Drop entries below threshold
	# message("converting to sparse matrix, abs_cor_min: ", opts$abs_cor_min)
	# cmat[cmat < opts$abs_cor_min] = 0
	# gc()

	# # Convert to sparse matrix, for increased efficiency
	# cmat = Matrix(cmat, sparse=TRUE)
	# message("non-zero coefficients: ", nnzero(cmat))
	# gc()

	# Returns scale-free distribution fits and connectivity vectors for each tissue
	# at the power series.
	fits = sapply(opts$powers, function(pow) {
		message("power: ", pow)
		# message("\tadjusting weights...")
		cmat_power = cmat^pow
		gc()
		# message("\tcalculating connectivity")
		# Calculate connectivity for nodes in both tissues
		k1 = apply(cmat_power, 1, sum)
		k2 = apply(cmat_power, 2, sum)

		# message("\tfitting scale-free connectivity model")
		fit = scaleFreeFitIndex(c(k1, k2), nBreaks=opts$n_breaks, removeFirst=TRUE)
		gc()
		return(fit)
		# return(list(fit=fit, k1=k1, k2=k2))
	})

	con_eval_pairs[[i]] = fits
	gc()
}
# names(con_eval_pairs) = apply(paired_tissue, 2, paste, collapse="_")

names(con_eval_pairs) = apply(paired_tissue, 2, function(fac) {
	paste(levels(paired_tissue)[fac], collapse="_")
})


dir.create("output")
save(opts, con_eval, con_eval_pairs, file="output/con_eval.RData")
# save(con_eval_pairs, file="output/thresh_eval_pairs.RData")
