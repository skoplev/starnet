#!/usr/bin/env Rscript

# Detect modules using blockwise WGCNA and heterogenous beta values.
#
# Input file must be an R object containing either:
#   - mat; matrix with data matrix, and row_meta; table with gene meta data.
#   - expr_recast; data frame with reshaped data containing metadata in 2 first columns followed by gene expression matrix.

rm(list=ls())

library(data.table)
library(reshape2)
library(plyr)
library(parallel)
library(RColorBrewer)

library(compiler)
enableJIT(3)

# Load and setup WGCNA
library(WGCNA)
options(stringsAsFactors = FALSE)

enableWGCNAThreads(nThreads=2)  # assuming 2 threads per core


# Config
# --------------------------
opts = list()

# Files and folders

# # Local
# opts$data_dir = "~/DataProjects/cross-tissue"  # root of data directory
# opts$emat_file = "STARNET/gene_exp_norm_reshape/expr_recast.RData"

# Minerva
opts$data_dir = "~/links/STARNET/koples01/data"  # root of data directory
# opts$emat_file = "cross-tissue/gene_exp_norm_batch_imp/all.RData"
opts$emat_file = "cross-tissue/gene_exp_norm_reshape/expr_recast.RData"


opts$out_folder = "modules"  # output folder in data directory, WARNING: overwrites modules in this folder!
opts$project_root = ".."  # relative path, for loading additional libraries
opts$beta_mat_file = "determinePower/output/beta.csv"


opts$beta = 3.0  # agerage from cross-tissue

opts$beta_within = 5.2
ots$beta_between = 2.7

# methods: 
#	single, uses a single beta value.
#	between_within, uses two different beta values for within and between tissue.
#	complete, uses a complete specification of beta values.
opts$method = "single"
# opts$method = "between_within"
# opts$method = "complete"

opts$max_block_size = 40000
opts$min_module_size = 30

# Test case samples a subset of the expression matrix
# opts$test = TRUE
opts$test = FALSE


# Parse options
# ------------------------------------------------
# setwd(opts$project_root)

# Load optimal tissue-specific and between tissue beta values
opts$beta_mat = read.table(opts$beta_mat_file, row.names=1, header=TRUE, sep=",")
opts$beta_mat = data.matrix(opts$beta_mat)

if (!all(rownames(opts$beta_mat) == colnames(opts$beta_mat))) {
	stop("Mismatch between rownames and column names of beta matrix.")
}

if (! nrow(opts$beta_mat) > 1) {
	stop("Invalid beta matrix dimension.")
}

if (opts$method == "between_within") {
	# Change the beta matrix
	diag(opts$beta_mat) = opts$beta_within
	opts$beta_mat[lower.tri(opts$beta_mat) | upper.tri(opts$beta_mat)] = opts$beta_between
}


# Load modified WGCNA code
# ------------------------------------------------

# Overwrites blockwiseModules() function from WGCNA
# to avoid overflow errors for inputs expression matrices with many genes.
# line 738 (and other lines) changed to:
# if (sum(as.numeric(modGenes)) - sum(as.numeric(reassign)) < minModuleSize) 
source(file.path(opts$project_root, "lib/WGCNA/R/blockwiseModulesC.R"))
environment(blockwiseModules) = asNamespace("WGCNA")  # attach function to namespace
environment(projectiveKMeans) = asNamespace("WGCNA")  # attach function to namespace
environment(TOMsimilarity) = asNamespace("WGCNA")  # attach function to namespace

# Parse gene expression data
# --------------------------------------
# Loads expr_recast data frame with tissue-specific expression
load(file.path(opts$data_dir, opts$emat_file))

if (exists("expr_recast")) {
	# Table with both expression and meta data in first two columns
	# Get expression matrix from recasted data frame with all tissue 
	emat = expr_recast[, 3:ncol(expr_recast)]
	meta_genes = expr_recast[, 1:2]
	meta_genes = as.data.frame(meta_genes)

	meta_genes$gene_symbol = sapply(strsplit(as.character(meta_genes$transcript_id), "_"), function(x) x[1])
	meta_genes$ensembl = sapply(strsplit(as.character(meta_genes$transcript_id), "_"), function(x) x[2])

	# rownames(emat) = expr_recast$transcript_id
	colnames(emat) = colnames(expr_recast)[3:ncol(expr_recast)]
	patient_ids = colnames(emat)
} else if (exists("mat") & exists("row_meta")) {
	# Separate sata matrix and meta
	emat = mat
	meta_genes = row_meta
	patient_ids = colnames(emat)
} else {
	stop("Invalid data format in provided R object")
}

# Missing data
message("missing data fraction: ", sum(is.na(emat))/(nrow(emat) * ncol(emat)))

# Standardize expression data
# emat = t(emat)
emat = scale(t(emat))
emat[is.na(emat)] = 0.0  # impute to average

# Random sample of data, for running on subset of data.
# For test purposes only.
if (opts$test) {
	feature_idx = sample(1:ncol(emat), 1000)
	# feature_idx = sample(1:ncol(emat), 20000)
	emat = emat[, feature_idx]
	meta_genes = meta_genes[feature_idx, ]
}


if (opts$method == "single") {
	# Blockwise modules with single beta value
	dir.create(file.path(opts$data_dir, "TOMs"))  # directory for TOM data

	bwnet = blockwiseModules(emat,
		power=opts$beta,
		randomSeed=42000,
		maxBlockSize=opts$max_block_size,
		nThreads=2,
		saveTOMs=FALSE,
		saveTOMFileBase=file.path(opts$data_dir, "TOMs/TOM"),
		verbose=2)
} else if (opts$method == "between_within" | opts$method == "complete") {
	# Blockwise modules with different beta values

	# Preclustering centers
	n_preclust_centers = as.integer(min(
		ncol(emat) / 20,
		100*ncol(emat) / opts$max_block_size))

	# K-means clustering on the cross-tissue genes determinging blocks
	clust = projectiveKMeans(emat,
		preferredSize=opts$max_block_size, 
		checkData=FALSE,
		sizePenaltyPower=5,
		nCenters=n_preclust_centers, 
		verbose=2,
		indent=1)

	gene_blocks = .orderLabelsBySize(clust$clusters)

	# Loop over each block, identifying modules using WGCNA
	modules = rep(NA, ncol(emat))
	for (block in unique(gene_blocks)) {
		block_gene_idx = gene_blocks == block
		block_tissues = meta_genes$tissue[block_gene_idx]

		# Calculate Pearson's co-expression correlations
		netw_mat = cor(emat[, block_gene_idx])
		netw_mat = abs(netw_mat)  # adjacency matrix

		# Adjust weights based on specificed within- and cross-tissue beta values.
		# Loop over all tissue combinations
		for (i in 1:nrow(opts$beta_mat)) {
			for (j in 1:ncol(opts$beta_mat)) {
				tissue_i = rownames(opts$beta_mat)[i]
				tissue_j = colnames(opts$beta_mat)[j]

				# Find coordinates of the tissue combination
				idx1 = block_tissues == tissue_i
				idx2 = block_tissues == tissue_j

				if (all(!idx1) | all(!idx2)) next
				# otherwise, adjust
				netw_mat[idx1, idx2] = netw_mat[idx1, idx2] ^ opts$beta_mat[i, j ]
			}
		}

		# Topological overalap matrix
		# Overwrites for increased memory efficiency
		netw_mat = 1 - TOMsimilarity(netw_mat)

		# Hierarcical clustering
		hc_tree = hclust(as.dist(netw_mat), method="average")

		# Cut tree determining the modules
		dynamic_modules = cutreeDynamic(hc_tree,
			distM=netw_mat,
			deepSplit=2,
			pamRespectsDendro=FALSE,
			minClusterSize=opts$min_module_size)

		dynamic_colors = labels2colors(dynamic_modules)

		if (length(unique(dynamic_colors)) == 1) {
			# Block is single module, store as module ID
			warning("Block: ", block, " is considered as single module.")
			module_id = paste0(dynamic_colors, "_", block)
		} else {
			# Merge modules with similar expression profiles
			ME_list = moduleEigengenes(emat[, block_gene_idx],
				colors=dynamic_colors)

			ME_dist = 1 - cor(ME_list$eigengenes)

			ME_tree = hclust(as.dist(ME_dist), method="average")

			ME_merged = mergeCloseModules(emat[, block_gene_idx],
				dynamic_colors,
				cutHeight=0.2)

			# Block module IDs
			module_id = paste0(ME_merged$colors, "_", block)
		}

		# Store block module IDs
		modules[block_gene_idx] = module_id

		gc()
	}  # end of block computation

	# Calculate eigengenes for all modules, across blocks
	bwnet = moduleEigengenes(emat, colors=modules)
	bwnet$colors = modules
	bwnet$gene_blocks = gene_blocks
} else {
	stop("Invalid method: ", opts$method)
}

# Store modules in data directory
dir.create(file.path(opts$data_dir, opts$out_folder))

# Store module data
save(bwnet, meta_genes, patient_ids, opts, file=file.path(opts$data_dir, opts$out_folder, paste0(opts$method, "-cross-tissue.RData")))


# # Load topological overlap matrix (TOM) into environment.
# # a dist object.
# tom1 = new.env()
# load(file.path(data_dir, "TOMs/TOM-block.1.RData"), tom1)

# # library(igraph)
# library(emdbook)  # for lseq
# library(magicaxis)  # for magplot

# # adjacency matrix
# adj_mat = as.matrix(tom1$TOM)

# # Calculate degree
# deg = apply(adj_mat, 1, sum)
# hist(deg, breaks=50)

# # Empirical CDF function
# deg_cdf = ecdf(deg)

# # log-log plots of degree distributions
# k = lseq(10, max(deg), length.out=20)
# par(mfrow=c(1, 2))
# scaleFreePlot(deg)
# magplot(log10(k), log10(1 - deg_cdf(k)),
# 	pch=16,
# 	xlab=expression("log"[10] * "k"),
# 	ylab=expression("log"[10] * " p(k)"),
# 	unlog="xy")

