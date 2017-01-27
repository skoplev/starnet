# Detect modules using blockwise WGCNA and heterogenous beta values.

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
opts$data_dir = "/Users/sk/DataProjects/cross-tissue"  # root of data directory
opts$emat_file = "STARNET/gene_exp_norm_reshape/expr_recast.RData"
opts$beta_mat_file = "co-expression/determinePower/output/beta.csv"
opts$project_root = "/Users/sk/Google Drive/projects/cross-tissue"

opts$beta = 3.0  # agerage from cross-tissue

opts$beta_within = 5.2
opts$beta_between = 2.7


# methods: 
#	single, uses a single beta value
#	between_within, uses two different beta values for within and between tissue
#	complete, uses a complete specification of beta  values
opts$method = "single"
opts$method = "between_within"
opts$method = "complete"

opts$max_block_size = 20000

# Test case samples a subset of the expression matrix
opts$test = TRUE
# opts$test = FALSE


# Parse options
# ------------------------------------------------
setwd(opts$project_root)

# Load optimal tissue-specific and between tissue beta values
opts$beta_mat = read.table(, row.names=1, header=TRUE, sep=",")
opts$beta_mat = data.matrix(opts$beta_mat)

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
source("lib/WGCNA/R/blockwiseModulesC.R")
environment(blockwiseModules) = asNamespace("WGCNA")  # attach function to namespace
environment(projectiveKMeans) = asNamespace("WGCNA")  # attach function to namespace
environment(TOMsimilarity) = asNamespace("WGCNA")  # attach function to namespace

# Parse gene expression data
# --------------------------------------
# Loads expr_recast data frame with tissue-specific expression
load(file.path(opts$data_dir, opts$emat_file))

# Get expression matrix from recasted data frame with all tissue 
emat = expr_recast[, 3:ncol(expr_recast)]
meta_gene = expr_recast[, 1:2]
meta_gene = as.data.frame(meta_gene)

# rownames(emat) = expr_recast$transcript_id
colnames(emat) = colnames(expr_recast)[3:ncol(expr_recast)]
patient_ids = colnames(emat)

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
	emat = emat[, feature_idx]
	meta_gene = meta_gene[feature_idx, ]
}

dir.create(file.path(opts$data_dir, "TOMs"))  # directory for data

if (opts$method == "single") {
	# Blockwise modules with single beta value

	bwnet = blockwiseModules(emat,
		power=opts$beta,
		randomSeed=42000,
		maxBlockSize=opts$max_block_size,
		nThreads=2,
		saveTOMs=TRUE,
		saveTOMFileBase=file.path(opts$data_dir, "TOMs/TOM"),
		verbose=2
	)
else if (opts$method == "between_within" | opts$methods == "complete") {
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
		verbose=2, indent=1)

	gene_blocks = .orderLabelsBySize(clust$clusters)

	# Loop over each block, identifying modules using WGCNA
	modules = rep(NA, ncol(emat))
	for (block in unique(gene_blocks)) {
		block_gene_idx = gene_blocks == block
		block_tissues = meta_gene$tissue[block_gene_idx]

		# Calculate Pearson's co-expression correlations
		adj_mat = cor(emat[, block_gene_idx])
		adj_mat = abs(adj_mat)  # adjacency matrix

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
				adj_mat[idx1, idx2] = adj_mat[idx1, idx2] ^ opts$beta_mat[i, j ]
			}
		}

		# Topological overalp matrix
		tom_dist = 1 - TOMsimilarity(adj_mat)

		#
		hc_tree = hclust(as.dist(tom_dist), method="average")

		min_module_size = 30

		dynamic_modules = cutreeDynamic(hc_tree,
			distM=tom_dist,
			deepSplit=2, pamRespectsDendro=FALSE,
			minClusterSize=min_module_size)

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
	bwnet$colors=modules
} else {
	stop("Invalid method: " method)
}

# Store modules in data directory
dir.create(file.path(opts$data_dir, "modules"))

# Store module data
save(bwnet, meta_gene, patient_ids, opts, file=file.path(opts$data_dir, "modules", paste0(opts$method, "-cross-tissue.RData")))


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

