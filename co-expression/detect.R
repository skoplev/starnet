# Detect modules
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

setwd("/Users/sk/Google Drive/projects/cross-tissue/co-expression")

# source("../src/parse.r")

# Loads expr_recast data frame with tissue-specific expression
load("/Users/sk/DataProjects/cross-tissue/STARNET/gene_exp_norm_reshape/expr_recast.RData")

beta = 2.41  # power determined for each individual tissue

# # Count missing values 
# missing_values = apply(expr_recast, 1, function(row) sum(is.na(row)))

# # Subsample expression matrix
# k = 10000
# sub_mat = expr_recast[sample(nrow(expr_recast), k), 3:ncol(expr_recast)]
# sub_mat = as.matrix(sub_mat)

# Get expression matrix from recasted data frame with all tissue 
mat = expr_recast[, 3:ncol(expr_recast)]
row_meta = expr_recast[, 1:2]
# rownames(mat) = expr_recast$transcript_id
colnames(mat) = colnames(expr_recast)[3:ncol(expr_recast)]

# Missing data
message("missing data fraction: ", sum(is.na(mat))/(nrow(mat) * ncol(mat)))

# Standardize expression data
mat = scale(t(mat))
mat[is.na(mat)] = 0.0  # impute to average


# Random sample of data, for running on subset of data
# feature_idx = sample(1:ncol(mat), 20000)

dir.create("TOMs")  # directory for data
bwnet = blockwiseModules(
	mat,
	# mat[,feature_idx],
	power=beta,
	maxBlockSize=20000,
	nThreads=2,
	saveTOMs=TRUE,
	saveTOMFileBase="TOMs/TOM",
	verbose=5
)

# Load TOM into environment
tom1 = new.env()
load("TOMs/TOM-block.2.RData", tom1)

# library(igraph)
library(emdbook)  # for lseq
library(magicaxis)  # for magplot

# adjacency matrix
adj_mat = as.matrix(tom1$TOM)

# Calculate degree
deg = apply(adj_mat, 1, sum)
hist(deg, breaks=50)

# Empirical CDF function
deg_cdf = ecdf(deg)

# log-log plots of degree distributions
k = lseq(10, max(deg), length.out=20)
par(mfrow=c(1, 2))
scaleFreePlot(deg)
magplot(log10(k), log10(1 - deg_cdf(k)),
	pch=16,
	xlab=expression("log"[10] * "k"),
	ylab=expression("log"[10] * " p(k)"),
	unlog="xy")


# bwnet$blockGenes
# Gene block tissue composition
block_tissue_comp = lapply(bwnet$blockGenes, function(idx) {
	idx = feature_idx[idx]  # translate to subset
	tissue = row_meta$tissue[idx]
	return(table(tissue))
})

modules = as.integer(factor(bwnet$colors))

table(modules)

# i = 24
module_tissue_comp = lapply(1:length(unique(modules)), function(i) {
	# row_meta[feature_idx[modules == i], ]
	tissue = row_meta$tissue[feature_idx[modules == i]]
	return(table(tissue))
})

# Heatmap of eigen
library(gplots)
x = as.matrix(bwnet$MEs)
heatmap.2(
	x,
	trace="none",
	col=colorRampPalette(brewer.pal(9, "RdBu"))(100)
)

# bwLabels = matchLabels(bwnet$colors, moduleLabels)
# plotDendroAndColors(bwnet$dendrograms[[1]], bwModuleColors[bwnet$blockGenes[[1]]],
# 	"Module colors", main = "Gene dendrogram and module colors in block 1",
# 	dendroLabels = FALSE, hang = 0.03,
# 	addGuide = TRUE, guideHang = 0.05)

unique_transcripts_clust = sapply(1:k, function(i) {
	length(unique(expr_recast$transcript_id[clust[[k]]$cluster == i]))
})
