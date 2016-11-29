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


# User input
beta = 2.41  # power determined for each individual tissue
data_dir = "/Users/sk/DataProjects/cross-tissue"  # root of data directory
setwd("/Users/sk/Google Drive/projects/cross-tissue")


# Overwrites blockwiseModules() function from WGCNA
# to avoid overflow errors for inputs expression matrices with many genes.
# line 738 (and other lines) changed to:
# if (sum(as.numeric(modGenes)) - sum(as.numeric(reassign)) < minModuleSize) 
source("lib/WGCNA/R/blockwiseModulesC.R")
environment(blockwiseModules) = asNamespace("WGCNA")  # attach function to namespace

setwd("co-expression")

# Loads expr_recast data frame with tissue-specific expression
load(file.path(data_dir, "STARNET/gene_exp_norm_reshape/expr_recast.RData"))

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
patient_ids = colnames(mat)

# Missing data
message("missing data fraction: ", sum(is.na(mat))/(nrow(mat) * ncol(mat)))

# Standardize expression data

# mat = t(mat)
mat = scale(t(mat))
mat[is.na(mat)] = 0.0  # impute to average

# Random sample of data, for running on subset of data
# feature_idx = sample(1:ncol(mat), 1000)

dir.create(file.path(data_dir, "TOMs"))  # directory for data
bwnet = blockwiseModules(
	mat,
	# mat[,feature_idx],
	power=beta,
	randomSeed=42000,
	maxBlockSize=20000,
	nThreads=2,
	saveTOMs=TRUE,
	saveTOMFileBase=file.path(data_dir, "TOMs/TOM"),
	verbose=2
)

# Store modules in data directory
dir.create(file.path(data_dir, "modules"))

# Store module data
save(bwnet, row_meta, patient_ids, file=file.path(data_dir, "modules", "cross-tissue.RData"))


# Load topological overlap matrix (TOM) into environment.
# a dist object.
tom1 = new.env()
load(file.path(data_dir, "TOMs/TOM-block.1.RData"), tom1)

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




# bwLabels = matchLabels(bwnet$colors, moduleLabels)
# plotDendroAndColors(bwnet$dendrograms[[1]], bwModuleColors[bwnet$blockGenes[[1]]],
# 	"Module colors", main = "Gene dendrogram and module colors in block 1",
# 	dendroLabels = FALSE, hang = 0.03,
# 	addGuide = TRUE, guideHang = 0.05)

unique_transcripts_clust = sapply(1:k, function(i) {
	length(unique(expr_recast$transcript_id[clust[[k]]$cluster == i]))
})
