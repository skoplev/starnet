# Analysis for basic distributions for the modules detected using 3 different
# methods.
#
rm(list=ls())

# Analysis of modules
library(reshape2)
library(RColorBrewer)
library(gplots)
library(magicaxis)

data_dir = "/Users/sk/DataProjects/cross-tissue"  # root of data directory
setwd("/Users/sk/Google Drive/projects/cross-tissue")

source("src/base.R")

# Load gwnet and row_meta tables
# load(file.path(data_dir, "modules/cross-tissue.RData"))

load(file.path(data_dir, "modules/between_within-cross-tissue.RData"), verbose=TRUE)

# tissue_col = brewer.pal(12, "Set3")
# tissue_col = brewer.pal(9, "Pastel1")
tissue_col = brewer.pal(9, "Set1")[-6]

modules = as.integer(factor(bwnet$colors))  # the modules detected


# Calculate module frequency statistics
# Gene block tissue composition
# block_tissue_comp = lapply(bwnet$blockGenes, function(idx) {
block_tissue_comp = lapply(unique(bwnet$gene_blocks), function(block_id) {
	idx = bwnet$gene_blocks == block_id
	# tissue2 = meta_genes$tissue[idx]
	# return(table(tissue2))
	return(table(meta_genes$tissue[idx]))
})

module_tissue_comp = lapply(1:length(unique(modules)), function(i) {
	# tissue = meta_genes$tissue[modules == i]
	# return(table(tissue))
	return(table(meta_genes$tissue[modules == i]))
})

block_tissue = countMat(block_tissue_comp)

module_tissue = countMat(module_tissue_comp)
# module_tissue = module_tissue[, apply(module_tissue, 2, sum) < 500]

module_size = apply(module_tissue, 2, sum)
# frac_secondary = 1 - apply(module_tissue, 2, max) / apply(module_tissue, 2, sum)


# Calculate tissue frequencies
module_tissue_freq = sweep(module_tissue, 2, apply(module_tissue, 2, sum), "/")

# Calculate tissue entropy
tissue_entropy = apply(module_tissue_freq, 2, entropy)

hc = hclust(dist(t(module_tissue_freq)))

min_entro = 0.1
exclude = module_size > 5000
sel_modules = tissue_entropy > min_entro & !exclude

sum(tissue_entropy > min_entro)

# Hierarcical cluster of tissue frequencies for all modules
# Used to organize barplots
hc_sel = hclust(dist(t(module_tissue_freq[,sel_modules])))


# svg("co-expression/plots/cross-tissue-overview2.svg")
pdf("co-expression/plots/cross-tissue-overview2.pdf")
par(mfcol=c(2, 2), lwd=1.0)
# plot(density(tissue_entropy), main="Module tissue entropy")  # not used
plot(density(tissue_entropy), xlim=c(0, 1), main="Module tissue entropy")  # not used

abline(v=min_entro, col="grey", lty=3)

# WCGNA block tissue composition
barplot(block_tissue, col=tissue_col,
	main="K-means partitioned WGCNA", ylab="Genes", xlab="Blocks")
# legend("topright", legend=rownames(block_tissue), pch=15, col=tissue_col)


# Plot of tissue entropy and size
mod_cols = rep("black", length(tissue_entropy))
mod_cols[sel_modules] = "red"
magplot(tissue_entropy, log10(module_size),
	col=mod_cols,
	main="Cross-tissue modules",
	xlab="Tissue entropy",
	ylab="Genes",
	unlog="y"
)

# cross_tissue_modules = module_tissue[, order(tissue_entropy, decreasing=TRUE)[2:k]]
cross_tissue_modules = module_tissue[, sel_modules]
cross_tissue_modules = cross_tissue_modules[, hc_sel$order]  # rerder based on hclust of tissue frequencies

# cross_tissue_modules = cross_tissue_modules[,
# 	order(apply(cross_tissue_modules, 2, sum), decreasing=TRUE)]

# par(lwd=0.4)
barplot(cross_tissue_modules,
	col=tissue_col,
	main=paste0("Cross-tissue modules, n=", sum(sel_modules)),
	ylab="Genes",
	xlab="By tissue composition",
	border=NA, space=0
)
legend("topright", legend=rownames(block_tissue),
	pch=22,
	pt.bg=tissue_col,
	cex=0.7
)
dev.off()

# Heatmap of eigengenes
x = as.matrix(bwnet$MEs)

# heatmap.2(
# 	x,
# 	trace="none",
# 	col=colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
# 	breaks=seq(-alpha, alpha, length.out=101)
# )

pdf("plots/cross-tissue-module-eigengenes.pdf", height=5)
alpha = 0.05  # eigengene color range, +-
x = as.matrix(bwnet$MEs)
x = x[, sel_modules]
heatmap.2(
	x,
	trace="none",
	col=colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
	breaks=seq(-alpha, alpha, length.out=101)  # cap of coloring 
)
dev.off()
