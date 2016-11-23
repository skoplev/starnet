library(reshape2)
library(RColorBrewer)

# Load gwnet and row_meta tables
load(file.path(data_dir, "modules/cross-tissue.RData"))

# tissue_col = brewer.pal(12, "Set3")
tissue_col = brewer.pal(9, "Pastel1")
# tissue_col = brewer.pal(9, "Set1")


# List of counts to matrix of counts
# Note that tables must not be named, which can happen implicitly.
countMat = function(count_list) {
	mat = dcast(melt(count_list), L1~Var1, value.var="value")
	mat = mat[, 2:ncol(mat)]
	mat = as.matrix(mat)
	mat[is.na(mat)] = 0  # zero counts
	mat = t(mat)
	return(mat)
}


# Gene block tissue composition
block_tissue_comp = lapply(bwnet$blockGenes, function(idx) {
	# tissue2 = row_meta$tissue[idx]
	# return(table(tissue2))
	return(table(row_meta$tissue[idx]))
})

module_tissue_comp = lapply(1:length(unique(modules)), function(i) {
	# tissue = row_meta$tissue[modules == i]
	# return(table(tissue))
	return(table(row_meta$tissue[modules == i]))
})

block_tissue = countMat(block_tissue_comp)

module_tissue = countMat(module_tissue_comp)
# module_tissue = module_tissue[, apply(module_tissue, 2, sum) < 500]

module_size = apply(module_tissue, 2, sum)
frac_secondary = 1 - apply(module_tissue, 2, max) / apply(module_tissue, 2, sum)

entropy = function(x) {
	x = x[x != 0]  # exclude zero
	H = -sum(x * log10(x))
	return(H)
}


# Calculate tissue entropy
tissue_entropy = apply(module_tissue_freq, 2, entropy)


min_entro = 0.3
exclude = module_size > 5000
sel_modules = tissue_entropy > min_entro & !exclude

mod_cols = rep("black", length(tissue_entropy))
mod_cols[sel_modules] = "red"

par(mfrow=c(2, 2), lwd=1.0)

plot(density(tissue_entropy), main="Module tissue entropy")

# WCGNA block tissues
barplot(block_tissue, col=tissue_col, main="K-means WCGNA blocks", ylab="Genes")
# legend("topright", legend=rownames(block_tissue), pch=15, col=tissue_col)


magplot(tissue_entropy, log10(module_size),
	col=mod_cols,
	main="Cross-tissue modules",
	xlab="Tissue entropy",
	ylab="Genes",
	unlog="y"
)

# cross_tissue_modules = module_tissue[, order(tissue_entropy, decreasing=TRUE)[2:k]]
cross_tissue_modules = module_tissue[, sel_modules]
cross_tissue_modules = cross_tissue_modules[,
	order(apply(cross_tissue_modules, 2, sum), decreasing=TRUE)]
par(lwd=0.4)
barplot(cross_tissue_modules,
	col=tissue_col,
	ylab="Genes"
)
legend("topright", legend=rownames(block_tissue), pch=22, pt.bg=tissue_col)


# Heatmap of eigengenes
library(gplots)
x = as.matrix(bwnet$MEs)

alpha = 0.05  # eigengene color range
heatmap.2(
	x,
	trace="none",
	col=colorRampPalette(brewer.pal(9, "RdBu"))(100),
	breaks=seq(-alpha, alpha, length.out=101)
)


# Calculate tissue frequencies
# module_tissue_freq = sweep(module_tissue, 2, apply(module_tissue, 2, sum), "/")


plot(frac_secondary, log10(module_size))

barplot(module_tissue, col=tissue_col)


barplot(block_tissue_comp_mat)
	
barplot(block_tissue_comp)

# do.call(rbind, block_tissue_comp)
# merge(block_tissue_comp[[1]], block_tissue_comp[[2]])

modules = as.integer(factor(bwnet$colors))
# table(modules)
