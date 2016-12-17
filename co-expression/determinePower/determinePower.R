#!/usr/bin/env Rscript
# Determines power of modules

rm(list=ls())

library(data.table)
library(reshape2)
library(plyr)
library(parallel)
library(RColorBrewer)
library(impute)

library(compiler)
enableJIT(3)

# Load and setup WGCNA
library(WGCNA)
options(stringsAsFactors = FALSE)

enableWGCNAThreads(nThreads=8)  # assuming 2 threads per core

setwd("/Users/sk/Google Drive/projects/cross-tissue")
source("src/base.R")


setwd("/Users/sk/Google Drive/projects/cross-tissue/co-expression/determinePower")
# setwd("/sc/orga/projects/STARNET/koples01/cross-tissue/co-expression/determinePower")

# Loads expr_recast data frame with tissue-specific expression
load("/Users/sk/DataProjects/cross-tissue/STARNET/gene_exp_norm_reshape/expr_recast.RData")

# STARNET phenotype data
pheno = fread(file.path(
	"/Volumes/SANDY/phenotype_data",
	"STARNET_main_phenotype_table.cases.Feb_29_2016.tbl"
))


# Get expression matrix from recasted data frame with all tissue 
mat = expr_recast[, 3:ncol(expr_recast)]
mat = data.matrix(mat)
row_meta = expr_recast[, 1:2]
colnames(mat) = colnames(expr_recast)[3:ncol(expr_recast)]

# Remove samples with more than 50% missing data
missing_frac = apply(mat, 2, function(col) {
		sum(is.na(col)) / length(col)
})
mat = mat[, missing_frac < 0.5]

# Missing data message, remining samples
message("missing data fraction: ", sum(is.na(mat))/(nrow(mat) * ncol(mat)))

# Impute missing data using k-nearest neighbors
mat = impute.knn(mat, k=20, colmax=0.5)$data

# Standardize expression data
mat_scaled = scale(t(mat))

# mat[is.na(mat)] = 0.0  # impute to average

# Combined tSNE plot for all tissues

# Match phenotype data to selected gene expression matrix

pheno_matched = pheno[match(colnames(mat), pheno$starnet.ID), ]

# syntax_cor = cor(mat_scaled, pheno_matched$syntax_score, use="pairwise.complete")
# starnet_cor = cor(mat_scaled, as.numeric(pheno_matched$starnet.ID), method="spearman", use="pairwise.complete")

library(tsne)
library(amap)

# Combined distance matrix
# dmat = dist(t(mat))
# Parallelized calculation of distance matrix
# dmat = Dist(mat_scaled[, 1:1000], method="euclidean", nbproc=4)
# dmat = Dist(mat_scaled, method="euclidean", nbproc=6)
dmat = Dist(mat_scaled[, row_meta$tissue != "BLOOD"], method="euclidean", nbproc=6)

# dmat = Dist(mat_scaled[, abs(syntax_cor) > 0.1], method="euclidean", nbproc=6)
# dmat = Dist(mat_scaled[, abs(starnet_cor) < 0.05], method="euclidean", nbproc=6)

# table(row_meta$tissue[abs(starnet_cor) < 0.05])
# table(row_meta$tissue[abs(syntax_cor) > 0.1])

embed = tsne(dmat, max_iter=2000, perplexity=15)

col = colorGradient(pheno_matched$syntax_score,
	gradlim=c(0, 100),
	colors=rev(brewer.pal(9, "Spectral")))

col = colorGradient(pheno_matched$DUKE,
	gradlim=c(0, 100),
	colors=rev(brewer.pal(9, "Spectral")))

col = colorGradient(pheno_matched$ndv,
	gradlim=c(0, 3),
	colors=brewer.pal(9, "YlGnBu"))

col = colorGradient(pheno_matched$lesions,
	gradlim=range(pheno$lesions, na.rm=TRUE),
	colors=brewer.pal(9, "YlGnBu"))

col = colorGradient(pheno_matched$starnet.ID,
	gradlim=range(as.numeric(pheno_matched$starnet.ID), na.rm=T),
	# as.numeric(pheno_matched$starnet.ID),
	# as.numeric(covar_matched$subject),
	colors=rev(brewer.pal(9, "Spectral")))


plot(embed[,1], embed[,2],
	col=col,
	pch=16,
	cex=1.5,
	xlab="", ylab=""
)

plotColorBar(
	colorGradient(seq(0, 100, length.out=100),
		colors=brewer.pal(9, "YlGnBu")),
	min=0, max=100, title="SYNTAX")



# Picking a soft power threshold to get scale-free correlation networks from each tissue alone
powers = seq(1, 4, length.out=50)

# Evaluate sequence of soft network cutoffs
thresh_eval = list()
for (tissue in unique(row_meta$tissue)) {
	message("Estimating power for: ", tissue)
	submat = mat[, row_meta$tissue == tissue]

	thresh_eval[[tissue]] = pickSoftThreshold(submat,
		powerVector=powers,
		verbose=5,
		blockSize=4000
	)

	print(thresh_eval[[tissue]])
}


dir.create("output")
save(thresh_eval, file="output/thresh_eval.RData")


# Average best slope
deg_slopes = sapply(thresh_eval, function(x) x$fitIndices$slope)
mean_deg_slopes = apply(deg_slopes, 1, mean)

best_pow = which.min(abs(1 + mean_deg_slopes))
message("Best average beta for tissues: ", powers[best_pow])


dir.create("plots")

svg("plots/scale-free.svg", width=3.7, height=4)
colors = brewer.pal(9, "Set1")[c(1:5, 7:9)]
plot(powers, thresh_eval[[1]]$fitIndices$slope, type="n",
	main="Scale-free co-expression networks",
	ylim=c(-0.5, -1.5),
	xlab=expression(beta),
	ylab="log-log slope"
)

abline(h=-1.0, col="grey", lty=3, lwd=1.5)
abline(v=powers[best_pow], col="grey", lwd=1.5)
for (i in 1:length(thresh_eval)) {
	lines(powers, thresh_eval[[i]]$fitIndices$slope, col=colors[i], lwd=1.2)
	points(powers, thresh_eval[[i]]$fitIndices$slope, pch=16, col=colors[i], cex=0.4)
}
legend("bottomright", legend=names(thresh_eval), col=colors, pch=16, cex=0.5)
dev.off()

plot(powers, apply(deg_slopes, 1, mean))

# par(mfrow())
# lapply(thresh_eval, function(sft) {
# 	plot(powers, -sft$fitIndices$slope)
# 	abline(h=1)
# })

# # cex1=0.9
# plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
# 	xlab="Soft Threshold (power)",
# 	ylab="Scale Free Topology Model Fit,signed R^2",
# 	# type="n",
# 	main = paste("Scale independence"),
# 	pch=16
# )

# text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
# labels=powers,cex=cex1,col="red");