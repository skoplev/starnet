#!/usr/bin/env Rscript
# Determines power of modules

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

enableWGCNAThreads(nThreads=8)  # assuming 2 threads per core

setwd("/Users/sk/Google Drive/projects/cross-tissue")
source("src/base.R")

data_dir = "/Users/sk/DataProjects/cross-tissue"

setwd("/Users/sk/Google Drive/projects/cross-tissue/co-expression/determinePower")
# setwd("/sc/orga/projects/STARNET/koples01/cross-tissue/co-expression/determinePower")


# STARNET phenotype data
pheno = fread(file.path(
	"/Volumes/SANDY/phenotype_data",
	"STARNET_main_phenotype_table.cases.Feb_29_2016.tbl"
))

# Load imputed recast gene expression matrix
load(file.path(data_dir, "STARNET/gene_exp_norm_batch_imp/all.RData"))



# Standardize expression data
# mat_scaled = scale(t(mat))
# mat_scaled = t(mat)  # no scale


# row_sd = apply(mat, 1, sd)

# sum(row_sd > 1)
# mat_scaled = t(mat[row_sd > 1,])

# mat[is.na(mat)] = 0.0  # impute to average

# Combined tSNE plot for all tissues

# Match phenotype data to selected gene expression matrix

pheno_matched = pheno[match(colnames(mat), pheno$starnet.ID), ]

library(tsne)
library(amap)

tsnePlot = function(embed, ...) {
	plot(embed[,1], embed[,2],
		xlab="", ylab="",
		xaxt="n", yaxt="n",
		...
	)
}


# emat is genes x samples
tsneSyntaxCmp = function(emat, pheno_matched, row_meta) {

	# Calculate correlations with SYNTAX score
	syntax_cor = cor(t(emat), pheno_matched$syntax_score, use="pairwise.complete")

	# Select highly correlated transcripts
	sel_transcripts = which(abs(syntax_cor) > 0.1)

	message("Selecting: ", length(sel_transcripts), " transcripts.")
	print(table(row_meta$tissue[sel_transcripts]))

	# Parallelized calculation of distance matrix
	dmat = Dist(t(emat), method="euclidean", nbproc=6)

	# Run tSNE
	embed = tsne(dmat, max_iter=2000, perplexity=20)

	# Plot results
	par(mfrow=c(2, 2), mar=c(2, 2, 2, 3))
	col = colorGradient(pheno_matched$syntax_score,
		gradlim=c(0, 100),
		colors=rev(brewer.pal(9, "Spectral")))
	tsnePlot(embed,
		col=col,
		pch=16,
		main="SYNTAX"
	)
	legendCol(colorRampPalette(rev(brewer.pal(9, "Spectral")))(20), c(0, 100))
	# text(embed[,1], embed[,2], labels=pheno_matched$syntax_score, cex=0.3)

	col = colorGradient(pheno_matched$DUKE,
		gradlim=c(0, 100),
		colors=rev(brewer.pal(9, "Spectral")))
	tsnePlot(embed,
		col=col,
		pch=16,
		main="DUKE"
	)
	# legendCol(rev(brewer.pal(9, "Spectral")), c(0, 100))
	legendCol(colorRampPalette(rev(brewer.pal(9, "Spectral")))(20), c(0, 100))
	# text(embed[,1], embed[,2], labels=pheno_matched$DUKE, cex=0.3)

	col = colorGradient(pheno_matched$ndv,
		# gradlim=range(pheno$ndv, na.rm=TRUE),
		gradlim=c(0, 3),
		colors=rev(brewer.pal(9, "Spectral")))
	tsnePlot(embed,
		col=col,
		pch=16,
		main="ndv"
	)
	legendCol(colorRampPalette(rev(brewer.pal(9, "Spectral")))(20), c(0, 3))
	# text(embed[,1], embed[,2], labels=pheno_matched$ndv, cex=0.3)

	col = colorGradient(pheno_matched$lesions,
		gradlim=c(0, 9),
		colors=rev(brewer.pal(9, "Spectral")))
	tsnePlot(embed,
		col=col,
		pch=16,
		main="lesions"
	)
	legendCol(colorRampPalette(rev(brewer.pal(9, "Spectral")))(20), c(0, 9))

	return(embed)
}

# Combined distance matrix
# dmat = dist(t(mat))
# dmat = Dist(mat_scaled[, 1:1000], method="euclidean", nbproc=4)
# dmat = Dist(mat_scaled, method="euclidean", nbproc=6)
# dmat = Dist(mat_scaled[, row_meta$tissue != "BLOOD"], method="euclidean", nbproc=6)
# dmat = Dist(mat_scaled, method="euclidean", nbproc=6)

# Calculate correlations
syntax_cor = cor(mat_scaled, pheno_matched$syntax_score, use="pairwise.complete")
duke_cor = cor(mat_scaled, pheno_matched$DUKE, use="pairwise.complete")
# starnet_cor = cor(mat_scaled, as.numeric(pheno_matched$starnet.ID), method="spearman", use="pairwise.complete")

sel_transcripts = which(abs(syntax_cor) > 0.1)
# sel_transcripts = which(abs(duke_cor) > 0.1)
# sel_transcripts = which(abs(duke_cor) > 0.15)


pdf("plots/syntax_cor_all_tSNE.pdf", width=7, height=6)
embed = tsneSyntaxCmp(mat, pheno_matched, row_meta)
dev.off()

pdf("plots/syntax_cor_all_tSNE_null.pdf", width=7, height=6)
embed = tsneSyntaxCmp(mat[, sample(ncol(mat))], pheno_matched, row_meta)
dev.off()


# table(row_meta$tissue[abs(starnet_cor) < 0.05])
# table(row_meta$tissue[abs(syntax_cor) > 0.1])


# col = colorGradient(pheno_matched$syntax_score,
# 	gradlim=c(0, 100),
# 	colors=rev(brewer.pal(9, "Spectral")))

# col = colorGradient(pheno_matched$syntax_score,
# 	gradlim=c(0, 100),
# 	colors=brewer.pal(9, "Blues"))

# col = colorGradient(pheno_matched$DUKE,
# 	gradlim=c(0, 100),
# 	colors=rev(brewer.pal(9, "Spectral")))

# col = colorGradient(pheno_matched$ndv,
# 	gradlim=c(0, 3),
# 	colors=brewer.pal(9, "YlGnBu"))

# col = colorGradient(pheno_matched$lesions,
# 	gradlim=range(pheno$lesions, na.rm=TRUE),
# 	colors=brewer.pal(9, "YlGnBu"))

col = colorGradient(pheno_matched$starnet.ID,
	gradlim=range(as.numeric(pheno_matched$starnet.ID), na.rm=T),
	# as.numeric(pheno_matched$starnet.ID),
	# as.numeric(covar_matched$subject),
	colors=rev(brewer.pal(9, "Spectral")))

col = colorGradient(factor(pheno_matched$Sex),
	gradlim=range(as.numeric(factor(pheno_matched$Sex)), na.rm=T),
	# as.numeric(pheno_matched$starnet.ID),
	# as.numeric(covar_matched$subject),
	colors=rev(brewer.pal(9, "Spectral")))


plot(embed[,1], embed[,2],
	col=col,
	pch=16,
	cex=1.5,
	xlab="", ylab=""
)

# text(embed[,1], embed[,2], labels=pheno_matched$syntax_score, cex=0.3)








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