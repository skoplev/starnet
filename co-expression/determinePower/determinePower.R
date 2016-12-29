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
load(file.path(data_dir, "STARNET/gene_exp_norm_batch_imp/all.RData"), verbose=TRUE)

if (any(is.na(mat))) {
	stop("Gene expression matrix contains missing data.")
}

# Standardize expression data
# mat_scaled = scale(t(mat))

# Match phenotype data to selected gene expression matrix
pheno_matched = pheno[match(colnames(mat), pheno$starnet.ID), ]

# Picking a soft power threshold to get scale-free correlation networks from each tissue alone
powers = seq(1, 10, length.out=20)
# powers = seq(1, 5, length.out=50)
# powers = seq(0.1, 10, length.out=50)
block_size = 4000

n_breaks = 50  # binning of connectivity

# abs_cor_min = 0.05  # minimum correlation to be considered
abs_cor_min = 0.1  # minimum correlation to be considered

# Evaluate sequence of soft network cutoffs
con_eval = list()
for (tissue in unique(row_meta$tissue)) {
	message("Estimating power for: ", tissue)
	# submat = mat[row_meta$tissue == tissue, ]

	# Symmetric correlation matrix within-tissue
	cmat = corFast(t(mat[row_meta$tissue == tissue, ]),
		nThreads=8  # parallel computation not working
	)

	# con_eval[[tissue]] = pickSoftThreshold.fromSimilarity(cmat,
	# 	# t(submat),
	# 	powerVector=powers,
	# 	verbose=5,
	# 	blockSize=block_size
	# )

	# Convert to adjacency matrix
	message("making adjacency matrix")
	cmat = abs(cmat)

	# Drop entries below threshold
	message("converting to sparse matrix")
	cmat[cmat < abs_cor_min] = 0.0

	# Convert to sparse matrix, for increased efficiency
	cmat = Matrix(cmat, sparse=TRUE)
	message("non-zero coefficients: ", sum(cmat > 0))
	gc()

	# Returns scale-free distribution fits and connectivity vectors for each tissue
	# at the power series.
	fits = lapply(powers, function(pow) {
		message("power: ", pow)
		message("\tadjusting weights...")
		cmat_power = cmat^pow
		message("\tcalculating connectivity")
		# Calculate connectivity for nodes in both tissues
		k = apply(cmat_power, 1, sum) - 1  # subtracting diagonal

		message("\tfitting scale-free connectivity model")
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
	# message("Estimating power for: ", paired_tissue[1, i], ", ", paired_tissue[2, i])
	# submat = mat[row_meta$tissue %in% paired_tissue[, i], ]

	# Cross-correlation only
	message("calculating cross-tissue correlation: ", paired_tissue[1, i], ", ", paired_tissue[2, i])
	cmat = corFast(
		t(mat[row_meta$tissue == paired_tissue[1, i], ]),
		t(mat[row_meta$tissue == paired_tissue[2, i], ]),
		nThreads=8  # might not work?
	)
	gc()  # garbage collection

	# con_eval_pairs[[i]] = pickSoftThreshold.fromSimilarity(cmat,
	# con_eval_pairs[[i]] = pickSoftThreshold(cmat,
	# 	dataIsExpr=FALSE,
	# 	# t(submat),
	# 	powerVector=powers,
	# 	verbose=2,
	# 	blockSize=block_size
	# )

	# Convert to adjacency matrix
	message("making adjacency matrix")
	cmat = abs(cmat)

	# Drop entries below threshold
	message("converting to sparse matrix")
	cmat[cmat < abs_cor_min] = 0.0

	# Convert to sparse matrix, for increased efficiency
	cmat = Matrix(cmat, sparse=TRUE)
	message("non-zero coefficients: ", sum(cmat > 0))
	gc()

	# Returns scale-free distribution fits and connectivity vectors for each tissue
	# at the power series.
	fits = lapply(powers, function(pow) {
		message("power: ", pow)
		message("\tadjusting weights...")
		cmat_power = cmat^pow
		message("\tcalculating connectivity")
		# Calculate connectivity for nodes in both tissues
		k1 = apply(cmat_power, 1, sum)
		k2 = apply(cmat_power, 2, sum)

		message("\tfitting scale-free connectivity model")
		fit = scaleFreeFitIndex(c(k1, k2), nBreaks=n_breaks, removeFirst=TRUE)
		gc()
		return(list(fit=fit, k1=k1, k2=k2))
	})

	con_eval_pairs[[i]] = fits
	gc()
}

names(thresh_eval_pairs) = apply(paired_tissue, 2, paste, collapse="_")
# thresh_eval_pairs = thresh_eval  # fix

dir.create("output")
save(thresh_eval, file="output/thresh_eval.RData")
save(thresh_eval_pairs, file="output/thresh_eval_pairs.RData")


i = 17

plot(thresh_eval_pairs[[i]]$fitIndices$Power,
	thresh_eval_pairs[[i]]$fitIndices$SFT.R.sq,
	main=names(thresh_eval_pairs)[i],
	type="l"
)


# i = 2
par(mfrow=c(3, 3))

for (i in 1:length(thresh_eval)) {

	tissue = names(thresh_eval)[i]
	plot(thresh_eval[[i]]$fitIndices$Power,
		thresh_eval[[i]]$fitIndices$SFT.R.sq,
		main=names(thresh_eval)[i],
		type="l",
		lwd=2.0,
		ylim=c(0, 1)
	)

	# Find matching

	pairs_idx = which(apply(paired_tissue, 2, function(col) tissue %in% col))
	for (idx in pairs_idx) {
		lines(thresh_eval_pairs[[idx]]$fitIndices$Power,
			thresh_eval_pairs[[idx]]$fitIndices$SFT.R.sq,
			main=names(thresh_eval_pairs)[idx],
			type="l", col="grey"
		)
	}

}


sapply(thresh_eval_pairs, function(x) x$powerEstimate)

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