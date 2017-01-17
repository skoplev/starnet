library(RColorBrewer)
library(reshape2)
library(gplots)

setwd("/Users/sk/Google Drive/projects/cross-tissue/co-expression/determinePower")

# Calculates optimal beta values from matrix of linear regression fits
# using the scale-free R^2 criterion from WGCNA.
# FInds lowest beta above threhshold
getOptimalBeta = function(con_eval, powers, beta_min=0.85) {
	betas = lapply(con_eval, function(x) {
		# Exclude powers below 1.
		include_powers = powers > 0.999
		pow = powers[include_powers]
		beta = unlist(x[1, include_powers])

		# first beta above threshold
		best_pow = pow[min(which(beta > beta_min))]
		# or highest beta
		if (is.na(best_pow)) {
			# best_pow = pow[which.max(beta)]
			best_pow = pow[min(which(beta / max(beta, na.rm=TRUE) > 0.975), na.rm=TRUE)]
		}
		return(best_pow)
	})
	return(unlist(betas))
}


# Load outout
load('~/DataProjects/cross-tissue/STARNET/determine_power/con_eval2.RData')


# Load tissue specification. Not strictly necesary
load(opts$emat_file, verbose=TRUE)

data_dir = "/Users/sk/DataProjects/cross-tissue"

# Load imputed recast gene expression matrix
opts$emat_file = file.path(data_dir, "STARNET/gene_exp_norm_batch_imp/all.RData")

# Load reshaped gene expression matrix
load(opts$emat_file, verbose=TRUE)

# Get tissue names
tissues = unique(row_meta$tissue)

names(con_eval) = tissues

# R^2 powers plot for each tissue
i = 1


tissue_col = brewer.pal(9, "Set1")[c(1:5, 7:9)]

# Separate tissue codes
tissue_pairs = lapply(
	names(con_eval_pairs),
	function(s) strsplit(s, "_")[[1]]
)


pdf("plots/beta_series.pdf", height=9, width=5)
par(mfrow=c(4, 2))
for (i in 1:length(con_eval)) {
	tissue = names(con_eval)[i]

	# Plot within-tissue R^2 curve
	plot(opts$powers,
		con_eval[[i]][1,],
		type="l",
		main=tissue,
		ylim=c(0, 1),
		lwd=2,
		ylab=expression(R^2),
		xlab=expression(beta)
	)

	# find pairwise ids involving ith tissue
	idx = which(sapply(tissue_pairs, function(x) tissue %in% x))

	for (j in idx) {
		to_tis = tissue_pairs[[j]][tissue_pairs[[j]] != tissue]

		lines(opts$powers,
			con_eval_pairs[[j]][1,],
			col=tissue_col[match(to_tis, tissues)],
			lwd=1.5
		)
		tissue_pairs
	}
}

plot(0, 0, type="n")
legend("bottomleft", legend=tissues, col=tissue_col, pch=15, cex=1.0)
dev.off()

# Get matrix of optimal beta values

tissue_opt = data.frame(Reduce(rbind, tissue_pairs))
colnames(tissue_opt) = c("tissueA", "tissueB")

tissue_opt$beta = getOptimalBeta(con_eval_pairs, opts$powers)


# Dublicate for each tissue combination, for making the symmetrical parts of the beta matrix
tissue_inv = tissue_opt
tissue_inv[, c(1, 2)] = tissue_opt[, c(2, 1)]  # swap

tissue_opt = rbind(tissue_opt, tissue_inv)  # combine


# within tissue beta balues
tissue_within_opt = data.frame(
	tissueA=tissues, tissueB=tissues,
	beta=getOptimalBeta(con_eval, opts$powers))

tissue_opt = rbind(tissue_opt, tissue_within_opt)

# Reshape data into beta matrix
beta_df = dcast(tissue_opt, tissueA~tissueB)
beta_mat = data.matrix(beta_df[, 2:ncol(beta_df)])
rownames(beta_mat) = beta_df[, 1]

# Plot beta matrix as heatmap
pdf("plots/betas.pdf", width=5, height=5)
heatmap.2(beta_mat,
	cellnote=beta_mat,
	trace="none",
	col=colorRampPalette(rev(brewer.pal(9, "Spectral")))(21),
	notecol="black"
)
dev.off()

# --------------------------

plot(opts$powers,
	con_eval_pairs[[i]][1,],
	type="l",
	main=names(con_eval_pairs)[i]
)

plot(con_eval_pairs[[i]]$fitIndices$Power,
	con_eval_pairs[[i]]$fitIndices$SFT.R.sq,
	main=names(con_eval_pairs)[i],
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