# Parses STARNET count matrices.
# Normalizes using DESeq2 size factors, writes matrices to
# Collapses by HUGO gene symbols using sum for multiple mappings
# to the same gene symbol.

rm(list=ls())

library(data.table)
library(org.Hs.eg.db)
library(DESeq2)
library(mgcv)
library(sva)
library(plyr)
library(MASS)
library(parallel)
library(penalized)
library(hexbin)


library(compiler)
enableJIT(3)

# library(parallel)


# Collapse rows of matrix by vector of row ids, by summing multiple rows with the same ids
# Returns matrix 
collapseMatSum = function(mat, row_ids) {
	mat = as.matrix(mat)

	unique_row_ids = unique(row_ids)
	unique_row_ids = unique_row_ids[!is.na(unique_row_ids)]

	# Index row numbers
	message("Indexing gene symbols")
	row_index = sapply(unique_row_ids, function(id) {
		return(which(id == row_ids))
	}, USE.NAMES=TRUE)

	message("Aggregating expression matrix by genes")
	mat_collapse = sapply(unique_row_ids, function(id) {
		submat = mat[row_index[[id]], , drop=FALSE]

		if (nrow(submat) == 1) {
			# No aggregation
			return(submat)
		} else {
			# Multiple rows, return column-wise sum
			return(
				colSums(submat, na.rm=TRUE)
			)
		}
		return(NA)
	})
	mat_collapse = t(mat_collapse)
	colnames(mat_collapse) = colnames(mat)

	return(mat_collapse)
}


# Constructs model matrix from comma-separated vector of multinomial variables.
# Each row is weighted to add to 1.
multiCatModelMatrix = function(catvec)  {
	catvec_sep = strsplit(catvec, ",")

	catvec_uni = unique(unlist(catvec_sep))
	catvec_uni = catvec_uni[!is.na(catvec_uni)]

	model_mat = matrix(0, ncol=length(catvec_uni), nrow=length(catvec))
	colnames(model_mat) = catvec_uni
	for (i in 1:length(catvec_sep)) {
		model_mat[i, ] = catvec_uni %in% catvec_sep[[i]]
	}

	multiplicity = apply(model_mat, 1, sum)
	multiplicity[multiplicity == 0] = 1  # avoid division by zero
	model_mat = sweep(model_mat, 1, multiplicity, "/")

	return(data.frame(model_mat))
}

# Replaces missing entries in each column with median or most frequenct value.
# 
imputeMedianModel = function(model) {
	model = as.data.frame(model)
	for (i in 1:ncol(model)) {
		# print(i)
		# count missing values
		n_missing = sum(is.na(model[, i]))
		if (n_missing > 0 & n_missing != nrow(model)) {
			# message(colnames(model)[i], " imputed for ", n_missing, " samples.")
			if (is.numeric(model[, i])) {
				# Impute to median
				model[is.na(model[, i]), i] = median(model[, i], na.rm=TRUE)
			} else {
				# Calculate frequencies
				freq = sort(table(model[, i]), decreasing=TRUE)
				# Impute to most frequent value
				model[is.na(model[, i]), i] = names(freq)[1]
			}
		}
	}
	if (sum(is.na(model)) > 0) {
		warning("imputeMedianModel() failed")
	}
	return(model)
}

# Finds full rank set of matrix columns.
# Recursive function
fullRankSubmat = function(mat) {
	if (qr(mat)$rank == ncol(mat)) {
		# Base case: full rank
		return(mat)
	} else {
		# Calc matrix ranks when removing each column
		rank_holdout = sapply(1:ncol(mat),
			function(i) qr(mat[, -i])$rank)

		remove_col = which.max(rank_holdout)
		message("Removing matrix column: ", colnames(mat)[remove_col])
		return(fullRankSubmat(mat[, -remove_col]))
	}
}


data_dir = "/Users/sk/DataProjects/cross-tissue"

setwd("/Users/sk/GoogleDrive/projects/STARNET/cross-tissue")



# Directory containing expression count matrices
emat_dir = file.path(data_dir, "STARNET/gene_exp/matrices")

# Load phenotype and covariate data for batch correction
pheno = fread(file.path(
	"/Volumes/SANDY/phenotype_data",
	"STARNET_main_phenotype_table.cases.Feb_29_2016.tbl"
))
pheno$Smoking.Years[pheno$Smoking.Years == "" | is.na(pheno$Smoking.Years)] = 0  # assumes that missing infor => no smoking
pheno$Smoking.Years[pheno$Smoking.Years == "25-30"] = "30"  # worst case assumption
pheno$Smoking.Years = as.numeric(pheno$Smoking.Years)

# Impute missing phenotype data. For sva analysis only.
pheno_imp = mice(pheno)

pheno = complete(pheno_imp)
# pheno = imputeMedianModel(pheno)

# First covariate table
covar = fread(file.path(
	"/Volumes/SANDY/phenotype_data",
	"covariates.cases_and_controls.April_12_2016.txt"
))
covar = rename(covar, c("sample"="id"))

# Second covariate table
covar2 = fread(file.path(
	"/Volumes/SANDY/phenotype_data",
	"covariates.tbl"
))

# Merge covariate tables
covar_merged = merge(covar, covar2,
	by=c("id", "sex", "age"),
	all=TRUE)
# Fix discrepant read_lenghts. Default read_length to first covar table.
covar_merged = rename(covar_merged, c("read_length.x"="read_length"))
# Fill missing read_length fields from covar2
covar_merged$read_length[is.na(covar_merged$read_length)] = covar_merged$read_length.y[is.na(covar_merged$read_length)]
# Drop read_length.y
covar_merged = covar_merged[, names(covar_merged) != "read_length.y", with=FALSE]
# Drop incomplete fields
covar_merged = covar_merged[, !names(covar_merged) %in% c("subject", "batch"), with=FALSE]



# Load count matrices
# gene symbol, exp
expr_files = list.files(emat_dir, "*.mat")

# Load all data matrices
expr_mats = lapply(expr_files, function(file_name) {
	d = fread(
		file.path(emat_dir, file_name))
	return(d)
})
names(expr_mats) = expr_files

lapply(expr_mats, dim)





# Normalize using size factors from DESeq2
# ---------------------------------------------------------------
# converts $id to rownames of returned matrices
message("Correcting for size factors.")
expr_mats_norm = sapply(expr_mats, function(emat) {
	# get numerical matrix
	mat = emat[,2:ncol(emat), with=FALSE]
	rownames(mat) = emat$id

	size_factors = estimateSizeFactorsForMatrix(as.matrix(mat))

	# divide each column by size factor
	norm_mat = sweep(mat, 2, size_factors, "/")

	return(norm_mat)
})

# # Write normalized matrices as .tsv files
# for (i in 1:length(expr_mats_norm)) {
# 	message("writing", names(expr_mats_norm)[i])

# 	write.table(expr_mats_norm[[i]],
# 		file.path(data_dir, "STARNET/gene_exp_norm", names(expr_mats_collapsed)[i]),
# 		sep="\t",
# 		quote=FALSE, col.names=NA
# 	)
# }
# save(expr_mats_norm, file=file.path(data_dir, "STARNET/gene_exp_norm/all.RData"))
# load(file.path(data_dir, "STARNET/gene_exp_norm/all.RData"))

# mat = expr_mats_norm[[1]]

# lapply(expr_mats_norm, dim)
frac_detect = 0.1  # fraction of samples that must detect 
message("Filtering transcripts based on detection limit")
expr_mats_norm = lapply(expr_mats_norm, function(mat) {

	# Calculate the number of samples containing some RNA
	rna_detected = apply(mat > 5, 1, sum)
	# Number of samples for each RNA
	n = apply(mat, 1, function(row) sum(!is.na(row)))

	# Filter transcript based on fraction of sample where the RNA is above detection limit
	rna_include = (rna_detected / n) > frac_detect

	mat = mat[rna_include,]

	return(mat)
})


# min_sd = 0.5  # at the count level
# # Variance filter and pseudo-log2 transform
# message("Filtering transcripts based on standard deviation")
# expr_mats_norm = lapply(expr_mats_norm, function(mat) { 
# 	# Filter genes based on standard deviation
# 	mat = mat[apply(mat, 1, sd, na.rm=T) > min_sd, ]

# 	return(mat)
# })


# Log transform
expr_mats_norm = lapply(expr_mats_norm,
	function(mat) log2(mat + 1))

# Rename expr mats
names(expr_mats_norm) = sapply(
	strsplit(names(expr_mats_norm), "[.]"),
	function(x) x[4]
)

# Exclude matrices
expr_mats_norm = expr_mats_norm[names(expr_mats_norm) != "COR"]

# Linear regression normalization.
# Uses SVD FLow Cell correction, covariates, surrogate variables (sva), and
# penalized regression to account for covariates.
# Phenotype data is imputed using MICE.
# ----------------------------------------------------------------

# Parameters:
# include flow cells with more than this number of samples
flow_cells_ignore_nsamples = 5  

# Remove flow cell eigenvector based on eigenvalue cutoff
flow_cells_min_eigen = 4.0

# Number of surrogate variables
nsurrogate_vars = 4
# nsurrogate_vars = 8

# Adjustment penality for regression model
adjust_l2_penalty = 1.0
# adjust_l2_penalty = 2.0

expr_mats_batch = lapply(seq_along(expr_mats_norm), function(i) {
	mat = expr_mats_norm[[i]]
	mat = as.matrix(mat)
	tissue = names(expr_mats_norm)[i]

	message("Normalizing ", tissue)

	# Match phenotype data
	patient_ids = sapply(
		strsplit(colnames(mat), "_"),
		function(x) x[2]
	)

	pheno_matched = pheno[match(patient_ids, pheno$starnet.ID), ]
	pheno_matched = as.data.frame(pheno_matched)

	# Match covariates
	covar_matched = covar_merged[match(colnames(mat), covar_merged$id), ]

	# Phenotype model
	# -----------------------------------------

	# Exclude phenotypes
	pheno_model = pheno_matched[, !colnames(pheno_matched) %in% c("starnet.ID", "Sex", "Age")]
	# pheno_model = pheno_model[, colnames(pheno_model) %in% c("BMI", "LDL", "syntax_score", "DUKE", "lesions", "ndv")]

	pheno_model = imputeMedianModel(pheno_model)

	# Remove singular covariates
	n_factors = apply(pheno_model, 2, function(x) length(levels(factor(x))))
	if (any(n_factors < 2)) {
		warning("Excluding singular covariates: ", colnames(pheno_model)[n_factors < 2])
	}
	pheno_model = pheno_model[, n_factors > 1]
	# pheno_model[is.na(pheno_model)] = 0


	# Covariate model
	# ------------------------------------------

	covar_model = covar_matched[, c("sex", "age", "read_length", "lab", "protocol"), with=FALSE]

	# Impute missing values from all model components to median value
	covar_model = imputeMedianModel(covar_model)

	covar_model$sex = factor(covar_model$sex)
	covar_model$age = as.numeric(covar_model$age)
	covar_model$read_length = factor(covar_model$read_length)
	covar_model$lab = factor(covar_model$lab)
	covar_model$protocol = factor(covar_model$protocol)

	# Remove singular covariates
	n_factors = apply(covar_model, 2, function(x) length(levels(factor(x))))
	if (any(n_factors < 2)) {
		warning("Excluding singular covariates: ", colnames(covar_model)[n_factors < 2])
	}
	covar_model = covar_model[, n_factors > 1]

	# covar_model[is.na(covar_model)] = 0

	# Flow Cell model
	# --------------------------------------------------------------
	# covar_merged_matched = covar_merged
	flow_model = multiCatModelMatrix(covar_matched$flowcell)

	# Remove flows with few samples
	flow_model = flow_model[, apply(flow_model, 2, sum) > flow_cells_ignore_nsamples]

	# Remove linear dependencies from flow_model using singular value decomposition
	flow_svd = svd(flow_model)

	flow_model = flow_svd$u[ , flow_svd$d > flow_cells_min_eigen ]
	flow_model = data.frame(flow_model)

	# Combine models and estimate surrogate covariates
	# taking into account the known covariates and the phenotype--to maintain.
	# ------------------------------------------------------------
	message("covar model size: ", ncol(covar_model))
	message("flow model size: ", ncol(flow_model))
	message("pheno model size: ", ncol(pheno_model))

	null_model = model.matrix(~., data=cbind(covar_model, flow_model))
	# null_model[is.na(null_model)] = 0

	# Remove linearly dependent columns of model matrix
	null_model = fullRankSubmat(null_model)


	# Identify surrogate model
	all_model = model.matrix(~., data=cbind(pheno_model, covar_model, flow_model))
	# all_model[is.na(all_model)] = 0  # zero remaining missing values

	all_model = fullRankSubmat(all_model)  # intercept removal

	# Run sva to identify surrogate variables
	message("Identifying surrogate variables")
	latent_model = sva(mat,
		mod=all_model,
		mod0=null_model,
		n.sv=nsurrogate_vars
	)

	# Plot surrogate variable againts STARNET ID 
	pdf(paste0("norm/plots/sur_var_", tissue, ".pdf"))
	par(mfrow=c(2, 2))
	# for (i in 1:min(25, ncol(latent_model$sv))) {
	for (i in 1:nsurrogate_vars) {
		plot(pheno_matched$starnet.ID, latent_model$sv[,i],
			xlab="STARNET ID", ylab="Surrogate var.")
	}
	dev.off()


	# Adjust regression matrix using penalized regression with known and identified covariates.
	# ----------------------------------------------------------

	# Construct adjustement model matrix
	adjust_model = model.matrix(~ 0 + .,  # no intercept, intercept assumed by penalized()
		data=cbind(covar_model, flow_model, latent_model$sv))
	adjust_model[is.na(adjust_model)] = 0
	fullRankSubmat(adjust_model)

	# rescale age to [0, 1] to have more balanced priors
	age_col = which(colnames(adjust_model) == "age")
	adjust_model[, age_col] = adjust_model[, age_col] / 100

	# L2 penalized regression model.
	# log caputes excessive newline output
	message("Fitting penalized regresssion models with ", ncol(adjust_model), " parameters.")
	log = capture.output({
		fits = apply(mat, 1, function(x) {
			penalized(x, adjust_model, lambda2=adjust_l2_penalty)
		})
	})

	# Get residuals for regression fits
	residual_mat = sapply(fits, function(x) x@residuals)

	# Get intercepts
	intercepts = sapply(fits, function(x) coefficients(x)[1])

	# plot(density(residual_mat))
	# add intercepts to adjusted mat to maintaine absolute gene expression
	adjust_mat = sweep(t(residual_mat), 1, intercepts, "+")

	message("Adjusted data, correlation with original: ",
		cor(as.vector(adjust_mat), as.vector(mat)))

	return(adjust_mat)

	# Interactive plots of the effects of adjustment on correlations.
	#
	# ---------------------------------------------------------------


	# Get linear regression coefficients
	# coefs = sapply(fits, coefficients)

	# Test to see if correlational structure is maintained
	idx = sample(nrow(mat), 1000)
	# idx = 1:1000

	cmat = cor(t(mat[idx, ]))
	cmat_adj = cor(t(adjust_mat[idx, ]))

	# Hexplot comparing correlations with adjusted correlations with beta calibration curves.
	betas = c(1.0, 1.5, 2.0, 2.5)

	# Linear for of correlation and adjusted correlations
	fit = lm(cmat[lower.tri(cmat)]~cmat_adj[lower.tri(cmat_adj)])

	pdf("norm/plots/adjust/aor_1000genes_cor_adj.pdf")

	hb = hexbin(cmat[lower.tri(cmat)], cmat_adj[lower.tri(cmat_adj)],
		xbins=100)

	p = plot(hb,
		xlab="r",
		ylab="Adjusted r",
		main=paste("R^2=", format(summary(fit)$r.squared, digits=3),
			"beta=", paste(betas, collapse=",")),
		colramp=colorRampPalette(rev(brewer.pal(11, 'Spectral'))))

	pushHexport(p$plot.vp)

	# beta calibration lines
	x = seq(0, 1, length.out=50)

	for (i in 1:length(betas)) {
		beta = betas[i]
		line_col = gray.colors(length(betas))[i]

		grid.lines(x, x^beta, gp=gpar(col=line_col, lwd=2.0), default.units="native")
		grid.lines(-x, -x^beta, gp=gpar(col=line_col, lwd=2.0), default.units="native")
	}
	upViewport()
	dev.off()


	# Density estimate plots
	pdf("norm/plots/adjust/aor_1000genes_cor_distr.pdf", width=4, height=3)
	cols = c("black", brewer.pal(9, "Set1")[1])
	plot(density(cmat_adj[lower.tri(cmat_adj)]),
		lwd=2, col=cols[2],
		xlab="r")
	lines(density(cmat[lower.tri(cmat)]), lwd=2, col=cols[1])
	legend("topright", legend=c("Ori.", "Adj."), col=cols, pch=15)
	dev.off()


	# Boxplot comparison of 
	pdf("norm/plots/adjust/aor_1000genes_boxplot_r05.pdf", width=4, height=5)
	cor_thresh = 0.5
	par(bty="n")
	boxplot(
		list(
			Ori.=cmat[lower.tri(cmat) & cmat <= -cor_thresh],
			Adj.=cmat_adj[lower.tri(cmat_adj) & cmat <= -cor_thresh],
			Ori.=cmat[lower.tri(cmat) & cmat > cor_thresh],
			Adj.=cmat_adj[lower.tri(cmat_adj) & cmat > cor_thresh]
		),
		las=2,
		cex=0.2,
		ylim=c(-1, 1),
		# col=brewer.pal(3, "Set1")[1:2]
		# col=gray.colors(2)
		col=brewer.pal(6, "Paired")[c(2, 1, 6, 5)],
		ylab="r"
	)
	abline(h=0, lty=2, col=gray.colors(5)[4])
	dev.off()

	# Find genes that correlate among the subsampled genes
	high_cor = which(abs(cmat) > 0.8, arr.ind=TRUE)
	head(high_cor, 50)

	high_cor_adj = which(abs(cmat_adj) > 0.7, arr.ind=TRUE)
	head(high_cor_adj, 50)

	# n = 2231
	# m = 35

	n = idx[491]
	m = idx[17]

	# pch = as.numeric(factor(pheno_matched$Sex))
	pdf("norm/plots/adjust/aor_1000genes_example.pdf", width=3)
	colors = c(brewer.pal(9, "Set1"), brewer.pal(8, "Set2"), brewer.pal(12, "Set3"))

	pch = as.numeric(factor(covar_matched$lab)) + 15
	cex = (pheno_matched$Age - 30) / 80
	col = colors[as.numeric(factor(covar_matched$flowcell))]

	# range(pheno_matched$Age )
	n_name = strsplit(rownames(mat)[n], "_")[[1]][1]
	m_name = strsplit(rownames(mat)[m], "_")[[1]][1]

	# Correlation estimates
	cor_adj = cor.test(adjust_mat[n, ], adjust_mat[m, ])
	cor_ori = cor.test(mat[n, ], mat[m, ])

	par(mfrow=c(2, 1), bty="o")
	plot(mat[n, ], mat[m, ],
		pch=pch,
		cex=cex,
		col=col,
		xlab=paste(n_name, expression(log(x + 1))),
		ylab=paste(m_name, expression(log(x + 1))),
		main=paste0("Original, r=", format(cor_ori$estimate, digits=3))
	)

	plot(adjust_mat[n, ], adjust_mat[m, ],
		pch=pch,
		cex=cex,
		col=col,
		xlab=paste(n_name, expression(log(x + 1))),
		ylab=paste(m_name, expression(log(x + 1))),
		# ylab=expression(log(x + 1) * (adjusted)),
		main=paste0("Adjusted, r=", format(cor_adj$estimate, digits=3))
	)
	dev.off()

})
names(expr_mats_batch) = names(expr_mats_norm)

lapply(expr_mats_batch, dim)

# Write normalized matrices as .tsv files
dir.create(file.path(data_dir, "STARNET/gene_exp_norm_batch"))
for (i in 1:length(expr_mats_batch)) {
	message("writing ", names(expr_mats_batch)[i])

	write.table(expr_mats_batch[[i]],
		file.path(data_dir, "STARNET/gene_exp_norm_batch", 
			paste0(names(expr_mats_batch)[i], ".mat")),
		sep="\t",
		quote=FALSE, col.names=NA
	)
}
# save(expr_mats_batch, file=file.path(data_dir, "STARNET/gene_exp_norm_batch/all.RData"))
save(expr_mats_batch, file=file.path(data_dir, "STARNET/gene_exp_norm_batch/all.RData"))

# OLD 
# ------------------------------

# load(file.path(data_dir, "STARNET/gene_exp_norm_batch/all.RData"))
# expr_mats_batch = lapply(expr_mats_batch, t)

# adj_mat = scale(t(mat[1:2000,]))

# heatmap(fit$coef)
# fit = rlm(t(mat[1:2000, ]) ~ adjust_model)

# adj_mat = fit$residuals
# adj_mat = scale(adj_mat)

# x = coefs
# x[is.na(x)] = 0
# heatmap.2(x,
# 	col=colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
# 	breaks=seq(-4, 4, length.out=100 + 1),
# 	trace="none")

# library(gplots)
# heatmap.2(adj_mat,
# 	col=colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
# 	breaks=seq(-4, 4, length.out=100 + 1),
# 	trace="none")


# # Batch correction using
# # Batch correction, using read lengths 50 and 100 bp, along with flowcell.
# # Technical batch correction and adjustments.
# # Also filters transcripts based on standard deviation
# # ------------------------------------------------------
# expr_mats_combatch = lapply(expr_mats_norm, function(mat) {

# 	# Match technical covariates
# 	covar_matched = covar[match(colnames(mat), covar$id), ]

# 	# Identify batches
# 	batch = covar_matched$read_length
# 	if (all(is.na(batch))) {
# 		batch[is.na(batch)] = 100  # default read length if all are missing
# 	}
# 	batch[is.na(batch)] = median(batch, na.rm=T)  # impute missing batches
# 	batch = factor(batch)

# 	# subj_range = range(as.numeric(covar$subject), na.rm=TRUE)
# 	# # batch = cut(as.numeric(covar_matched$subject), 20)
# 	# batch = cut(as.numeric(covar_matched$subject), seq(subj_range[1], subj_range[2], length.out=11))
# 	# batch[is.na(batch)] = levels(batch)[10]

# 	# Correct batch effects
# 	# Model matrix taking into acount primariy covariates
# 	options(na.action="na.pass")  # keep NA rows in model matrix

# 	# modcombat = model.matrix(~1 + flowcell, data=covar_matched)

# 	# Default model
# 	modcombat = model.matrix(~1, data=covar_matched)

# 	tryCatch({
# 		flow_model = multiCatModelMatrix(covar_matched$flowcell)

# 		# Remove flows with few samples
# 		flow_model = flow_model[, apply(flow_model, 2, sum) > 5]

# 		# Remove linear combinations of the batch and flow model
# 		dependent_col = fixDependence(
# 			data.matrix(as.numeric(batch)),
# 			data.matrix(flow_model))

# 		flow_model = flow_model[, !1:ncol(flow_model) %in% dependent_col]

# 		# Remove linear dependencies from flow_model using singular value decomposition
# 		# epsilon = 10^-6  # small value definition
# 		flow_svd = svd(flow_model)

# 		# flow_model = flow_svd$u[ , flow_svd$d > epsilon]
# 		flow_model = flow_svd$u[ , flow_svd$d > 4.0]
# 		flow_model = data.frame(flow_model)

# 		# Age and sex covariates
# 		# covar_model = covar_matched[, c("sex", "age"), with=FALSE]
# 		covar_model = covar_matched[, c("sex", "age", "subject"), with=FALSE]

# 		# Impute to median for missing covariates
# 		covar_model$sex[is.na(covar_model$sex)] = median(covar_model$sex, na.rm=TRUE)
# 		covar_model$age[is.na(covar_model$age)] = median(covar_model$age, na.rm=TRUE)
# 		covar_model$subject[covar_model$subject == "117A"] = "117"
# 		covar_model$subject[is.na(covar_model$subject)] = median(covar_model$subject, na.rm=TRUE)
# 		covar_model$subject = as.numeric(covar_model$subject)
# 		# covar_model$read_length[is.na(covar_model$read_length)] = "100"  # default

# 		# Construct covariate object from flow_model
# 		# modcombat = model.matrix(~1 + ., data=flow_model)  # sets modcombat
# 		modcombat = model.matrix(~1 + ., data=cbind(covar_model, flow_model))  # sets modcombat
# 	}, error=function(e) {
# 		warning(e)
# 		warning("1st sample: ", colnames(mat)[1])
# 	})

# 	# Batch corrected expression matrix
# 	tryCatch({
# 		batch_mat = ComBat(mat,
# 			batch=batch,
# 			mod=modcombat  # model with first column maintained and rest subtracted
# 		)
# 		return(batch_mat)
# 	}, error=function(e) {
# 		warning(e)
# 		warning("1st sample: ", colnames(mat)[1])
# 		# Use expression matrix as batch corrected data
# 		return(NA)
# 	})
# })
# lapply(expr_mats_batch, dim)

# # Write normalized matrices as .tsv files
# dir.create(file.path(data_dir, "STARNET/gene_exp_norm_batch"))
# for (i in 1:length(expr_mats_batch)) {
# 	message("writing ", names(expr_mats_batch)[i])

# 	write.table(expr_mats_batch[[i]],
# 		file.path(data_dir, "STARNET/gene_exp_norm_batch", names(expr_mats_batch)[i]),
# 		sep="\t",
# 		quote=FALSE, col.names=NA
# 	)
# }
# save(expr_mats_batch, file=file.path(data_dir, "STARNET/gene_exp_norm_batch/all.RData"))



# # Collapse normalized gene expression matrices. Does not use the batch correction.
# # -----------------------------------------------------------
# # Rename row names to gene symbols, aggregate by summing multiple matches to gene symbol.
# # assumes that transcripts have associated gene symbols
# expr_mats_collapsed = lapply(expr_mats_norm, function(mat) {
# 	# Ensembl IDs
# 	ensembl_versioned = sapply(
# 		strsplit(rownames(mat), "_"),
# 		function(vec) vec[length(vec)]  # last element
# 	)

# 	# get base ENSEMBL IDs
# 	ensembl_ids = sapply(strsplit(ensembl_versioned, "[.]"),
# 		function(x) x[1]
# 	)

# 	# rename IDs
# 	gene_symbols = mapIds(org.Hs.eg.db, keys=ensembl_ids, column="SYMBOL", keytype="ENSEMBL")

# 	mat_collapse = collapseMatSum(mat, gene_symbols)

# 	return(mat_collapse)
# })


# # Write collapsed normalized matrices
# for (i in 1:length(expr_mats_collapsed)) {
# 	message("Writing: ", names(expr_mats_collapsed)[i])

# 	write.table(expr_mats_collapsed[[i]],
# 		file.path(data_dir, "STARNET/gene_exp_norm_collapsed", names(expr_mats_collapsed)[i]),
# 		sep="\t",
# 		quote=FALSE, col.names=NA
# 	)
# }
