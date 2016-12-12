# Parses STARNET count matrices.
# Normalizes using DESeq2 size factors, writes matrices to
# Collapses by HUGO gene symbols using sum for multiple mappings
# to the same gene symbol.

rm(list=ls())

library(data.table)
library(org.Hs.eg.db)
library(DESeq2)

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

	model_mat = matrix(0, ncol=length(catvec_uni), nrow=ncol(mat))
	colnames(model_mat) = catvec_uni
	for (i in 1:length(catvec_sep)) {
		model_mat[i, ] = catvec_uni %in% catvec_sep[[i]]
	}

	multiplicity = apply(model_mat, 1, sum)
	multiplicity[multiplicity == 0] = 1  # avoid division by zero
	model_mat = sweep(model_mat, 1, multiplicity, "/")

	return(data.frame(model_mat))
}


data_dir = "/Users/sk/DataProjects/cross-tissue"

# Directory containing expression count matrices
emat_dir = file.path(data_dir, "STARNET/gene_exp/matrices")


# gene symbol, exp
expr_files = list.files(emat_dir)

# Load all data matrices
expr_mats = lapply(expr_files, function(file_name) {
	d = fread(
		file.path(emat_dir, file_name))
	return(d)
})
names(expr_mats) = expr_files

# Load phenotype and covariate data for batch correction

pheno = fread(file.path(
	"/Volumes/SANDY/phenotype_data",
	"STARNET_main_phenotype_table.cases.Feb_29_2016.tbl"
))

covar = fread(file.path(
	"/Volumes/SANDY/phenotype_data",
	"covariates.cases_and_controls.April_12_2016.txt"
))


# Normalize using size factors from DESeq2
# ---------------------------------------------------------------
# converts $id to rownames of returned matrices
expr_mats_norm = sapply(expr_mats, function(emat) {
	# get numerical matrix
	mat = emat[,2:ncol(emat), with=FALSE]
	rownames(mat) = emat$id

	size_factors = estimateSizeFactorsForMatrix(as.matrix(mat))

	# divide each column by size factor
	norm_mat = sweep(mat, 2, size_factors, "/")

	return(norm_mat)
})

# Write normalized matrices as .tsv files
for (i in 1:length(expr_mats_norm)) {
	message("writing", names(expr_mats_norm)[i])

	write.table(expr_mats_norm[[i]],
		file.path(data_dir, "STARNET/gene_exp_norm", names(expr_mats_collapsed)[i]),
		sep="\t",
		quote=FALSE, col.names=NA
	)
}
# save(expr_mats_norm, file=file.path(data_dir, "STARNET/gene_exp_norm/all.RData"))
# load(file.path(data_dir, "STARNET/gene_exp_norm/all.RData"))

# Batch correction, using read lengths 50 and 100 bp, along with flowcell.
# Technical batch correction and adjustments.
# Also filters transcripts based on standard deviation
# ------------------------------------------------------
# Variance filter and batch correction
min_sd = 0.5
expr_mats_batch = lapply(expr_mats_norm, function(mat) {
	# Filter genes based on standard deviation
	mat = mat[apply(mat, 1, sd, na.rm=T) > min_sd, ]

	# Match phenotype data to selected gene expression matrix
	patient_ids = sapply(
		strsplit(colnames(mat), "_"),
		function(x) x[2]
	)

	pheno_matched = pheno[match(patient_ids, pheno$starnet.ID), ]
	covar_matched = covar[match(colnames(mat), covar$sample), ]


	# Identify batches
	batch = covar_matched$read_length
	if (all(is.na(batch))) {
		batch[is.na(batch)] = 100  # default read length if all are missing
	}
	batch[is.na(batch)] = median(batch, na.rm=T)  # 

	# Correct batch effects
	# Model matrix taking into acount primariy covariates
	options(na.action="na.pass")  # keep NA rows in model matrix

	# modcombat = model.matrix(~syntax_score + BMI + LDL + Age, data=pheno_matched)
	modcombat = model.matrix(~1 + flowcell, data=covar_matched)

	flow_model = multiCatModelMatrix(covar_matched$flowcell)

	# Remove flows with few samples
	flow_model = flow_model[, apply(flow_model, 2, sum) > 5]

	# Remove linear dependencies from flow_model using singular value decomposition
	epsilon = 10^-6  # small value definition
	flow_svd = svd(flow_model)

	flow_model = flow_svd$u[ , flow_svd$d > epsilon]
	flow_model = data.frame(flow_model)

	# # Impute model input to median
	# modcombat[is.na(modcombat[, 2]), 2] = median(modcombat[, 2], na.rm=T)  # syntax score
	# modcombat[is.na(modcombat[, 3]), 3] = median(modcombat[, 3], na.rm=T)  # BMI
	# modcombat[is.na(modcombat[, 4]), 4] = median(modcombat[, 4], na.rm=T)  # BMI
	# modcombat[is.na(modcombat[, 5]), 5] = median(modcombat[, 5], na.rm=T)  # BMI

	# modcombat = model.matrix(~1, data=pheno_matched)

	# Remove linear dependencies between batch and flow_model
	flow_model = flow_model[, cor(as.numeric(batch), flow_model) < 1 - epsilon]

	# Construct covariate object from flow_model
	modcombat = model.matrix(~1 + ., data=flow_model)

	# all_model = cbind(batch, flow_model)
	# mod = model.matrix(~1 + ., data=all_model)

	# Batch corrected expression matrix
	tryCatch({
		batch_mat = ComBat(mat,
			batch=batch,
			mod=modcombat  # model with first column maintained and rest subtracted
		)
		# batch_mat = sva(mat,
		# 	# batch=batch,
		# 	mod=mod
		# )

		return(batch_mat)
	}, error=function(e){
		warning(e)
		# Use expression matrix as batch corrected data
		return(mat)
	})
})

# Write normalized matrices as .tsv files
dir.create(file.path(data_dir, "STARNET/gene_exp_norm_batch"))
for (i in 1:length(expr_mats_batch)) {
	message("writing ", names(expr_mats_batch)[i])

	write.table(expr_mats_batch[[i]],
		file.path(data_dir, "STARNET/gene_exp_norm_batch", names(expr_mats_batch)[i]),
		sep="\t",
		quote=FALSE, col.names=NA
	)
}
save(expr_mats_batch, file=file.path(data_dir, "STARNET/gene_exp_norm_batch/all.RData"))



# Correct for clinical covariates
# --------------------------------------------------------------
lapply(expr_mats_batch, function(mat) {
	
})


# Collapse normalized gene expression matrices. Does not use the batch correction.
# -----------------------------------------------------------
# Rename row names to gene symbols, aggregate by summing multiple matches to gene symbol.
# assumes that transcripts have associated gene symbols
expr_mats_collapsed = lapply(expr_mats_norm, function(mat) {
	# Ensembl IDs
	ensembl_versioned = sapply(
		strsplit(rownames(mat), "_"),
		function(vec) vec[length(vec)]  # last element
	)

	# get base ENSEMBL IDs
	ensembl_ids = sapply(strsplit(ensembl_versioned, "[.]"),
		function(x) x[1]
	)

	# rename IDs
	gene_symbols = mapIds(org.Hs.eg.db, keys=ensembl_ids, column="SYMBOL", keytype="ENSEMBL")

	mat_collapse = collapseMatSum(mat, gene_symbols)

	return(mat_collapse)
})


# Write collapsed normalized matrices
for (i in 1:length(expr_mats_collapsed)) {
	message("Writing: ", names(expr_mats_collapsed)[i])

	write.table(expr_mats_collapsed[[i]],
		file.path(data_dir, "STARNET/gene_exp_norm_collapsed", names(expr_mats_collapsed)[i]),
		sep="\t",
		quote=FALSE, col.names=NA
	)
}
