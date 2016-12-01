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

# Normalize using size factors from DESeq2
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
