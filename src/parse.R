# Functions for parsing gene expression matrices.

# Test that the gene names are the same across matrices
testExprRownames = function(mats) {
	for (i in 2:length(mats)) {
		if (!all(rownames(mats[[1]]) == rownames(mats[[i]]))) {
			stop("Gene symbol mismatch")
		}
	}
}


# Loads normalized gene expression matrices from folder.
# Filters genes that varies less than provided standard deviation.
# Returns single data frame with data
loadNormData = function(data_dir, min_sd=0, exclude_files) {
	require(data.table)
	require(reshape2)
	require(plyr)

	# Load expression data
	expr_files = list.files(
		data_dir,
		pattern="*.mat")

	expr_files = expr_files[!expr_files %in% exclude_files]

	# Load all expression data from each tissue
	expr_mats = sapply(expr_files, function(file_name) {
		message("Loading ", file_name)

		data_table = fread(
			file.path(data_dir, file_name),
			sep="\t"
		)

		mat = data_table[, 2:ncol(data_table), with=FALSE]
		mat = as.matrix(mat)
		rownames(mat) = data_table[[1]]  # first column

		if (ncol(mat) > nrow(mat)) {
			warning("Gene expression matrix has more columns, than rows.")
		}

		return(mat)
	})

	# testExprRownames(expr_mats)

	# Prefilter genes
	if (min_sd > 0) {
		message("Filtering data")
		expr_mats = sapply(expr_mats, function(mat) {
			# Gene-wise standard deviation
			std_devs = apply(mat, 1, sd)

			return(mat[std_devs > min_sd, ])
		})
	}

	# Reshape filtered expression matrices
	message("Melting data")
	expr_melted = lapply(expr_mats, function(mat) {
		melted = melt(t(mat))  # data frame of melted expression matrix, efficient melt

		melted = rename(melted, c(
			"Var1"="sample_id",
			"Var2"="transcript_id"))

		# Get sample info
		tissue = sapply(
			strsplit(colnames(mat), "_"),
			function(x) x[1]  # 1st elem
		)
		patient_id = sapply(
			strsplit(colnames(mat), "_"),
			function(x) x[2]  # 2nd elem
		)

		# Annotate sample
		# Get sample number of melted data,
		sample_idx = match(melted[,1], colnames(mat))

		melted = data.frame(
			tissue=tissue[sample_idx],
			patient_id=patient_id[sample_idx],
			melted)

		return(melted)
	})
	names(expr_melted) = names(expr_mats)

	# Combine melted data into single data frame
	expr_melted_all = rbindlist(expr_melted)

	# Recast into table with patients in columns and tissue-specific transcripts in rows
	message("Recasting data")
	expr_recast = dcast(expr_melted_all, tissue + transcript_id ~ patient_id)

	# clean other elements from memory
	rm(expr_melted_all, expr_melted, expr_mats, scaled_mat)
	gc()  # garbage collection

	return(expr_recast)
}
