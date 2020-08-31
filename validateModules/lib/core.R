computeBlockwiseCmatTOM = function(expr, gene_blocks, beta_mat, out_dir) {

	dir.create(out_dir)

	# Adopted from detectModules.R
	#
	for (block in unique(gene_blocks)) {
		message("Block: ", block)
		block_gene_idx = gene_blocks == block
		block_tissues = expr$meta_row$tissue[block_gene_idx]

		message("Computing correlation coefficients...")
		# Calculate Pearson's co-expression correlations
		netw_mat = cor(expr$mat[, block_gene_idx])

		message("Writing matrix...")
		saveRDS(netw_mat, file.path(out_dir, paste0("cmat_block", block, ".rds")))
		gc()

		message("Computing absolute values...")
		netw_mat = abs(netw_mat)  # adjacency matrix

		message("Adjusting adjacency coefficients...")
		# Adjust weights based on specificed within- and cross-tissue beta values.
		# Loop over all tissue combinations
		for (i in 1:nrow(beta_mat)) {
			for (j in 1:ncol(beta_mat)) {
				tissue_i = rownames(beta_mat)[i]
				tissue_j = colnames(beta_mat)[j]

				# Find coordinates of the tissue combination
				idx1 = block_tissues == tissue_i
				idx2 = block_tissues == tissue_j

				# Exclude missing tissues
				idx1 = na.omit(idx1)
				idx2 = na.omit(idx2)

				if (all(!idx1) | all(!idx2)) next
				# otherwise, adjust
				netw_mat[idx1, idx2] = netw_mat[idx1, idx2] ^ beta_mat[i, j ]
			}
		}

		message("Computing topological overlap matrix...")
		# Topological overlap matrix
		# Overwrites for increased memory efficiency
		netw_mat = 1 - TOMsimilarity(netw_mat)

		message("Writing matrix...")
		saveRDS(netw_mat, file.path(out_dir, paste0("tom_block", block, ".rds")))
		gc()
	}
}

betaMat = function(within=5.2, between=2.7) {
	# Beta coefficient config 
	beta_mat = matrix(0, ncol=7, nrow=7)
	diag(beta_mat) = within
	beta_mat[lower.tri(beta_mat) | upper.tri(beta_mat)] = between

	tissues = c("SF", "VAF", "LIV", "SKLM", "BLOOD", "AOR", "MAM")
	rownames(beta_mat) = tissues
	colnames(beta_mat) = tissues

	return(beta_mat)
}