# Learn Bayesian networks for each module based on standardized expression data.
# Input: 
#		modules is a list of module (data.frames with colums) subset of 
#   	emat is a tissue:transcript expression matrix with corresponding meta data
learnBayesNets = function(modules, meta_genes, emat) {
	require(rcausal)

	if (nrow(meta_genes) != nrow(emat)) {
		stop("Expression matrix and transcript meta data mismatch.")
	}

	if (! "id" %in% colnames(meta_genes)) {
		stop("Transcript meta data does not contain id column")
	}

	bayes_nets = lapply(modules, function(mod) {
		# Tissue_transcipt IDs for module
		mod_ids = paste(mod$tissue, mod$transcript_id, sep="_")
		mod_gene_symbols = mod$gene_symbol

		# Print first 10 transcript IDs
		message("Module: ", paste(mod_ids[1:2], collapse=", "), "...(n=", length(mod_ids), ")")

		idx = match(mod_ids, meta_genes$id)

		# Get matrix
		submat = t(emat[idx, ])
		colnames(submat) = mod_ids

		# Standardize and impute missing to mean
		submat = scale(submat)
		submat[is.na(submat)] = 0.0
		submat = as.data.frame(submat)

		# bn = hc(submat)  # greedy hill climbing, BIC score
		bn = fges(submat, maxDegree=100)  # Fast-greedy equivalence search

		return(bn)
	})

	return(bayes_nets)
}


# Different interface than v1, based on clust vector.
# Max size of module to infer Bayesian network.
learnBayesNets.2 = function(clust, meta_genes, emat, max_size=3000) {
	require(rcausal)

	if (nrow(meta_genes) != nrow(emat)) {
		stop("Expression matrix and transcript meta data mismatch.")
	}

	if (length(clust) != nrow(meta_genes)) {
		stop("Clust vector length differ from meta_genes.")
	}

	bayes_nets = lapply(1:max(clust, na.rm=TRUE), function(k) {
		row_idx = which(clust == k)

		# Check module size
		if (length(row_idx) > max_size) {
			return(NA)
		}

		if (length(row_idx) < 10) {
			return(NA)
		}

		# Tissue_transcipt IDs for module
		mod_ids = paste(meta_genes$tissue[row_idx], meta_genes$transcript_id[row_idx], sep="_")
		# mod_gene_symbols = meta_genes$gene_symbol[row_idx]

		# Print first 10 transcript IDs
		message("Module: ", paste(mod_ids[1:2], collapse=", "), "...(n=", length(mod_ids), ")")


		# Get matrix
		submat = t(emat[row_idx,])

		# Exclude samples that have no observations
		col_idx = apply(submat, 1, function(col) {
			return(!all(is.na(col)))
		})

		submat = submat[col_idx, ]

		# Effective sample size
		message("Sample size: ", sum(col_idx))


		colnames(submat) = mod_ids

		# Standardize and impute missing to mean
		submat = scale(submat)
		submat[is.na(submat)] = 0.0
		# submat = as.data.frame(submat)

		# bn = hc(submat)  # greedy hill climbing, BIC score
		
		bn = fges(submat, maxDegree=100
		)  # Fast-greedy equivalence search


		bn = NULL

		gc()
		return(bn)
	})

	return(bayes_nets)
}
