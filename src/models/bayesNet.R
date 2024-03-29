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

		bn = NULL

		# bn = hc(submat)  # greedy hill climbing, BIC score

		bn = fges(submat, maxDegree=100
		)  # Fast-greedy equivalence search

		gc()
		return(bn)
	})

	return(bayes_nets)
}


# Uses bnlearn and priors: outgoing edges required to start from TFs.
learnBayesNets.3 = function(clust, meta_genes, emat, source_symbols, max_size=3000) {
	require(bnlearn)

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
		node_ids = paste(meta_genes$tissue[row_idx], meta_genes$transcript_id[row_idx], sep="_")
		# node_ids = gsub("-", ".", node_ids)  # replaces '-' with '.' for compatability with bnlearn blacklist

		# Black list for non-TFs
		is_source = meta_genes$gene_symbol[row_idx] %in% source_symbols
		exclude_edges = expand.grid(node_ids[!is_source], node_ids)

		# Print first 10 transcript IDs
		message("Module: ", paste(node_ids[1:2], collapse=", "), "...(n=", length(node_ids), ")")


		# Get matrix
		submat = t(emat[row_idx,])

		# Exclude samples that have no observations
		col_idx = apply(submat, 1, function(col) {
			return(!all(is.na(col)))
		})

		submat = submat[col_idx, ]

		# Effective sample size
		message("Sample size: ", sum(col_idx))


		colnames(submat) = node_ids

		# Standardize and impute missing to mean
		submat = scale(submat)
		submat[is.na(submat)] = 0.0

		bn = NULL  # default

		# bn = hc(data.frame(submat), blacklist=exclude_edges)  # greedy hill climbing, BIC score

		# tetradrunner.getAlgorithmDescription(algoId='fges')
		# tetradrunner.getAlgorithmParameters(algoId='fges', scoreId='fisher-z')
		# tetradrunner.getAlgorithmParameters(algoId='fges', scoreId='sem-bic')

		# convert data frame to list of vectors (from, to) for the prior specification
		# prior = priorKnowledge(forbiddirect=unlist(apply(exclude_edges, 1, list), recursive=FALSE))

		message("Constructing priors...")
		# Silent
		invisible(capture.output(
			prior <- priorKnowledge(forbiddirect=unlist(apply(exclude_edges, 1, list), recursive=FALSE))
		))

		message("Running FGES search...")
		# Compute FGES search
		bn = tetradrunner(algoId='fges',
			df=submat,
			# scoreId='fisher-z',
			scoreId='sem-bic',
			dataType='continuous',
			alpha=0.05,
			faithfulnessAssumed=TRUE,
			maxDegree=-1,
			verbose=FALSE,
			priorKnowledge=prior
		)

		gc()
		return(bn)
	})

	return(bayes_nets)
}


# Bayesian inference using
learnBayesNets.4 = function(clust, meta_genes, emat, ref_netw, eqtl_netw, max_size=3000) {
	# require(bnlearn)

	if (nrow(meta_genes) != nrow(emat)) {
		stop("Expression matrix and transcript meta data mismatch.")
	}

	if (length(clust) != nrow(meta_genes)) {
		stop("Clust vector length differ from meta_genes.")
	}


	# Some formatting
	# ---------------------------

	# from to ids of cis-trans-eQTL network
	eqtl_netw$from = paste(eqtl_netw$tissue, eqtl_netw$gene.x, sep="_")
	eqtl_netw$to = paste(eqtl_netw$tissue, eqtl_netw$gene.y, sep="_")

	# matching tissue-ensembl IDs
	meta_genes$tissue_ensembl = paste(meta_genes$tissue, meta_genes$ensembl, sep="_")

	# Run Bayesian inference for each co-expression module
	bayes_nets = lapply(1:max(clust, na.rm=TRUE), function(k) {
		row_idx = which(clust == k)

		message("Module ", k)

		# Check module size
		if (length(row_idx) > max_size) {
			return(NA)
		}

		if (length(row_idx) < 10) {
			return(NA)
		}

		# Some ID handles on selected nodes
		node = data.frame(
			ids=meta_genes$tissue_transcript_id[row_idx],
			tissue_ensembl=meta_genes$tissue_ensembl[row_idx],
			symbols=meta_genes$gene_symbol[row_idx],
			tissue=meta_genes$tissue[row_idx]
		)

		# eQTL subnetwork with relevance for kth co-expression module
		# ----------------------------------
		eqtl_idx = eqtl_netw$from %in% node$tissue_ensembl & eqtl_netw$to %in% node$tissue_ensembl

		eqtl_netw_sub = eqtl_netw[eqtl_idx, ]


		eqtl_netw_include = data.frame(
			from=node$ids[match(eqtl_netw_sub$from, node$tissue_ensembl)],
			to=node$ids[match(eqtl_netw_sub$to, node$tissue_ensembl)]
		)

		# Remove multiples --  signals from multiple SNPs
		eqtl_netw_include = unique(eqtl_netw_include)

		message("cis-trans-eQTL constraints: ", nrow(eqtl_netw_include))

		# if (sum(eqtl_idx) > 0) {
		# 	stop("cis-trans-eQTL interaction detected, but not implemented...")
		# }


		# Interactions from known regulatory networks
		# ----------------------------------

		ref_netw_include = lapply(unique(node$tissue), function(tissue) {
			node_symbols_by_tissue = node$symbols[node$tissue == tissue]

			ref_idx = ref_netw$to %in% node_symbols_by_tissue &
				ref_netw$from %in% node_symbols_by_tissue

			netw_sub = ref_netw[ref_idx, ]

			# Translate into data frame with node IDs
			netw_sub_ids = data.frame(
				from=node$ids[match(netw_sub$from, node$symbols)],
				to=node$ids[match(netw_sub$to, node$symbols)]
			)

			return(netw_sub_ids)
		})

		ref_netw_include = rbindlist(ref_netw_include)

		message("Known regulatory interactions included: ", nrow(ref_netw_include))

		# node_ids = gsub("-", ".", node_ids)  # replaces '-' with '.' for compatability with bnlearn blacklist

		# Create white list
		# -----------------------------------

		include_edges = rbind(eqtl_netw_include, ref_netw_include)


		# Black list for non-TFs
		# is_source = meta_genes$gene_symbol[row_idx] %in% source_symbols
		# exclude_edges = expand.grid(node_ids[!is_source], node_ids)

		# Print first 10 transcript IDs
		message("Module: ", paste(node$ids[1:2], collapse=", "), "...(n=", length(node$ids), ")")


		# Get matrix
		submat = t(emat[row_idx,])

		# Exclude samples that have no observations
		col_idx = apply(submat, 1, function(col) {
			return(!all(is.na(col)))
		})

		submat = submat[col_idx, ]

		# Effective sample size
		message("Sample size: ", sum(col_idx))


		colnames(submat) = node$ids

		# Standardize and impute missing to mean
		submat = scale(submat)
		submat[is.na(submat)] = 0.0

		bn = NULL  # default

		# bn = hc(data.frame(submat), blacklist=exclude_edges)  # greedy hill climbing, BIC score

		# tetradrunner.getAlgorithmDescription(algoId='fges')
		# tetradrunner.getAlgorithmParameters(algoId='fges', scoreId='fisher-z')
		# tetradrunner.getAlgorithmParameters(algoId='fges', scoreId='sem-bic')

		# convert data frame to list of vectors (from, to) for the prior specification
		# prior = priorKnowledge(forbiddirect=unlist(apply(exclude_edges, 1, list), recursive=FALSE))

		message("Constructing priors...")

		# Silent
		invisible(capture.output(
			prior <- priorKnowledge(requiredirect=unlist(apply(include_edges, 1, list), recursive=FALSE))
			# prior <- priorKnowledge(forbiddirect=unlist(apply(exclude_edges, 1, list), recursive=FALSE))
		))

		message("Running FGES search...")
		# Compute FGES search
		bn = tetradrunner(algoId='fges',
			df=submat,
			# scoreId='fisher-z',
			scoreId='sem-bic',
			dataType='continuous',
			alpha=0.05,
			faithfulnessAssumed=TRUE,
			maxDegree=100,
			verbose=FALSE,
			priorKnowledge=prior
		)

		bn$ref_netw_include = ref_netw_include
		bn$eqtl_netw_include = eqtl_netw_include

		gc()
		return(bn)
	})

	return(bayes_nets)
}



# no priors, Tetrad on selection of co-expression modules
learnBayesNets.5 = function(clust, meta_genes, emat, max_size=3000) {
	# require(bnlearn)

	if (nrow(meta_genes) != nrow(emat)) {
		stop("Expression matrix and transcript meta data mismatch.")
	}

	if (length(clust) != nrow(meta_genes)) {
		stop("Clust vector length differ from meta_genes.")
	}

	bayes_nets = lapply(na.omit(unique(clust)), function(k) {
		row_idx = which(clust == k)

		# Check module size
		if (length(row_idx) > max_size) {
			return(NA)
		}

		if (length(row_idx) < 10) {
			return(NA)
		}

		# Tissue_transcipt IDs for module
		node_ids = paste(meta_genes$tissue[row_idx], meta_genes$transcript_id[row_idx], sep="_")
		# node_ids = gsub("-", ".", node_ids)  # replaces '-' with '.' for compatability with bnlearn blacklist

		# Print first 10 transcript IDs
		message("Module: ", paste(node_ids[1:2], collapse=", "), "...(n=", length(node_ids), ")")


		# Get matrix
		submat = t(emat[row_idx,])

		# Exclude samples that have no observations
		col_idx = apply(submat, 1, function(col) {
			return(!all(is.na(col)))
		})

		submat = submat[col_idx, ]

		# Effective sample size
		message("Sample size: ", sum(col_idx))


		colnames(submat) = node_ids

		# Standardize and impute missing to mean
		submat = scale(submat)
		submat[is.na(submat)] = 0.0

		# bn = NULL  # default
		# bn = hc(data.frame(submat))  # greedy hill climbing, BIC score

		# bn = hc(data.frame(submat), blacklist=exclude_edges)  # greedy hill climbing, BIC score

		# tetradrunner.getAlgorithmDescription(algoId='fges')
		# tetradrunner.getAlgorithmParameters(algoId='fges', scoreId='fisher-z')
		# tetradrunner.getAlgorithmParameters(algoId='fges', scoreId='sem-bic')

		# convert data frame to list of vectors (from, to) for the prior specification
		# prior = priorKnowledge(forbiddirect=unlist(apply(exclude_edges, 1, list), recursive=FALSE))

		# message("Constructing priors...")
		# # Silent
		# invisible(capture.output(
		# 	prior <- priorKnowledge(forbiddirect=unlist(apply(exclude_edges, 1, list), recursive=FALSE))
		# ))



		message("Running FGES search...")
		# Compute FGES search
		bn = tetradrunner(
			algoId='fges',
			df=submat,
			# scoreId='fisher-z',
			scoreId='sem-bic',
			dataType='continuous',
			# alpha=0.05,
			faithfulnessAssumed=TRUE,
			# maxDegree=-1,
			maxDegree=100,
			# numberResampling = 10,
			# resamplingEnsemble = 2,
			verbose=FALSE
			# priorKnowledge=prior
		)

		gc()
		return(bn)
	})
	names(bayes_nets) = na.omit(unique(clust))

	return(bayes_nets)
}

