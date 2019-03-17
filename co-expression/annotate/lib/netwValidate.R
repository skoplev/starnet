loadMarbachNetw = function(file_name) {
	tab = fread(file_name)
	names(tab) = c("from", "to", "weight")

	# From-to edge IDs, for looking up if edges are present
	tab$edge_id = paste(tab$from, tab$to, sep="_")

	return(tab)
}

# Sensitivity and precision analysis
cmpNetw = function(bnet, modules, mod_ids, rnetw, method="global") {
	# Test each STARNET subnetwork
	results = sapply(mod_ids, function(k) {
		idx = modules$clust == k
		node_ids = modules$node_id[idx]
		node_symbols = modules$gene_symbol[idx]

		if (method == "local") {
			sub_bnet = bnet[bnet$HEAD %in% node_ids & bnet$TAIL %in% node_ids, ] 
			sub_ref_netw = rnetw[rnetw$from %in% node_symbols & rnetw$to %in% node_symbols, ]
		} else if (method == "global") {
			sub_bnet = bnet[bnet$HEAD %in% node_ids | bnet$TAIL %in% node_ids, ]
			sub_ref_netw = rnetw[rnetw$from %in% node_symbols | rnetw$to %in% node_symbols, ]
		} else {
			stop("Unrecognized method: ", method)
		}

		sub_bnet$STARNET = TRUE
		sub_ref_netw$marbach = TRUE

		# Combine networks, keep all
		netw_merge = merge(
			sub_bnet[ , c("edge_id", "STARNET")],
			sub_ref_netw[, c("edge_id", "marbach", "weight")],  # V3 is the score
			by="edge_id", all=TRUE)

		netw_merge$STARNET[is.na(netw_merge$STARNET)] = FALSE  # not found in STARNET
		netw_merge$marbach[is.na(netw_merge$marbach)] = FALSE

		netw_merge$weight[is.na(netw_merge$weight)] = 0  # Marbach weight

		out = list()
		TP = sum(netw_merge$STARNET & netw_merge$marbach)  # true positive
		P = sum(netw_merge$marbach)  # positive
		FP = sum(netw_merge$STARNET & !netw_merge$marbach)

		# Sensitivity (validation rate)
		out$sensitivity =  TP / P
		out$precision = TP / (TP + FP)
		out$TP = TP

		return(out)
	})

	# Some formating
	results = t(results)
	rownames(results) = mod_ids

	results = data.matrix(as.data.frame(results))

	return(results)
}

# Convenient way of getting matrix from list of statistics
getStat = function(tests, stat) {
	i = which(colnames(tests[[1]]) == stat)
	values = sapply(tests, function(x) x[, i])
	values[is.na(values)] = 0

	values = data.matrix(data.frame(values))
	return(values)
}

