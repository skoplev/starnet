# Remove GO terms that are prevalent across modules
removePrevalentTerms = function(enrich_tab, max_preval) {
	mod_counts = table(enrich_tab$termID)
	include_idx = mod_counts < max_preval
	message("Excluding GO terms: ", sum(!include_idx))
	include_terms = names(which(include_idx))
	return(enrich_tab[enrich_tab$termID %in% include_terms, ])
}

# Filter out GO terms from WGCNA GO enrichment analysis
filterCommonGOTerms = function(go, max_preval=50) {
	go$bestPTerms$BP$enrichment = removePrevalentTerms(go$bestPTerms$BP$enrichment, max_preval)
	go$bestPTerms$CC$enrichment = removePrevalentTerms(go$bestPTerms$CC$enrichment, max_preval)
	go$bestPTerms$MF$enrichment = removePrevalentTerms(go$bestPTerms$MF$enrichment, max_preval)
	return(go)
}

enrichmentGO = function(modules) {
	require(org.Hs.eg.db)
	require(WGCNA)
	# Map ensembl IDs
	ensembl = sapply(strsplit(modules$meta_genes$ensembl, "[.]"), function(x) x[1])

	# Load map of Ensembl -> ENTREX IDs
	entrez_map = select(org.Hs.eg.db, ensembl, "ENTREZID", "ENSEMBL")

	entrez = entrez_map$ENTREZID[match(ensembl, entrez_map$ENSEMBL)]

	go_enrich = GOenrichmentAnalysis(between$clust,
		entrez,
		# between$meta_genes$ensembl_nosuf,
		organism="human",
		removeDuplicates=FALSE,
		nBestP=100,
		# pCut=0.05  # doesn't work
	)

	return(go_enrich)
}

hyperGeometricModuleTest = function(env, genes) {
	require(qvalue)
	mod_stats = countModuleTissueStat(env)
	gene_bool = env$meta_genes$gene_symbol %in% genes
	# sum(gene_bool)

	# Aggregate gene symbols by module
	genes_module = sapply(1:max(env$clust), function(k) {
		idx = env$clust == k & gene_bool

		if (sum(idx) == 0) {
			return(c())
		} else {
			return(env$meta_genes[idx, ])
		}
	})

	# Count number of found GWAS-associated genes
	tab = data.frame(n_genes=sapply(genes_module, function(df) max(0, nrow(df))))
	tab$genes = sapply(genes_module, function(df) paste(df$gene_symbol, collapse=";"))

	# Hypergeometric test of CAD-associated transcripts in each module
	m_mRNA = sum(env$meta_genes$gene_symbol %in% genes)
	n_non_mRNA = sum(! env$meta_genes$gene_symbol %in% genes)

	p_hyper = sapply(1:max(env$clust), function(k) {
		module_size = mod_stats$size[k]
		module_mRNA = max(0, nrow(genes_module[[k]]))

		p = 1 - phyper(module_mRNA, m_mRNA, n_non_mRNA, module_size)
		return(p)
	})

	tab$pval = p_hyper
	try({
		tab$qvalue = qvalue(tab$pval)$qvalue
	})

	return(tab)
}

# Multiple hypergeometric tests for enrichment in co-expression modules.
# clust_genes is a vector of gene IDs.
# clust is a vector of integers encoding the assigned clusters
# genes is the set of genes tested
hyperGeometricModuleTestCore = function(clust_genes, clust, genes) {
	require(qvalue)
	stopifnot(length(clust_genes) == length(clust))

	# mod_stats = countModuleTissueStat(env)
	gene_bool = clust_genes %in% genes
	# sum(gene_bool)

	# Calculate module sizes
	clust_size = table(clust)

	# Aggregate gene symbols by module
	genes_module = sapply(1:max(clust), function(k) {
		idx = clust == k & gene_bool

		if (sum(idx) == 0) {
			return(c())
		} else {
			return(clust_genes[idx])
		}
	})

	# Count number of found GWAS-associated genes
	tab = data.frame(n_genes=sapply(genes_module, function(genes) max(0, length(genes))))
	tab$genes = sapply(genes_module, function(genes) paste(genes, collapse=";"))

	# Hypergeometric test of CAD-associated transcripts in each module
	m_mRNA = sum(clust_genes %in% genes)
	n_non_mRNA = sum(!clust_genes %in% genes)

	p_hyper = sapply(1:max(clust), function(k) {
		module_size = clust_size[k]
		module_mRNA = max(0, length(genes_module[[k]]))

		p = 1 - phyper(module_mRNA, m_mRNA, n_non_mRNA, module_size)
		return(p)
	})

	tab$pval = p_hyper
	try({
		tab$qvalue = qvalue(tab$pval)$qvalue
	})

	return(tab)
}
