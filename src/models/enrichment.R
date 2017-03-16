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
