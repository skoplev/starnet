findMouseHomologue = function(genes, human_homology, mouse_homology) {
	# Map to Homologene IDs
	homologene_ids = human_homology[["HomoloGene ID"]][
		match(genes, human_homology$Symbol)
	]

	# Find matching mouse gene symbols
	mouse_symbol = mouse_homology$Symbol[
		match(homologene_ids, mouse_homology[["HomoloGene ID"]])
	]

	return(mouse_symbol)
}

# Matched correlation matrix.
# Filters out control probes, missing gene symbols, and all genes not in allowed list
matchCor = function(mat, allowed) {
	# Missing gene symbols
	exclude = which(is.na(colnames(mat)))

	# Only genes with homologues in STARNET
	exclude = c(exclude, which(!colnames(mat) %in% allowed))

	cmat = cor(mat[, ! 1:ncol(mat) %in% exclude])

	return(cmat)
}

aveAbsCor = function(cmat) mean(abs(cmat[lower.tri(cmat)]), na.rm=TRUE)

# Correlation permutation test
corPermTest = function(cmat, genes, m,
	mat_group=rep("all", ncol(expr_mat)),
	genes_group=rep("all", length(genes)))
{

	# User input test
	stopifnot(length(mat_group) == ncol(expr_mat))
	stopifnot(length(genes_group) == length(genes))

	message("n groups: ", length(levels(factor(genes_group))))

	idx = match(genes, colnames(cmat))
	message("Perm test. Found genes: ",  sum(!is.na(idx)), " / ", length(genes))

	# Include only gene groups that was found in reference matrix
	genes_group = genes_group[!is.na(idx)]

	idx = na.omit(idx)

	# Test if all groups are found in reference
	stopifnot(all(genes_group %in% mat_group))

	sub_cor = cmat[idx, idx]

	# stat = mean(abs(sub_cor[lower.tri(sub_cor)]))
	stat = aveAbsCor(sub_cor)

	stat_perm = sapply(1:m, function(i) {
		# Random sample
		# idx_rand = sample(ncol(cmat), length(idx))

		# Sample per group
		idx_rand = sapply(unique(genes_group), function(group) {
			sample(which(mat_group == group), sum(genes_group == group))
		})
		idx_rand = unlist(idx_rand)

		rand_cor = cmat[idx_rand, idx_rand]

		# rand_cor_mean = mean(abs(rand_cor[lower.tri(rand_cor)]))
		rand_cor_mean = aveAbsCor(rand_cor)
		return(rand_cor_mean)
	})

	p = mean(stat < stat_perm)
	p = max(p, 1/m)

	return(list(stat=stat, stat_perm=stat_perm, p=p, ngenes=length(idx)))
}


# Correlation test based on expression table.
# Calculates correlation coefficients when needed.
# expr_mat is formatted samples x genes (rows x cols).
# genes is general and can be used with tissue_gene IDs.
corPermTestExprMat = function(expr_mat, genes, m,
	mat_group=rep("all", ncol(expr_mat)),
	genes_group=rep("all", length(genes)))
{
	require(WGCNA)  # for fast correlation

	# User input test
	stopifnot(length(mat_group) == ncol(expr_mat))
	stopifnot(length(genes_group) == length(genes))

	message("n groups: ", length(levels(factor(genes_group))))

	# Find gene IDs in reference matrix
	idx = match(genes, colnames(expr_mat))
	message("Perm test. Found genes: ",  sum(!is.na(idx)), " / ", length(genes))

	# Include only gene groups that was found in reference matrix
	genes_group = genes_group[!is.na(idx)]

	idx = na.omit(idx)

	# Test if all groups are found in reference
	# message(genes_group[!genes_group %in% mat_group])
	stopifnot(all(genes_group %in% mat_group))

	# Calculate correlation matrix of input genes
	target_cor = cor(expr_mat[, idx], use="pairwise.complete.obs")

	message("Missing cor fraction: ", sum(is.na(target_cor))/length(target_cor))

	# Main permutation tatistic, average absolute correlation coefficient
	stat = aveAbsCor(target_cor)

	stat_perm = sapply(1:m, function(i) {
		if (i %% 50 == 0) message("iter ", i, " out of ", m)

		# Sample per group
		idx_rand = sapply(unique(genes_group), function(group) {
			sample(which(mat_group == group), sum(genes_group == group))
		})
		idx_rand = unlist(idx_rand)

		rand_cor = cor(expr_mat[, idx_rand], use="pairwise.complete.obs")

		# message("Missing rand cor fraction: ", sum(is.na(rand_cor))/length(rand_cor))

		rand_cor_mean = aveAbsCor(rand_cor)
		return(rand_cor_mean)
	})

	p = mean(stat < stat_perm)
	p = max(p, 1/m)

	return(list(stat=stat, stat_perm=stat_perm, p=p, ngenes=length(idx)))
}


makePermTable = function(perm_tests) {
	p_vals = sapply(perm_tests, function(x) x$p)
	mean_connect = sapply(perm_tests, function(x) x$stat)
	n_genes = sapply(perm_tests, function(x) x$ngenes)

	df = data.frame(names(perm_tests), p_vals, mean_connect, n_genes)
	return(df)
}


plotPermuteTest = function(test, ...) {
	hist(test$stat_perm,
		# xlim=c(0, 0.5),
		xlim=range(0, test$stat_perm, test$stat, 0.3),
		col="black",
		xlab=expression("Mean connectivity " * bar(c)),
		# main=paste0(mice, ", module ", mod),
		breaks=50,
		...
	)

	points(test$stat, 0, col=brewer.pal(9, "Set1")[1], pch=16)

	legend("topright", legend=paste("p <", test$p), bty="n")
}


# Correlation-based validation tests
# ---------------------------------------------------------
# Pairwise correlation comparisons for specific modules
# Methods:
# CT is cross-tissue correlations only, module_genes are assumed
# to encode tissue_gene.
# all: includes all genes
corModuleTest = function(mat1, mat2, module_genes, method="all") {
	require(WGCNA)

	found_idx1 = module_genes %in% rownames(mat1)
	message("Found genes mat1: ", sum(found_idx1), " / ", length(module_genes))

	found_idx2 = module_genes %in% rownames(mat2)
	message("Found genes mat2: ", sum(found_idx2), " / ", length(module_genes))

	found_both = found_idx1 & found_idx2
	message("Found genes mat1 + mat2: ", sum(found_both), " / ", length(module_genes))

	module_genes_found = module_genes[found_both]

	if (length(module_genes_found) < 2) return(NA)

	# Map found genes to rows
	idx1 = match(module_genes_found, rownames(mat1))
	idx2 = match(module_genes_found, rownames(mat2))

	# Calculate correlation matrices
	cmat1 = cor(t(mat1[idx1, ]), use="pairwise.complete.obs")
	cmat2 = cor(t(mat2[idx2, ]), use="pairwise.complete.obs")


	# Sets vectors of matching correlation coefficients
	r1 = NULL
	r2 = NULL
	if (method == "all") {
		r1 = cmat1[lower.tri(cmat1)]
		r2 = cmat2[lower.tri(cmat2)]
	} else if (method == "CT") {
		# module
		tissues = sapply(strsplit(module_genes_found, "_"), function(x) x[1])

		# Array index of lower triangular correlation matrix
		mat_index = which(lower.tri(cmat1), arr.ind=TRUE)
		mat_index = data.frame(mat_index)

		# Add tissues for each matrix index
		mat_index$row_tissue = tissues[mat_index$row]
		mat_index$col_tissue = tissues[mat_index$col]

		# Retain only matrix indices of corss-tissue correlations
		mat_index = mat_index[mat_index$row_tissue != mat_index$col_tissue, ]

		r1 = cmat1[as.matrix(mat_index[, 1:2])]
		r2 = cmat2[as.matrix(mat_index[, 1:2])]
	}

	if (length(r1) > 2) {
		cor_test = cor.test(r1, r2)
	} else {
		warning("Not enough correlation coefficient  for cor.test()")
		cor_test = NA
	}

	return(list(r1=r1, r2=r2,
		test=cor_test,
		method=method,
		ngenes=length(module_genes_found)
		# R2=cor_test$estimate^2)
		)
	)
}

plotCorTest = function(cor_test, main="", ... ) {
	require(hexbin)
	require(grid)
	hbin = hexbin(cor_test$r1, cor_test$r2,
		xbins=50
	)
	p = plot(hbin,
		colramp=colorRampPalette(rev(brewer.pal(11, "Spectral"))),
		main=paste0(main,
			" R2=", format(cor_test$test$estimate^2, digits=2),
			" p=", format(cor_test$test$p.value, digits=2)
		),
		...
	)

	# push plot viewport
	pushHexport(p$plot.vp)

	fit = lm(cor_test$r2~cor_test$r1)
	grid.abline(fit$coef[1], fit$coef[2])

	upViewport()
}
