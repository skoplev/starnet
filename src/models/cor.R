# Correlations between matched phenotypes and 
# mat: matrix with features in columns, such as eigengenes or 
# Returns a table of selected statistics from cor.test
phenoCorTest = function(mat, pheno_matched, phenotypes=NA) {
	require(qvalue)
	if (nrow(mat) != nrow(pheno_matched)) {
		stop("Dimension mismatch, row number of mat and pheno_matched should be identicial.")
	}

	if (is.na(phenotypes)) {
		phenotypes = colnames(pheno_matched)  # all phenotypes
	}

	# Name of each tested feature in mat
	feature = colnames(mat)
	if (is.null(feature)) {
		feature = rep(NA, ncol(mat))
	}

	pheno_cor = list()
	for (phenotype in phenotypes) {
		message(phenotype)

		# Correlation tests
		cor_tests = lapply(1:ncol(mat), function(k) {
			tryCatch({
				cor_test = cor.test(as.numeric(mat[,k]), as.numeric(pheno_matched[[phenotype]]), method="pearson")
				return(cor_test)
			}, error=function(e) {
				# Null
				cor_test = cor.test(c(0, 0, 0), c(0, 0, 0))
				return(cor_test)
			})
		})

		cor_pvals = sapply(cor_tests, function(x) x$p.value)
		cor_coef = sapply(cor_tests, function(x) x$estimate)

		# qvalue, local correction for multiple hypothesis testing
		cor_qvals = list()
		cor_qvals$qvalues = rep(NA, length(cor_pvals))
		tryCatch({
			cor_qvals = qvalue(p=cor_pvals)
		}, error=function(e) {
			warning(e)
		})


		tab = data.frame(
			# feature=feature,
			cor=cor_coef,
			# cor_tests=cor_tests,
			pval=cor_pvals,
			qval=cor_qvals$qvalues
		)
		# pheno_cor[[phenotype]] = tab

		rownames(tab) = feature

		pheno_cor[[phenotype]] = tab
		# Order correlation tables based on p-values
		# pheno_cor[[phenotype]] = tab[order(tab$pval), ]
		# lapply(tab, function(x) {
		# 	x[order(x$pval), ]
		# })

	}
	return(pheno_cor)
}


# Mean correlation coefficient permutation test.
# type is "rows" or "columns"
corPermuteTest = function(mat, null_mat, m=1000, type="rows") {
	require(WGCNA)  # for fast cor()
	# Test input
	if (! type %in% c("rows", "cols")) {
		stop("Invalid type specified: ", type)
	}

	if (type == "cols") {
		if (nrow(mat) != nrow(null_mat)) {
			stop("Matrix row mismatch for column perturbation")
		}
	}

	cmat = cor(
		t(mat),
		use="pairwise.complete.obs")
	cmat = abs(cmat)

	cor_coeff = cmat[lower.tri(cmat)]  # only include correlation coefficients once, and no self-correlation 
	cor_coeff_mean = mean(cor_coeff, na.rm=TRUE)

	# mat_gene_symbols = getGeneSymbols(rownames(mat))
	perm_cor_coeff_mean = sapply(1:m, function(i) {
		if (i %% 100 == 0) {
			message("Iter ", i, " out of ", m)
		}

		if (type == "rows") {
			# Pick random set of genes from null mat of same size as genes in input mat
			rand_genes = sample(
				x=nrow(null_mat),
				size=nrow(mat))

			mat_perm = null_mat[rand_genes, ]

		} else if (type == "cols") {
			rand_cols = sample(
				x=ncol(null_mat),
				size=ncol(mat))

			mat_perm = null_mat[, rand_cols]
		} else {
			stop("")  # redundant stop
		}

		cmat_perm = cor(
			t(mat_perm),
			use="pairwise.complete.obs")

		cmat_perm = abs(cmat_perm)

		cor_coeff_perm = cmat_perm[lower.tri(cmat_perm)]

		return(mean(cor_coeff_perm, na.rm=TRUE))
	})

	# Estimate p-value, or bound
	p_est = max(mean(cor_coeff_mean < perm_cor_coeff_mean), 1/m)

	return(list(
		cor_mean=cor_coeff_mean,
		cor_mean_null=perm_cor_coeff_mean,
		p_val=p_est)
	)
}
