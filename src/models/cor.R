# Correlations between matched phenotypes and 
# mat: matrix with features in columns, such as eigengenes or 
# Returns a table of selected statistics from cor.test
phenoCorTest = function(mat, pheno_matched, phenotypes=NA) {
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
