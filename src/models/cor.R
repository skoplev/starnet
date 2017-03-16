# Correlations between matched phenotypes and 
# mat: matrix with features in columns, such as eigengenes or 
# Returns a table of selected statistics from cor.test
phenoCorTest = function(mat, pheno_matched, phenotypes) {
	pheno_cor = list()
	for (phenotype in phenotypes) {
		message(phenotype)

		# Correlation tests
		cor_tests = lapply(1:ncol(mat), function(k) {
			tryCatch({
				cor_test = cor.test(mat[,k], pheno_matched[[phenotype]], method="pearson")
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
		cor_qvals = qvalue(p=cor_pvals)

		tab = data.frame(
			module=colnames(mat),
			cor=cor_coef,
			# cor_tests=cor_tests,
			pval=cor_pvals,
			qval=cor_qvals$lfdr
		)
		# pheno_cor[[phenotype]] = tab

		pheno_cor[[phenotype]] = tab
		# Order correlation tables based on p-values
		# pheno_cor[[phenotype]] = tab[order(tab$pval), ]
		# lapply(tab, function(x) {
		# 	x[order(x$pval), ]
		# })

	}
	return(pheno_cor)
}
