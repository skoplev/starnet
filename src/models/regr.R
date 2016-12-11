# Fits multivariate linear model for each phenotype based the eigengenes
# Input matrices are required to be sample matched.
# Samples in rows, features in columns
fitLinearEigenPheno = function(
	phenotype,
	module_eigengenes,
	exclude_phenotypes=c("starnet.ID", "id"))
{
	phenotype = as.data.frame(phenotype)  # if data.table
	module_eigengenes = as.data.frame(module_eigengenes)

	stopifnot(nrow(phenotype) == nrow(module_eigengenes))

	# Exclude phenotype columns
	phenotype = phenotype[, !colnames(phenotype) %in% exclude_phenotypes]

	# Fit linear models of eigengenes vs phenotype.
	fits = list()
	for (pfeat in colnames(phenotype)) {
		message("Fitting ", pfeat)

		# Make target~input data frame
		df = cbind(phenotype[[pfeat]], module_eigengenes)
		if (length(levels(factor(df[,1]))) == 2) {
			# {-1, 1} dummy encoding for binomial categorical variables
			df[,1] = as.integer(factor(df[,1])) * 2 - 3
		}

		# exclude samples with missing target variable
		df = df[!is.na(df[,1]), ]
		df = data.frame(df)  # ensure same length? (don't know why it is necessary)

		tryCatch({
			# Linear fit formula comparing
			first_rest_form = as.formula(
				paste(colnames(df)[1], "~", 
					paste(colnames(df)[2:ncol(df)], collapse="+")
				)
			)

			# first vs rest linear fit
			fits[[pfeat]] = lm(first_rest_form, data=df)
		}, error=function(e) {
			# Warning for failed linear fits
			warning(pfeat, " excluded.")
		})
	}

	return(fits)
}

# Get pmat matrix from list of linear fits
getPmat = function(fits, remove_intercept=TRUE) {
	# number of parameters
	pars = sapply(fits, function(fit) {
		nrow(summary(fit)$coef)
	})

	# exclude fits without the correct number of coefficients
	# 
	warnings("Excluding ", names(which(pars != median(pars))))
	fits = fits[pars == median(pars)]

	pmat = sapply(fits, function(fit) {
		summary(fit)$coef[,4]
	})

	if (remove_intercept) {
		pmat = pmat[!rownames(pmat) %in% "(Intercept)", ]
	}
	return(pmat)
}


# filter matrix such that it contains at least one row or column lower than alpha
# Used for p-value matrix filtering.
filterMatRowColMin = function(mat, alpha) {
	mat = as.matrix(mat)
	include_row = apply(mat, 1, min, na.rm=TRUE) < alpha
	include_col = apply(mat, 2, min, na.rm=TRUE) < alpha
	return(mat[include_row, include_col])
}
