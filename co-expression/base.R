entropy = function(x) {
	x = x[x != 0]  # exclude zero
	H = -sum(x * log10(x))
	return(H)
}

# List of counts to matrix of counts
# Note that tables must not be named, which can happen implicitly.
countMat = function(count_list) {
	mat = dcast(melt(count_list), L1~Var1, value.var="value")
	mat = mat[, 2:ncol(mat)]
	mat = as.matrix(mat)
	mat[is.na(mat)] = 0  # zero counts
	mat = t(mat)
	return(mat)
}

## Add an alpha value to a colour
addAlpha = function(col, alpha=1){
	if(missing(col))
		stop("Please provide a vector of colours.")
	apply(sapply(col, col2rgb)/255, 2, 
		function(x) 
		rgb(x[1], x[2], x[3], alpha=alpha))  
}

# Construct matrix of tissue colors by setting rgv alpha values based on frequency
freqTissueColor = function(tissue_freq, tissue_col) {
	col_mat = matrix(NA, ncol=ncol(tissue_freq), nrow=nrow(tissue_freq))
	for (i in 1:nrow(tissue_freq)) {
		for (j in 1:ncol(tissue_freq)) {
			col_mat[i, j] = addAlpha(tissue_col[i], tissue_freq[i, j])
		}
	}
	return(col_mat)
}
