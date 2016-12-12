# Calculate entropy based on frequency vector
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

# Horizontal error bars
errBarHoriz = function(x1, y1, x2, y2, width=0.1, col="black", lwd=1.0) {
	segments(x1, y1, x2, y2, col=col, lwd=lwd, xpd=TRUE)
	segments(x1, y1 - width, x1, y1 + width, col=col, lwd=lwd, xpd=TRUE)
	segments(x2, y1 - width, x2, y1 + width, col=col, lwd=lwd, xpd=TRUE)
}

# 
# positve specifices p
# returns correlation tests
plotPhenoCibersortCor = function(pheno_matched, cibersort_freq_matched,
	phenotype,
	ciber_i,
	k,
	positive=TRUE,
	bar_col="grey") {
	# Calculate all correlations
	freq_pheno_cor = cor(pheno_matched, cibersort_freq_matched[[ciber_i]],
		use="pairwise.complete.obs")

	freq_pheno_cor[rownames(freq_pheno_cor) == phenotype, ]

	# Pick kth highest correlated cell types
	cell_types_sel = names(
		sort(
			freq_pheno_cor[rownames(freq_pheno_cor) == phenotype, ],
			decreasing=positive
		)[1:k]
	)
	cell_types_sel = rev(cell_types_sel)  # for plotting order

	# Retest correlations with significance
	cor_test = lapply(cell_types_sel, function(cell) {
		cor.test(pheno_matched[[phenotype]], cibersort_freq_matched[[ciber_i]][[cell]],
			use="pairwise.complete.obs")
	})

	cor_coeff = sapply(cor_test, function(x) x$estimate)
	cor_conf = lapply(cor_test, function(x) x$conf.int)

	par(mar=c(4, 20, 4, 4))
	bar_centers = barplot(
		cor_coeff,
		names.arg=cell_types_sel,
		main=names(cibersort_freq_matched)[ciber_i],
		cex.main=0.7,
		horiz=TRUE,
		las=1,
		border=NA,
		col=bar_col ,
		space=1.5,
		xpd=TRUE,
		xlim=range(c(0, unlist(cor_conf))),
		xlab=paste(phenotype, " cor")
	)

	for (n in 1:k) {
		errBarHoriz(cor_conf[[n]][1], bar_centers[n], cor_conf[[n]][2], bar_centers[n])
	}

	return(cor_test)
}

# Linear color interpolation factory.
colorGradient = function(x, gradlim=NULL, colors=c("red","yellow","green"), colsteps=100) {
	pal = colorRampPalette(colors) (colsteps)  # the color palette
	if (!is.null(gradlim)) {
		return(
			pal[findInterval(x, seq(gradlim[1], gradlim[2], length.out=colsteps), all.inside=TRUE)]
		)
	} else {
		return(
			pal[findInterval(x, seq(min(x, na.rm=TRUE), max(x, na.rm=TRUE), length.out=colsteps), all.inside=TRUE)]
	 	)
	}
}

# Function to plot color bar
plotColorBar = function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') {
	scale = (length(lut)-1)/(max-min)

	# dev.new(width=1.75, height=5)
	plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
	axis(2, ticks, las=1)
	for (i in 1:(length(lut)-1)) {
		y = (i-1)/scale + min
		rect(0,y,10,y+1/scale, col=lut[i], border=NA)
 	}
}

