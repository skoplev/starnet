rm(list=ls())

setwd("~/Google Drive/projects/STARNET/cross-tissue")

library(data.table)


qqplotAnnot = function(x, y,
	probs=c(0.001, 0.01, 0.1, 0.9, 0.99, 0.999),
	...)
{
	qqplot(
		x,
		y,
		pch=16,
		cex=0.8,
		xlim=range(x, y),
		ylim=range(x, y),
		...
	)
	abline(0, 1, col="red")

	q_y = quantile(y, probs=probs, na.rm=TRUE)
	q_x = quantile(x, probs=probs, na.rm=TRUE)

	points(q_x, q_y, col="grey", cex=0.6)
	text(q_x, q_y, label=paste0(probs*100, "%"), pos=1, cex=0.8, col="grey")
}

# Main module table
mod_tab = fread("co-expression/tables/module_tab.csv")

sig_pmat = fread("pheno/tables/pheno_pval.csv")
sig_pmat = as.data.frame(sig_pmat)

# Remove zero pvalues
# smallest non-zero pvalue
# min_pval = min(sig_pmat[sig_pmat != 0])
# sig_pmat[sig_pmat == 0] = min_pval


sig_pmat = apply(sig_pmat, 2, function(p) {
	min_pval = min(p[p != 0 ])
	p[p == 0] = min_pval
	return(p)
})




cross_tissue = mod_tab$purity < 0.95
tissue_specific = mod_tab$purity >= 0.95


densPlot = function(x, y, cex_legend=1.0, ...) {
	densy = density(y)
	densx = density(x)

	wilcox_test = wilcox.test(x, y)

	plot(0,
		type="n",
		xlim=range(densx$x, densy$x),
		ylim=range(densx$y, densy$y),
		bty="l",
		ylab="Density",
		...
	)
	colors = c(
		rgb(228, 26, 28, 150, maxColorValue=255),
		rgb(55, 126, 184, 150, maxColorValue=255)
	)

	polygon(densx, col=colors[1])
	polygon(densy, col=colors[2])

	legend("topright",
		legend=c("Cross-tissue", "Tissue-specific",
			paste0("p=", format(wilcox_test$p.value, digits=4))
		),
		cex=cex_legend,
		fill=c(colors, "white"))
}


pdf("pheno/plots/module_size_density.pdf", height=2.9, width=4)
x = mod_tab$mod_size[cross_tissue]
y = mod_tab$mod_size[tissue_specific]
densPlot(log10(x), log10(y),
	cex_legend=0.7,
	xlab=expression("Module size (log"[10] * " n)"))
dev.off()


# max_mod_size = 1000
# cross_tissue = mod_tab$purity < 0.95 & mod_tab$mod_size < max_mod_size
# tissue_specific = mod_tab$purity >= 0.95 & mod_tab$mod_size < max_mod_size
# sum(cross_tissue)
# sum(tissue_specific)


features = colnames(sig_pmat)[-1]
pdf("pheno/plots/CAD_QQplot.pdf", height=2.9, width=8)
features = colnames(sig_pmat)[c(2, 3, 4)]
par(mfrow=c(1, 3))
for (feat in features) {
	x = -log10(sig_pmat[tissue_specific, feat])
	y = -log10(sig_pmat[cross_tissue, feat])

	wilcox_test = wilcox.test(x, y)

	qqplotAnnot(x, y,
		xlab=expression("Tissue-specific, -log" [10] * " p"),
		ylab=expression("Cross-tissue, -log"[10] * " p"),
		probs=NA,
		main=feat)

	legend("bottomright",
		legend=paste0("p=",
			format(wilcox_test$p.value, digits=4)
		),
		bty="n"
	)
}
dev.off()
