rm(list=ls())

library(data.table)
library(RColorBrewer)
library(magicaxis)
library(gplots)

setwd("~/GoogleDrive/projects/STARNET/cross-tissue")


# Load data
mod_tab = fread("co-expression/tables/module_tab.csv")


cross_tissue_idx = mod_tab$purity < 0.95

tissue_col = brewer.pal(9, "Set1")[-6][c(1, 5, 6, 4, 7, 3, 2)]
tissues = c("AOR", "MAM", "LIV", "VAF", "SF", "SKLM", "BLOOD")

fdr = 0.05


# Fraction of tissues
tissue_counts = data.frame(mod_tab)[, tissues]

tissue_frac = tissue_counts / apply(tissue_counts, 1, sum)
tissue_frac = as.matrix(tissue_frac)



# Load signature enrichment data
pheno_sig_pval = fread("pheno/tables/pheno_pval.csv")



min_pval = 1e-16
cad_pheno_features = c("case_control_DEG", "syntax_score", "DUKE")
cad_pheno_sig_mat = sapply(cad_pheno_features, function(feat) {
	pval = p.adjust(pheno_sig_pval[[feat]], method="BH")
	return(pval)
})
cad_pheno_sig_mat[cad_pheno_sig_mat < min_pval] = min_pval

# cad_pheno_sig_mat = cad_pheno_sig_mat * 1  # bool to (0, 1)
rownames(cad_pheno_sig_mat) = 1:nrow(cad_pheno_sig_mat)


# min_pval = 1e-64
pheno_features = c("BMI(kg/m2)", "CRP(mg/l)", "HbA1c(%)", "Waist/Hip", "P-Chol(mmol/l)", "fP-LDL-Chol(mmol/l)", "fP-HDL-Chol(mmol/l)", "fP-TG(mmol/l)")
pheno_sig_mat = sapply(pheno_features, function(feat) {
	pval = p.adjust(pheno_sig_pval[[feat]], method="BH")
	return(pval)
})
# pheno_sig_mat = pheno_sig_mat * 1  # bool to (0, 1)
pheno_sig_mat[pheno_sig_mat < min_pval] = min_pval
rownames(pheno_sig_mat) = 1:nrow(pheno_sig_mat)





scaled_mat = cbind(tissue_frac*(-log10(min_pval))/7, -log10(cad_pheno_sig_mat), -log10(pheno_sig_mat))

hc_cross_tissue = hclust(dist(scaled_mat[cross_tissue_idx, ]))
hc_tissue_specific = hclust(dist(scaled_mat[!cross_tissue_idx, ]))

dend_lim = c(0, 50)
pdf("co-expression/annotate/plots2/dend_cross_tissue.pdf", height=3, width = plot_width/2)
plot(as.dendrogram(hc_cross_tissue), ylim=dend_lim)
dev.off()

pdf("co-expression/annotate/plots2/dend_tissue_specific.pdf", height=3, width = plot_width/2)
plot(as.dendrogram(hc_tissue_specific), ylim=dend_lim)
dev.off()



module_idx = c(
	hc_tissue_specific$labels[hc_tissue_specific$order],
	hc_cross_tissue$labels[hc_cross_tissue$order]
)

module_idx = as.integer(module_idx)



# module_idx = order(cross_tissue_idx)
plot_width = 24

	
pdf("co-expression/annotate/plots2/tissue_distribution.pdf", height=3.0, width=plot_width)
mat = t(tissue_frac[module_idx, ])
par(lwd=0.2)  # line widths for barplot borders
barplot(mat,
	col=tissue_col,
	besides=TRUE,
	# border=NA,
	space=0.0,
	# density=20,
	legend.text=colnames(tissue_frac),
	names.arg=module_idx,
	las=2,
	cex.names=0.7,
	ylab="Tissue"
	# lwd=2.3
)
dev.off()


pdf("co-expression/annotate/plots2/module_size.pdf", height=3, width=plot_width)
magplot(mod_tab$mod_size[module_idx],
	pch=16,
	log="y",
	frame.plot=FALSE,
	grid=TRUE,
	ylab="Module size"
)

void = sapply(1:length(module_idx), function(i) {
	# x0=1 is a trick to get magplot to render
	segments(i, 1,
		i, mod_tab$mod_size[module_idx][i])
}) 
dev.off()



# CAD signatures
# ----------------------------------

n_colors = 100

pdf("co-expression/annotate/plots2/CAD_sig_enrichment.pdf", width=plot_width, height=5.0)
heatmap.2(
	-log10(t(cad_pheno_sig_mat[module_idx, ])),
	Colv=NA,
	# col=c("white", brewer.pal(8, "Dark2")[1]),
	# col=colorRampPalette(c("white", brewer.pal(8, "Dark2")[1]))(n_colors),
	# col=colorRampPalette(brewer.pal(9, "Greens"))(n_colors),
	col=colorRampPalette(c("white", brewer.pal(9, "Greens")))(n_colors),
	trace="none",
	rowsep=1:3,
	sepwidth=c(0.2, 0.2)
)
dev.off()


pdf("co-expression/annotate/plots2/Pheno_sig_enrichment.pdf", width=plot_width, height=5.0)
heatmap.2(
	-log10(t(pheno_sig_mat[module_idx, ])),
	Colv=NA,
	# col=colorRampPalette(c("white", brewer.pal(8, "Dark2")[3]))(n_colors),
	col=colorRampPalette(c("white", brewer.pal(9, "Blues")))(n_colors),
	trace="none",
	# rowsep=1:3,
	sepwidth=c(0.2, 0.2)
)
dev.off()


# SYNTAX, DUKE, DEG tissue composition
# ---------------------------------------------


# features = c("syntax_score", "DUKE", "case_control_DEG")
# pheno_padj = matrix(
# 	p.adjust(
# 		data.matrix(data.frame(pheno_sig_pval)[, features]),
# 		method="BH"),
# 	ncol=3)
# colnames(pheno_padj) = c("SYNTAX", "DUKE", "Case-Ctrl DEG")


# i = 3

# fdr = 0.01

# par(mfcol=c(2, 1))
# padj = pheno_padj[, i]
# names(padj) = 1:nrow(pheno_padj)

# padj = sort(padj)

# plot(-log10(padj[padj < fdr]),
# 	main=colnames(pheno_padj)[i],
# 	pch=16,
# 	bty="n")


# mat = t(tissue_frac[which(padj < fdr), ])

# par(lwd=0.2)  # line widths for barplot borders
# barplot(mat,
# 	col=tissue_col,
# 	# besides=TRUE,
# 	# border=NA,
# 	space=0.0,
# 	# density=20,
# 	legend.text=colnames(tissue_frac),
# 	# names.arg=module_idx,
# 	las=2,
# 	cex.names=0.7,
# 	ylab="Tissue"
# 	# lwd=2.3
# )


# pheno_padj[, i] < fdr

# idx = order(pheno_padj[, i])


# j = 1
# i = 3
# plot(tissue_frac[, j], -log10(pheno_padj[, i]))
