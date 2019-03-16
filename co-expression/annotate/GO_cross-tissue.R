#
rm(list=ls())

setwd("~/Google Drive/projects/STARNET/cross-tissue")

library(GO.db)
library(org.Hs.eg.db)
library(RColorBrewer)

install.packages("lib/anRichment", repos=NULL, lib=.Library, type="source")
library(anRichment)


# Load GO p-value matrix for all modules
go_pmat = read.table("co-expression/tables/mod_GO_pmat.tsv",
	check.names=FALSE)


# Load co-expression modules
# modules = fread("co-expression/tables/modules.csv")


# # Repeat GO enrichment with anRichment package
# # -------------------------------------------------
# # Map transcripts to entrez IDs
# # Map ensembl IDs
# ensembl = sapply(strsplit(modules$ensembl, "[.]"), function(x) x[1])

# # Load map of Ensembl -> ENTREX IDs
# entrez_map = select(org.Hs.eg.db, ensembl, "ENTREZID", "ENSEMBL")

# entrez = entrez_map$ENTREZID[match(ensembl, entrez_map$ENSEMBL)]


# GOcollection = buildGOcollection(organism="human")

# GOenrichment = enrichmentAnalysis(
# 	classLabels=modules$clust,
# 	identifiers=entrez,
# 	refCollection=GOcollection,
# 	useBackground="allOrgGenes",
# 	threshold=1e-4,
# 	thresholdType="Bonferroni",
# 	getOverlapEntrez=TRUE,
# 	getOverlapSymbols=TRUE)




# Main module annotation table
mod_tab = fread("co-expression/tables/module_tab.csv")

cross_tissue = mod_tab$purity < 0.95
tissue_specific = mod_tab$purity >= 0.95

sum(cross_tissue)
sum(tissue_specific)

# Set minimum non-zero pavalues, to ensure including all log p
min_pval = min(go_pmat[go_pmat != 0])
go_pmat[go_pmat == 0] = min_pval


# Wilcoxon test on log p values for cross-tissue and tissue specific
go_tests = apply(-log10(go_pmat), 2, function(log_pval) {
	# log_pval = -log10(pval)
	wilcox_test = wilcox.test(log_pval[cross_tissue], log_pval[tissue_specific], alternative="greater")
	wilcox_test$mean_diff = mean(log_pval[cross_tissue]) - mean(log_pval[tissue_specific])  # postive is associated with cross-tissue 

	return(wilcox_test)
})

go_tests[[1]]

# pval = sapply(go_tests, function(x) x$p.value)
# mean_diff = sapply(go_tests, function(x) x$mean_diff)

go_tab = data.frame(
	term_id=names(go_tests),
	pval=sapply(go_tests, function(x) x$p.value),
	mean_diff=sapply(go_tests, function(x) x$mean_diff)
)
go_tab$term = Term(as.character(go_tab$term_id))
go_tab$ontology = Ontology(as.character(go_tab$term_id))


go_tab$padj = p.adjust(go_tab$pval)

# p.adjust(pval) < 0.05
# sum(p.adjust(pval) < 0.05, na.rm=TRUE)

fdr = 0.01
# pval_adj = p.adjust(pval)
# sig_go_tab = go_tab[na.omit(go_tab$padj < fdr), ]
sig_go_tab = go_tab[which(go_tab$padj < fdr), ]
sig_go_tab = sig_go_tab[order(sig_go_tab$pval), ]

write.csv(sig_go_tab, "co-expression/tables/cross-tissue_GO.csv", row.names=FALSE)


sig_go_tab = sig_go_tab[order(sig_go_tab$mean_diff, decreasing=TRUE), ]

head(sig_go_tab, 100)

# plot(sig_go_tab$mean_diff, -log10(sig_go_tab$padj))
# plot(sig_go_pval$mean_diff, -log10(sig_go_pval$padj))

# plot(go_tab$mean_diff, go_tab$padj)

# sig_go = data.frame(term_id=names(sig_go_pval), padj=sig_go_pval)
# sig_go$term = Term(as.character(sig_go$term_id))
# head(sig_go)
# head(sig_go, 100)

sig_go_tab$term[grep("response", sig_go_tab$term)]
sig_go_tab$term[grep("secretion", sig_go_tab$term)]

k = 50
pdf("co-expression/annotate/plots2/cross-tissue_GO_dotchart50.pdf", width=12, height=9)

k = 25
pdf("co-expression/annotate/plots2/cross-tissue_GO_dotchart25.pdf", width=10, height=6)

idx = k:1
pvals = sig_go_tab$padj[idx]
mean_diff = sig_go_tab$mean_diff[idx]

names(pvals) = sig_go_tab$term[idx]
# names(mean_diff) = sig_go_tab$term[idx]

pts_shape = rep(15, k)
pts_shape[sig_go_tab$ontology[idx] == "BP"] = 16  # circle
pts_shape[sig_go_tab$ontology[idx] == "CC"] = 17  # triangle
pts_shape[sig_go_tab$ontology[idx] == "MF"] = 15  # function

colors = brewer.pal(9, "Set1")
pts_col = rep("black", k)
pts_col[grep("response", sig_go_tab$term[idx])] = colors[1]
pts_col[grep("secretion", sig_go_tab$term[idx])] = colors[2]

par(mfrow=c(1, 2),
	bty="n"
	# xpd=TRUE
)
dotchart(-log10(pvals),
	pch=pts_shape,
	xlim=range(0, -log10(pvals)),
	xlab=expression("-log"[10] * "p (BH)"),
	ylab="Cross-tissue GO terms",
	col=pts_col,
	bty="n"
)
abline(v=-log10(fdr), col="grey")
dotchart(mean_diff,
	# pch=16,
	pch=pts_shape,
	col=pts_col,
	xlim=range(0, mean_diff),
	xlab=expression("Ave. difference (log"[10] * "p)")
)

legend("bottomright",
	legend=c("Biological Process", "Cellular Component", "Molecular Function"),
	pch=c(16, 17, 15),
	bg="white"
)
dev.off()



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



pdf("co-expression/annotate/plots2/cross_tissue_GO_qqplot.pdf")
par(mfrow=c(5, 5), mar=c(2, 2, 2, 2))

for (i in 1:25) {
	term_id = sig_go_tab$term_id[i]
	pvals = go_pmat[,colnames(go_pmat) == term_id]

	qqplotAnnot(
		-log10(pvals[tissue_specific]),
		-log10(pvals[cross_tissue]),
		xlab="Tissue-specific", ylab="Cross-tissue",
		probs=NA
	)
	title(main=sig_go_tab$term[i], cex.main=0.8, xpd=TRUE)
}
dev.off()



pheno_pval = read.table("pheno/tables/pheno_pval.csv",
	sep=",",
	header=TRUE,
	check.names=FALSE
)


features = c("syntax_score", "DUKE", "case_control_DEG")
pheno_padj = matrix(
	p.adjust(
		data.matrix(pheno_pval[, features]),
		method="BH"),
	ncol=3)
colnames(pheno_padj) = c("SYNTAX", "DUKE", "Case-Ctrl DEG")

fdr = 0.01

n_cad_assoc = apply(pheno_padj < fdr, 1, sum)
# table(apply(pheno_padj < fdr, 1, sum))


# Scatterplot of GO enrichment
# ---------------------------------------

pdf("co-expression/annotate/plots2/cross_tissue_GO_scatterplot2.pdf", width=6, height=6)
# term1 = "secretion"
term1 = "protein secretion"
# term2 = "response to external stimulus"
term2 = "response to stimulus"
i = which(Term(colnames(go_pmat)) == term1)
j = which(Term(colnames(go_pmat)) == term2)

go_bp_terms = 28912
alpha = 0.05


pts_pch = rep(1, nrow(go_pmat))
pts_pch[cross_tissue] = 16


pts_cex = (n_cad_assoc + 1) * 0.5


pts_col = rep("grey", nrow(go_pmat))

sig_idx = go_pmat[, i] < alpha / go_bp_terms &
	go_pmat[, j] < alpha / go_bp_terms

pts_col[sig_idx] = "black"



plot(-log10(go_pmat[, i]), -log10(go_pmat[, j]),
	xlab=paste(term1, "(-log10 p)"),
	ylab=paste(term2, "(-log10 p)"),
	pch=pts_pch,
	cex=pts_cex,
	col=pts_col
)

# text(-log10(go_pmat[, i]), -log10(go_pmat[, j]),
# 	labels=1:nrow(go_pmat),
# 	pos=1
# )

lab_idx = cross_tissue & 
	go_pmat[, i] < alpha / go_bp_terms &
	go_pmat[, j] < alpha / go_bp_terms &
	n_cad_assoc >=2


text(-log10(go_pmat[lab_idx, i]), -log10(go_pmat[lab_idx, j]),
	labels=(1:nrow(go_pmat))[lab_idx],
	pos=1, cex=0.7, col="red"
)


# Bonferoni lines
abline(v=-log10(alpha / go_bp_terms), lty=3, col="grey")
abline(h=-log10(alpha / go_bp_terms), lty=3, col="grey")

legend("bottomright",
	legend=c("Tissue-specific", "Cross-tissue", "CAD associations: 0", "1", "2", "3"),
	pch=c(1, 16, 16, 16, 16, 16),
	pt.cex=c(1, 1, 0.5, 1, 1.5, 2.0),
	bg="white"
)
dev.off()

# plot(-log10(go_pmat[tissue_specific, i]), -log10(go_pmat[tissue_specific, j]))
# cor.test(-log10(go_pmat[tissue_specific, i]), -log10(go_pmat[tissue_specific, j]))

