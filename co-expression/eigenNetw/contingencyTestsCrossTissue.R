# Module CAD association tests, including chi squared contingency tests
# Plots Venn diagram of CAD-associated modules
# ----------------------------------

rm(list=ls())

library(limma)  # for Venn
library(RColorBrewer)

setwd("~/GoogleDrive/projects/STARNET/cross-tissue")

# Load module overview table
module_tab = read.table("co-expression/tables/module_tab.csv",
	sep=",",
	header=TRUE
)

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

# Venn diagram of CAD module criteria
pdf("co-expression/eigenNetw/v3/plots/DUKE_SYNTAX_DEG_venn.pdf")
membership = pheno_padj < fdr

counts = vennCounts(membership)
vennDiagram(counts, circle.col=brewer.pal(9, "Set1"))
dev.off()

which(apply(pheno_padj < fdr, 1, sum) >= 2)
# length(which(apply(pheno_padj < fdr, 1, sum) >= 2))


# Criteria for CAD-associated modules: at least two out of three
# ------------------------------------------------


# Chi square test
# n_assoc: minimum number of associations in pheno_padj to be classified as CAD module
contingCrossTissue = function(pheno_padj, fdr, n_assoc, purity, purity_thresh_CT=0.95) {
	stopifnot(nrow(pheno_padj) == length(purity))


	cad_modules_idx = apply(pheno_padj < fdr,
		1,
		function(x) sum(x) >= n_assoc)

	cad_modules = which(cad_modules_idx)
	length(cad_modules)


	# Counts of cross-tissue modules
	cross_tissue_idx = purity < purity_thresh_CT

	# Contingency table
	cad_ct_tab = table(cross_tissue_idx, cad_modules_idx)
	print(cad_ct_tab)

	test = chisq.test(cad_ct_tab, correct=FALSE)

	return(test)
}

contingCrossTissue(pheno_padj, fdr, n_assoc=1, purity=module_tab$purity)
contingCrossTissue(pheno_padj, fdr, n_assoc=2, purity=module_tab$purity)  # result from table in paper
contingCrossTissue(pheno_padj, fdr, n_assoc=3, purity=module_tab$purity)

# Exclude MAM, AOR primary tissues that are tissue specific
tissues = c("AOR", "BLOOD", "LIV", "MAM", "SKLM", "SF", "VAF")
primary_tissue = tissues[apply(module_tab[, tissues], 1, which.max)]
table(primary_tissue)

idx = !(primary_tissue %in% c("AOR", "MAM") & module_tab$purity >= 0.95)

contingCrossTissue(pheno_padj[idx, ], fdr, n_assoc=1, purity=module_tab$purity[idx])
contingCrossTissue(pheno_padj[idx, ], fdr, n_assoc=2, purity=module_tab$purity[idx])
contingCrossTissue(pheno_padj[idx, ], fdr, n_assoc=3, purity=module_tab$purity[idx])


# Collect statistics for plotting
# without TS primary artery modules
k_range = 1:3
pvals_no_artery = sapply(k_range, function(k) {
	result = contingCrossTissue(
		pheno_padj[idx, ],
		fdr,
		n_assoc=k,
		purity=module_tab$purity[idx])
	print(result)
	return(result$p.value)
})
names(pvals_no_artery) = k_range


# All modules
pvals_all = sapply(k_range, function(k) {
	result = contingCrossTissue(
		pheno_padj,
		fdr,
		n_assoc=k,
		purity=module_tab$purity)
	print(result)
	return(result$p.value)
})
names(pvals_all) = k_range



pdf("co-expression/eigenNetw/plots/contingencyCADmodulesCT.pdf", width=3, height=4)
barplot(-log10(cbind(pvals_all, pvals_no_artery)),
	las=2,
	ylab=expression("CT / TS (-log" [10] * " p)"),
	legend.text=rep(k_range, 2),
	col=c(brewer.pal(3, "Blues"), brewer.pal(3, "Reds")),
	beside=TRUE)

abline(h=-log10(0.05), lty=2)
abline(h=0)
dev.off()