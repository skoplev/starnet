# Estimate signatures for phenotypes using regression or classification

rm(list=ls())

library(RcppEigen)  # fastLm
# library(parallel)

library(data.table)
library(metap)

# library(glmnet)
library(WGCNA)
library(RColorBrewer)

library(gplots)
library(qvalue)

library(vioplot)



data_dir = "~/DataProjects/cross-tissue"

setwd("~/Google Drive/projects/STARNET/cross-tissue")
source("src/base.R")
source("src/parse.R")

pheno = fread(
	"~/Google Drive/projects/STARNET/phenotype/data/current/STARNET_main_phenotype_table.2017_12_03.tsv"
)

# Load module table
mod_tab = fread("co-expression/tables/module_tab.csv")


# Combined
load(file.path(data_dir, "STARNET/gene_exp_norm_reshape/expr_recast.RData"),
	verbose=TRUE)

# Parse expression matrix
# TODO: use parseExprTable() instead
mat = expr_recast[, 3:ncol(expr_recast)]

meta_row = expr_recast[, 1:2]
meta_row$ensembl = sapply(strsplit(as.character(meta_row$transcript_id), "_"), function(x) x[2])
meta_row$ensembl = sapply(strsplit(meta_row$ensembl, "[.]"), function(x) x[1])
meta_row$tissue_ensembl = paste(meta_row$tissue, meta_row$ensembl, sep="_")

rownames(mat) = paste(meta_row$tissue, meta_row$transcript_id, sep="_")

tmat = as.data.frame(t(mat))
colnames(tmat) = rownames(mat)



# Match phenotype data
pheno_match = pheno[match(colnames(mat), pheno$starnet.ID), ]

# Cleanup
rm(expr_recast)
rm(mat)


# Load module data
mod = fread("co-expression/tables/modules.csv")
mod$ensembl = sapply(strsplit(mod$ensembl, "[.]"), function(x) x[1])
mod$tissue_ensembl = paste(mod$tissue, mod$ensembl, sep="_")


# Load DEG tables
deg_tissues = c("AOR", "SKLM", "LIV", "VAF", "SF")

deg = lapply(deg_tissues, function(tissue) {
	d = fread(paste0("case-control/data/deseq/deseq_full_", tissue, ".csv"))
	d$tissue = tissue
	return(d)
})
names(deg) = deg_tissues

# Use AOR statistics for MAM tissue
deg$MAM = deg$AOR
deg$MAM$tissue = "MAM"

deg_all = rbindlist(deg)
colnames(deg_all)[1] = "ensembl"

deg_all$tissue_ensembl = paste(deg_all$tissue, deg_all$ensembl, sep="_")


# Evaluate DEG signatures for all tissues


# Strict signature criteria => 48 enriched modules
fdr = 0.05
base_mean_lim = 50
fc_lim = log2(1.3)  # +- 30%


# # Less strict criteria => 83 enriched modules
# fdr = 0.05
# base_mean_lim = 10
# fc_lim = log2(1.0)



idx_sig = deg_all$padj < fdr &
	deg_all$baseMean > base_mean_lim &
	abs(deg_all$log2FoldChange) > fc_lim

# sum(idx_sig)

deg_all_sig = deg_all[idx_sig, ]


source("src/models/enrichment.R")


deg_enrichment = hyperGeometricModuleTestCore(
	clust_genes=mod$tissue_ensembl,
	clust=mod$clust,
	genes=deg_all_sig$tissue_ensembl)
# Match DEG statistics to meta_row for all tissues

# sum(deg_enrichment$qvalue < 0.2)
# sum(deg_enrichment$pval < 0.1)

deg_enrichment$padj = p.adjust(deg_enrichment$pval)

sum(deg_enrichment$padj < 0.1)





min_pval = 1e-16
fdr = 0.1
d = deg_enrichment
d$clust = 1:nrow(deg_enrichment)

d = d[, colnames(d) != "genes"]  # drop concat gene column

d = cbind(d,
	mod_tab[,
		c("mod_size",
		"purity")
	]
)

d = d[order(d$pval, decreasing=TRUE), ]

d$padj[d$padj < min_pval] = min_pval

d = d[d$padj < fdr, ]


pdf("pheno/plots/signatures/DEG_enrichment_48.pdf", width=3, height=8)
# pdf("pheno/plots/signatures/DEG_enrichment_83.pdf", width=3, height=12)
colors = c(
	rgb(210, 210, 210, maxColorValue=255),
	brewer.pal(9, "Dark2")[5])
bar_cols = rep(colors[1], nrow(d))
bar_cols[d$purity < 0.95] = colors[2]  # cross-tissue

barplot(-log10(d$padj),
	names.arg=rownames(d),
	col=bar_cols,
	border=rgb(80, 80, 80, maxColorValue=255),
	las=2,
	horiz=TRUE,
	xlab=expression("DEG enrichment, -log"[10] * " p (BH)"),
	ylab="Co-expression modules",
	cex.names=0.8
)

legend("bottomright",
	legend=c("Tissue-specific", "Cross-tissue"),
	fill=colors,
	# pch=22
	cex=0.8,
	bty="n"
)
dev.off()



# Correlation-based signatures
# ------------------------------------------------------------------
tissues = unique(meta_row$tissue)

# feature = "syntax_score"
# feature = "DUKE"
# feature = "ndv"
# feature = "lesions"

features = c(
	"syntax_score",
	"DUKE",
	"ndv",
	"lesions",
	"BMI(kg/m2)",
	"CRP(mg/l)",
	"HbA1c(%)",
	"Waist/Hip",
	"P-Chol(mmol/l)",
	"fP-LDL-Chol(mmol/l)",
	"fP-HDL-Chol(mmol/l)",
	"fP-TG(mmol/l)"
)




# Including 
# deg_all_matched = deg_all[match(meta_row$tissue_ensembl, deg_all$tissue_ensembl), ]

# Calculate per-gene correlation statistics for each feature
fits = list()
# fits$DEG_fisher$p = deg_all_matched$pvalue  # p-values from DEGs

for (feature in features) {
	message(feature)
	fits[[feature]] = corAndPvalue(tmat, pheno_match[[feature]])
}

# sumlog(c(0.1, 0.1, 0.1))


# Aggregate p-values by module
n_modules = max(mod$clust)
# mod_unique = sort(unique(mod$clust))
pmat = sapply(fits, function(fit) {
	# Aggregate p-values by cluster
	aggregate_pvalue = sapply(1:n_modules, function(k) {
		# message(k)
		idx = mod$clust == k
		pvalues = fit$p[idx]
		return(sumlog(pvalues)$p)  # metaanalysis by Fisher's method
	})

	return(aggregate_pvalue)
})


# Make combined signature p-value table. Each row corresponds to a module.
sig_tab = cbind(
	data.table(clust=1:nrow(deg_enrichment), case_control_DEG=deg_enrichment$pval),
	data.table(pmat)
)
write.csv(sig_tab, "pheno/tables/pheno_pval.csv",
	row.names=FALSE)


pdf("pheno/plots/signatures/Fisher_meta.pdf", width=12, height=20)
min_pvalue = 1e-100
par(mfrow=c(ncol(pmat)/2, 2))
for (i in 1:ncol(pmat)) {
	pvals = pmat[, i]
	pvals[pvals < min_pvalue] = min_pvalue
	qvals = qvalue(pvals)$qvalue
	names(qvals) = 1:length(qvals)

	sorted_qvals = sort(qvals)

	cols = rep(rgb(150, 150, 150, maxColorValue=255), length(qvals))
	cols[names(sorted_qvals) %in% c(78)] = brewer.pal(9, "Set1")[1]
	cols[names(sorted_qvals) %in% c(98)] = brewer.pal(9, "Set1")[2]

	barplot(-log10(sorted_qvals),
		las=2,
		# space=-0.05,
		space=0,
		main=colnames(pmat)[i],
		cex.names=0.3,
		ylab="-log10 FDR (q-value)",
		col=cols,
		border=NA
	)
}
dev.off()



pdf("pheno/plots/signatures/Fisher_meta_heatmap.pdf")
pmat_adj = qvalue(pmat)$qvalue

min_pvalue = 1e-16
pmat_adj[pmat_adj < min_pvalue] = min_pvalue

heatmap.2(
	-log10(pmat_adj),
	trace="none",
	col=colorRampPalette(brewer.pal(9, "YlOrRd"))(100),
	ylab="Module ID",
	cexRow=0.2,
	mar=c(12, 20),
	key.title="",
	key.xlab="FDR (q-value)"
)
dev.off()


vioplotTest = function(x, y, col="grey", main="", ylab="") {
	par(bty="n")
	vioplot(x, y,
		col=col,
		names=c("-", "+")
	)
	# abline(h=0.05, col="grey", lty=3)
	title(ylab=ylab,
		main=paste0(main, " P=", format(wilcox.test(x, y)$p.value, digits=4)),
		cex.main=0.5
	)
}

# Violin plots of 
pdf("pheno/plots/signatures/cross-tissue-association2.pdf", height=5, width=12)

par(mfrow=c(3, 8), mar=c(2, 4, 2, 2))

tissue_order = c("AOR", "MAM", "LIV", "VAF", "SF", "SKLM", "BLOOD")
tissue_col = brewer.pal(9, "Set1")[c(1, 5, 7, 4, 8, 3, 2)]
names(tissue_col) = tissue_order
tissue_col = as.list(tissue_col)

features = c("DUKE", "syntax_score", "case_control_DEG")

fdr = 0.05


for (feat in features) {
	sig_idx =
		p.adjust(sig_tab[[feat]]) < fdr

	# sig_idx[is.na(sig_idx)] = FALSE

	for (t in tissue_order) {
		tissue_frac = mod_tab[[t]] / mod_tab$mod_size
		x = tissue_frac[!sig_idx]
		y = tissue_frac[sig_idx]

		vioplotTest(x, y, col=tissue_col[[t]], ylab=t, main=feat)
	}


	x = 1 - mod_tab$purity[!sig_idx]
	y = 1 - mod_tab$purity[sig_idx]
	vioplotTest(x, y,
		col=brewer.pal(9, "Set1")[9], main=feat, ylab="Cross-tissue")
}
dev.off()



fdr = 0.1
library(limma)
sig_idx = 
	p.adjust(pmat[, colnames(pmat) == "DUKE"]) < fdr &
	deg_enrichment$padj < fdr

sum(sig_idx)
which(sig_idx)

sig_idx = 
	p.adjust(pmat[, colnames(pmat) == "DUKE"]) < fdr |
	deg_enrichment$padj < fdr

sum(sig_idx)



pdf("pheno/plots/signatures/DEG_DUKE_enrichment_venn.pdf")
membership = data.frame(
	DEG=deg_enrichment$padj < fdr,
	DUKE=p.adjust(pmat[, colnames(pmat) == "DUKE"]) < fdr
	)

counts = vennCounts(membership)
vennDiagram(counts, circle.col=brewer.pal(9, "Set1"))
dev.off()

# OLD 
# -------------------------

i = 4
colnames(pmat)[i]

min_pvalue = 1e-16

x = pmat[, 1]
y = deg_enrichment$pval


# y =deg_enrichment$pval
# x = mod_tab$CAD_pval

x[x < min_pvalue] = min_pvalue
y[y < min_pvalue] = min_pvalue

plot(-log10(x), -log10(y))
cor.test(-log10(x), -log10(y))


head(sort(fits$syntax_score$p))
head(fits$syntax_score$p)

fit = list()
for (t in tissues) {
	idx = meta_row$tissue == t

	# Only complete observations
	# obs_idx = complete.cases(pheno_match$syntax) & complete.cases(tmat[, idx])

	fit[[t]] = corAndPvalue(tmat[, idx], pheno_match[[feature]])
}

lapply(fit, function(x) {
	sum(p.adjust(x$p, method="BH") < 0.1)
})


# # fit[[t]] = fastLm(pheno_match$syntax~., data=tmat[, idx][, 1:1000])
# # fit[[t]] = lm(pheno_match$syntax~., data=tmat[, idx][, 1:100])
# fit[[t]] = lm(pheno_match$syntax~., data=tmat[, idx][, 1:1000])

# fit[[t]] = glmnet(tmat[obs_idx, idx][, 1:10], pheno_match$syntax[obs_idx])

# fit[[t]] = glmnet(data.matrix(tmat[obs_idx, idx][, 1:10]), pheno_match$syntax[obs_idx])

# fit[[t]] = corAndPvalue(tmat[obs_idx, idx][, 1:10], pheno_match$syntax[obs_idx])

# fit[[t]] = corAndPvalue(tmat[, idx][, 1:1000], pheno_match$syntax)
# fit[[t]] = fastLm(pheno_match$syntax~., data=tmat[, idx])


sum(p.adjust(fit[[t]]$p) < 0.1)

apply(tmat, 1, function(row) sum(is.na(row)))


summary(fit[[t]])
