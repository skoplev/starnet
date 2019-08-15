rm(list=ls())

library(data.table)
library(RColorBrewer)

library(qvalue)
library(gplots)

library(igraph)
library(rcausal)
library(edgeR)


data_dir = "~/DataProjects/cross-tissue"

setwd("~/GoogleDrive/projects/STARNET/cross-tissue")
source("src/base.R")
source("src/parse.R")
source("src/models/enrichment.R")


# Load module table
mod_tab = fread("co-expression/tables/module_tab.csv")

# Combined
load(file.path(data_dir, "STARNET/gene_exp_norm_reshape/expr_recast.RData"),
	verbose=TRUE)


expr = parseExprTable(expr_recast)

# Cleanup
rm(expr_recast)



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


fdr = 0.05
base_mean_lim = 10
# base_mean_lim = 50
fc_lim = log2(1.3)  # +- 30%


idx_sig = deg_all$padj < fdr &
	deg_all$baseMean > base_mean_lim &
	abs(deg_all$log2FoldChange) > fc_lim

# sum(idx_sig)

deg_all_sig = deg_all[idx_sig, ]


# sort(table(deg_all_sig$gene_biotype))


# groups = c("antisense", "lincRNA", "sense_intronic", "snoRNA", "snRNA")
# groups = c("antisense", "lincRNA", "sense_intronic", "miRNA", "snoRNA", "snRNA")
groups = c("antisense", "lincRNA", "sense_intronic", "processed_transcript", "miRNA", "snoRNA", "snRNA")

deg_enrichment = lapply(groups, function(biotype) {
	# message(biotype)
	deg_enrichment = hyperGeometricModuleTestCore(
		clust_genes=mod$tissue_ensembl,
		clust=mod$clust,
		genes=deg_all_sig$tissue_ensembl[deg_all_sig$gene_biotype == biotype])

	# deg_enrichment$padj = p.adjust(deg_enrichment$pval)
	deg_enrichment$biotype = biotype

	deg_enrichment$module = 1:nrow(deg_enrichment)

	deg_enrichment = deg_enrichment[order(deg_enrichment$pval), ]

	return(deg_enrichment)
})
names(deg_enrichment) = groups

deg_enrichment = Reduce(rbind, deg_enrichment)
deg_enrichment$padj = p.adjust(deg_enrichment$pval)

fdr = 0.05

deg_enrichment_sig = deg_enrichment[deg_enrichment$padj < fdr, ]
deg_enrichment_sig


min_pval = 1e-16
deg_enrichment_sig$padj[deg_enrichment_sig$padj < min_pval] = min_pval

pdf("case-control/plots/noncoding_enrichment.pdf", height=3.0, width=10)
colors = brewer.pal(8, "Set2")

barplot(-log10(deg_enrichment_sig$padj),
	col=colors[as.integer(factor(deg_enrichment_sig$biotype))],
	names.arg=deg_enrichment_sig$module,
	las=2,
	xlab="Co-expression module",
	ylab=expression("-log"[10] * " p (BH)"),
	main=paste(groups, collapse=", ")
)
dev.off()


deg_enrichment_sig


# Write table of top-5 GO enrichment for modules
write.csv(mod_tab[deg_enrichment_sig$module, c("V1", "top_go_bp")], "case-control/noncodingModules/noncoding_modules_GO.csv", row.names=FALSE)

# deg_enrichment_sig[deg_enrichment_sig$module == 167, ]


k = 109  # module
k = 159  # module
k = 167  # module

data.frame(mod[mod$clust == k, ])

write(mod$gene_symbol[mod$clust == k],
	file=paste0("case-control/signatures/modules/mod_", k, ".txt")
)


k = 159  # module

k = 167  # module
# tissue = "VAF"
# tissue = "SKLM"
tissue = "AOR"

write(mod$gene_symbol[mod$clust == k & mod$tissue == tissue],
	file=paste0("case-control/signatures/modules/mod_", k, "_", tissue, ".txt")
)


# Write differential expression data for modules
# --------------------------------------------------

k = 167

# mod_tab

ensembl_ids = mod$ensembl[mod$clust == k & mod$tissue == "AOR"]

deg_aor_mod = deg$AOR[deg$AOR$V1 %in% ensembl_ids]

colnames(deg_aor_mod)[1] = "ensembl"

write.csv(deg_aor_mod, "case-control/signatures/modules/mod_167_AOR_DEG.csv", row.names=FALSE)



# Co-expression modules for differentially expressed genes
# -------------------------------------------------------

# Median gene expression in log2(n + 1), adjusted
median_expression = apply(expr$mat, 1, median, na.rm=TRUE)


k = 167  # module

# idx = mod$clust == k & mod$tissue_ensembl %in% deg_all_sig$tissue_ensembl

# Filter for absolute gene expression 
idx = mod$clust == k & mod$tissue_ensembl %in% deg_all_sig$tissue_ensembl & median_expression > log2(25)

# sum(idx)

# expr$meta_row$gene_biotype[idx]

sort(table(expr$meta_row$gene_biotype[idx]))


# Correlation matrix
cmat = cor(t(expr$mat[idx, ]), use="pairwise.complete.obs")

heatmap.2(cmat,
	trace="none",
	col=rev(colorRampPalette(brewer.pal(9, "RdBu"))(100)),
	cexRow=0.4,
	cexCol=0.4
)


# Bayesian network
# g = fges(t(expr$mat[idx, ]), maxDegree=20)
g = fges(t(expr$mat[idx, ]), maxDegree=30)


g = igraph.from.graphNEL(g$graphNEL)

el = as_edgelist(g)

# Subgraph including lincRNA

# lincRNA IDs
lincRNAs = expr$meta_row$tissue_transcript_id[expr$meta_row$gene_biotype == "lincRNA"]


el_lincRNA = el[el[, 1] %in% lincRNAs | el[, 2] %in% lincRNAs, ]

g_lincRNA = graph_from_edgelist(el_lincRNA)

# Annotate netork nodes
V(g_lincRNA)$tissue = sapply(strsplit(names(V(g_lincRNA)), "_"), function(x) x[1])
V(g_lincRNA)$gene = sapply(strsplit(names(V(g_lincRNA)), "_"), function(x) x[2])

V(g_lincRNA)$biotype = expr$meta_row$gene_biotype[match(names(V(g_lincRNA)), expr$meta_row$tissue_transcript_id)]

V(g_lincRNA)$degree_all = igraph::degree(g_lincRNA, mode="all")

write_graph(g_lincRNA, paste0("case-control/lincRNAnetw/mod_2_", k, ".gml"), "gml")


# Overall DEG enrichment and SYNTAX and DUKE
# ----------------------------------------------------

mod_ids = c(167, 109, 159)

# k = 109  # module
# k = 159  # module
# k = 167  # module


pheno_pmat = fread("pheno/tables/pheno_pval.csv")
pheno_pmat = as.data.frame(pheno_pmat)

pheno_pmat[mod_ids, ]

pheno_pmat_sub = pheno_pmat[mod_ids, c("case_control_DEG", "syntax_score", "DUKE")]
colnames(pheno_pmat_sub) = c("DEG", "SYNTAX", "DUKE")
rownames(pheno_pmat_sub) = mod_ids
min_pval = 1e-40

pheno_pmat_sub[pheno_pmat_sub < min_pval] = min_pval


# pdf("case-control/plots/network_enrichment_167_109_159.pdf", width=3, height=4)
pdf("case-control/plots/network_enrichment_167_109_159.pdf", width=2, height=3)
par(mar=c(6, 4, 2, 2))

# Init plot
plot(0, 0,
	type="n",
	xlim=c(1, 3),
	ylim=range(0, -log10(pheno_pmat_sub)),
	xlab="", ylab="",
	axes=FALSE
)

axis(1, at=1:3, labels=colnames(pheno_pmat_sub), las=2)
axis(2, pos=0.6, las=1)

mtext("Enrichment (-log10 p)", side=2, line=3, cex=1.1)


abline(h=-log10(0.05), lty=2,
	col="black",
	xpd=FALSE)

pts_pch = c(21, 23, 24)
pts_cols = brewer.pal(9, "Set1")

for (i in 1:length(mod_ids)) {
	x = 1:ncol(pheno_pmat_sub)
	y = -log10(pheno_pmat_sub[i, ])
	lines(x, y, col=pts_cols[i])
	points(x, y,
		pch=pts_pch[i],
		bg="white",
		col=pts_cols[i],
		lwd=1.3,
		xpd=TRUE
	)
}

legend("topright", legend=mod_ids, pch=pts_pch, col=pts_cols, pt.lwd=1.3)
dev.off()



# Eigengene correlations
# ----------------------------------------------------------

pheno = fread("~/GoogleDrive/projects/STARNET/phenotype/data/current/STARNET_main_phenotype_table.2017_12_03.tsv")

eig = fread("co-expression/tables/eigengene_mat.csv", header=TRUE)

eig_mat = eig[, -1]
eig_mat = data.matrix(eig_mat)
rownames(eig_mat) = eig$V1  # patient ID

pheno_matched = pheno[match(rownames(eig_mat), pheno$starnet.ID), ]


k = 167


pdf("case-control/plots/mod_167_eigengene_cor.pdf", height=3.5, width=6)

feats = c("fP-LDL-Chol(mmol/l)", "fP-HDL-Chol(mmol/l)")
# feats = c("SBP", "DBP")
par(mfrow=c(1, 2))

for (feat in feats) {

	if (feat %in% c("fP-LDL-Chol(mmol/l)", "fP-HDL-Chol(mmol/l)")) {
		outlier = pheno_matched[[feat]] > 30
	} else {
		outlier = rep(FALSE, nrow(pheno_matched))
	}
	message("outliers: ", sum(outlier, na.rm=TRUE))

	x = pheno_matched[[feat]]
	x = x[!outlier]

	y = eig_mat[, k]
	y = y[!outlier]


	fit = lm(y~x)

	cor_test = cor.test(x, y)
	plot(x, y,
		cex=0.4,
		pch=21,
		lwd=0.5,
		col="black",
		bg=rgb(150, 150, 150, maxColorValue=255),
		xlab=feat, ylab=paste("Eigengene", k)
	)

	abline(fit$coefficients,
		# col="red",
		col=brewer.pal(9, "Set1")[1],
		lwd=1.5
	)
	legend("topright",
		legend=c(
			paste0("P=", format(cor_test$p.value, digits=3)),
			paste0("r=", format(cor_test$estimate, digits=3))
		),
		bty="n",
		cex=0.8
	)
}

dev.off()






# Gender associations despite adjustment?


# Comparative boxplot
dotBoxplot = function(val, col, ...) {
	boxplot(val,
		ylim=range(val$male, val$female, na.rm=TRUE),
		frame=FALSE,
		outline=FALSE,
		...
	)

	jit_fac = 0.25
	points(
		jitter(rep(1, length(val[[1]])), amount=jit_fac),
		col=colors[1],
		cex=0.4,
		val$male, pch=16)

	points(
		jitter(rep(2, length(val[[2]])), amount=jit_fac),
		col=colors[2],
		cex=0.4,
		val$female, pch=16)
}


colors = brewer.pal(9, "Set1")


k = 167
gender = pheno_matched$Sex
gender[gender == "patient nr male7female"] = NA

contingency_tab = table(gender, eig_mat[, k] > 0.08)
# contingency_tab = table(gender, abs(eig_mat[, k]) > 0.01)
chisq.test(contingency_tab)



table(pheno_matched$Sex)

# # FTX expression 
# stopifnot(all(colnames(expr$mat) == pheno_matched$starnet.ID))

# val = list()
# j = which(expr$meta_row$gene_symbol == "FTX" & expr$meta_row$tissue == "AOR")
# val$male = expr$mat[j, pheno_matched$Sex == "male"]
# val$female = expr$mat[j, pheno_matched$Sex == "female"]

# j = which(expr$meta_row$gene_symbol == "FTX" & expr$meta_row$tissue == "AOR")
# i = which(expr$meta_row$gene_symbol == "XIST" & expr$meta_row$tissue == "AOR")

# plot(expr$mat[i, ], expr$mat[j, ])
# cor.test(expr$mat[i, ], expr$mat[j, ])



#
# aor_mat = fread("~/DataProjects/STARNET/oscar_mat/STARNET.AOR.exp.mat.F.gene_filt.EDAseq.gc.lm.or.RQN")
aor_mat = fread("~/DataProjects/STARNET/oscar_mat/STARNET.AOR.exp.mat")
aor_mat = as.data.frame(aor_mat)
rownames(aor_mat) = as.character(aor_mat$id)
aor_mat = aor_mat[, -1]
colnames(aor_mat) = sapply(strsplit(colnames(aor_mat), "_"), function(x) x[2])

aor_mat = data.matrix(aor_mat)

aor_mat = cpm(aor_mat)  # CPM normalization



# Match phenotype data
pheno_matched_aor = pheno[match(colnames(aor_mat), pheno$starnet.ID), ]

i = grep("ENSG00000230590", rownames(aor_mat))  # FTX row

ftx_val = list()
ftx_val$male = aor_mat[i, pheno_matched_aor$Sex == "male"]
ftx_val$female = aor_mat[i, pheno_matched_aor$Sex == "female"]

# dotBoxplot(ftx_val, col=colors, ylab="FTX (CPM)")

t.test(ftx_val$male, ftx_val$female)
wilcox.test(ftx_val$male, ftx_val$female)


i = grep("ENSG00000229807", rownames(aor_mat))  # XIST row

xist_val = list()
xist_val$male = aor_mat[i, pheno_matched_aor$Sex == "male"]
xist_val$female = aor_mat[i, pheno_matched_aor$Sex == "female"]




# Eigengene boxplot by gender
k = 167
val = list()
val$male = eig_mat[pheno_matched$Sex == "male", k]
val$female = eig_mat[pheno_matched$Sex == "female", k]

# pdf("case-control/plots/mod167/mod167_gender_boxplot.pdf", height=3.5, width=2.5)
pdf("case-control/plots/mod167/mod167_gender_boxplot.pdf", height=3.5, width=5.0)

par(mfrow=c(1, 3))
dotBoxplot(val, col=colors, ylab=paste("Eigengene", k))
dotBoxplot(ftx_val, col=colors, ylab="FTX (CPM)")
dotBoxplot(xist_val, col=colors, ylab="XIST (CPM)")

t.test(val$male, val$female)
wilcox.test(val$male, val$female)

dev.off()


# FTX correlations to phenotypes
# ------------------------------------

tissue = "AOR"
# tissue = "MAM"

pheno_matched = pheno[match(colnames(expr$mat), pheno$starnet.ID), ]



feats = c("fP-LDL-Chol(mmol/l)", "fP-HDL-Chol(mmol/l)")
# feats = c("syntax_score", "DUKE", "fP-LDL-Chol(mmol/l)", "fP-HDL-Chol(mmol/l)", "SBP", "DBP", "BMI(kg/m2)", "CRP(mg/l)", "HbA1c(%)", "P-Chol(mmol/l)", "fP-TG(mmol/l)")

# feats = c("SBP", "DBP")
# par(mfrow=c(3, 4))

pdf("case-control/plots/FTX_LDL_HDL_cor.pdf", height=6.5, width=6)

feats = c("syntax_score", "DUKE")
pdf("case-control/plots/FTX_SYNTAX_DUKE_cor.pdf", height=6.5, width=6)

par(mfrow=c(2, 2))
for (tissue in c("MAM", "AOR")) {
	i = which(expr$meta_row$tissue == tissue & expr$meta_row$gene_symbol == "FTX")

	for (feat in feats) {

		if (feat %in% c("fP-LDL-Chol(mmol/l)", "fP-HDL-Chol(mmol/l)")) {
			outlier = pheno_matched[[feat]] > 30
		} else {
			outlier = rep(FALSE, nrow(pheno_matched))
		}
		message("outliers: ", sum(outlier, na.rm=TRUE))

		x = pheno_matched[[feat]]
		x = x[!outlier]

		y = expr$mat[i, ]
		y = y[!outlier]


		fit = lm(y~x)

		cor_test = cor.test(x, y)
		plot(x, y,
			cex=0.4,
			pch=21,
			lwd=0.5,
			col="black",
			bg=rgb(150, 150, 150, maxColorValue=255),
			xlab=feat, ylab=paste0(expr$meta_row$gene_symbol[i], " expression ", tissue, ", adj. log(n + 1)")
		)

		abline(fit$coefficients,
			# col="red",
			col=brewer.pal(9, "Set1")[1],
			lwd=1.5
		)
		legend("topright",
			legend=c(
				paste0("P=", format(cor_test$p.value, digits=3)),
				paste0("r=", format(cor_test$estimate, digits=3))
			),
			bty="n",
			cex=0.8
		)
	}
}
dev.off()


# Key driver plots
# --------------------------------

# mod_tab
# expr$meta_row

kd = fread("co-expression/tables/modules.results.txt")


kd$biotype = expr$meta_row$gene_biotype[match(kd$NODE, as.character(expr$meta_row$tissue_transcript_id))]

kd$gene = sapply(strsplit(as.character(kd$NODE), "_"), function(x) x[2])
kd$tissue = sapply(strsplit(as.character(kd$NODE), "_"), function(x) x[1])
kd$ensembl = sapply(strsplit(as.character(kd$NODE), "_"), function(x) x[3])
kd$ensembl_base = sapply(strsplit(as.character(kd$ensembl), "[.]"), function(x) x[1])


kd_sub = list()
for (k in c(167, 159, 109)) {
	kd_sub[[k]] = kd[kd$MODULE == k, ]
}

pdf("case-control/plots/key_drivers_167_109_159.pdf", width=3.0)
colors = brewer.pal(9, "Set1")

m = 50

# Get key driver subset
x = kd_sub[[167]][1:m, ]


# Color key drivers based on biotype
pts_col = rep("white", nrow(x))
pts_col[x$biotype == "lincRNA"] = colors[1]
pts_col[x$biotype == "protein_coding"] = colors[2]

plot(-log10(x$P), 1:nrow(x),
	xlim=range(0, -log10(x$P)),
	ylim=rev(range(1, m)),  # reverse plotting order
	pch=21,
	bg=pts_col,
	bty="l",
	ylab="Ranking",
	xlab="Key driver (-log10 p)",
	main=paste(c("109", "167", "[159]"), collapse=",")
)

text(-log10(x$P), 1:nrow(x), labels=x$gene,
	pos=4,
	cex=0.7,
	xpd=TRUE
)

# Module 109
# x = -log10(kd_sub[[109]]$P[1:m])
x = kd_sub[[109]][1:m]

pts_col = rep("white", nrow(x))
pts_col[x$biotype == "lincRNA"] = colors[1]
pts_col[x$biotype == "protein_coding"] = colors[2]

points(-log10(x$P), 1:nrow(x),
	pch=21,
	bg=pts_col
)

text(-log10(x$P), 1:nrow(x),
	labels=x$gene,
	pos=4,
	cex=0.7,
	xpd=TRUE
)


# No key drivers for module 159...
x = -log10(kd_sub[[159]]$P[1:m])
x = rev(x)
points(x, 1:length(x))


legend("bottomright",
	legend=c("lincRNA", "Protein coding"), pch=21, pt.bg=colors,
	cex=0.8)
dev.off()


# Write list of key driver ensembl IDs of module 167
# For identifying eQTL
# -----------------------------------------------------
k = 167
kd_sub = kd[kd$MODULE == k, ]


table(kd_sub$tissue)
table(kd_sub$biotype)

kd[kd$MODULE == k, ]

write(kd_sub$ensembl_base, "case-control/signatures/modules/key_drivers/mod_167_kd.txt")

aor_eqtl = fread("case-control/signatures/modules/key_drivers/mod_167_kd_eQTL_AOR_format.csv")


# kd_sub$biotype 

aor_eqtl[aor_eqtl$gene_biotype == "lincRNA" & aor_eqtl$adj.p < 0.05, ]
unique(aor_eqtl$hgnc_symbol[aor_eqtl$gene_biotype == "lincRNA" & aor_eqtl$adj.p < 0.05])
unique(aor_eqtl$gene[aor_eqtl$gene_biotype == "lincRNA" & aor_eqtl$adj.p < 0.05])


aor_eqtl[aor_eqtl$hgnc_symbol == "FTX", ]
aor_eqtl[aor_eqtl$hgnc_symbol == "CALU", ]
aor_eqtl[aor_eqtl$hgnc_symbol == "MON2", ]
aor_eqtl[aor_eqtl$hgnc_symbol == "HSP90B1", ]
aor_eqtl[aor_eqtl$hgnc_symbol == "LIMA1", ]
aor_eqtl[aor_eqtl$hgnc_symbol == "RSRP1", ]
aor_eqtl[aor_eqtl$hgnc_symbol == "PDIA6", ]
aor_eqtl[aor_eqtl$hgnc_symbol == "CDKN1B", ]
aor_eqtl[aor_eqtl$hgnc_symbol == "LYST", ]
aor_eqtl[aor_eqtl$hgnc_symbol == "RUNX1T1", ]

kd[grep("ENSG00000273329.1", kd$NODE), ]


kd[kd$MODULE == k & kd$gene == "TAS2R14", ]



# Enrichment results for module 167 AOR, barplot
# ----------------------------------------

enrich_167_AOR = fread("case-control/signatures/modules/enrich/PANTHER_mod_167_AOR.txt", skip=10)

min_fold_enrichment = 4
m = 15  # top results

enrich_167_AOR = enrich_167_AOR[enrich_167_AOR[["upload_1 (fold Enrichment)"]] > min_fold_enrichment, ]

# Order by FDR
enrich_167_AOR = enrich_167_AOR[order(enrich_167_AOR[["upload_1 (FDR)"]]), ]

pdf("case-control/plots/enrichment_PANTHER_mod_167_AOR.pdf", width=6, height=5)
p = enrich_167_AOR[["upload_1 (FDR)"]]
names(p) = sapply(strsplit(enrich_167_AOR[["GO biological process complete"]], " \\("), function(x) x[1])
p = p[1:m]
p = rev(p)

par(mar=c(4, 24, 2, 2))
# barplot(-log10(p[1:m]), horiz=TRUE, las=2)
barplot(-log10(p), horiz=TRUE, las=2,
	col=rgb(220, 220, 220, maxColorValue=255),
	border=rgb(150, 150, 150, maxColorValue=255),
	cex.names=0.8,
	xlab="-log10 p (BH)"
)
dev.off()




