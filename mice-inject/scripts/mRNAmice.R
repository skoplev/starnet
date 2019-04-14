rm(list=ls())

library(data.table)
library(RColorBrewer)
library(DESeq2)
library(biomaRt)
library(limma)  # for vennDiagram
library(gplots)
library(UpSetR)


setwd("~/GoogleDrive/projects/STARNET/cross-tissue")

source("src/models/cor.R")

# counts = fread("mice-inject/data/feature_counts/gene_counts.txt")
counts = fread("mice-inject/data/feature_counts_all/gene_counts.txt")

# Format into count matrix
count_mat = counts[, -1:-6]
rownames(count_mat) = counts$Geneid

colnames(count_mat) = gsub(".fastq.*", "", colnames(count_mat))  # shorten to sample IDs

count_mat = data.matrix(count_mat)

# sample metadata
samples = data.frame(
	id=colnames(count_mat),
	day=sapply(strsplit(colnames(count_mat), "_|[.]"), function(x) x[2]),
	condition=factor(sapply(strsplit(colnames(count_mat), "_"), function(x) x[1]))
)

# Map gene names
mart = useMart("ensembl", dataset="mmusculus_gene_ensembl")

gene_map = getBM(attributes=c("mgi_symbol", "ensembl_gene_id", "gene_biotype"), mart=mart)

getGeneSymbol = function(ensembl) {
	symbol = gene_map$mgi_symbol[match(ensembl, gene_map$ensembl_gene_id)]
	return(symbol)
}


# Load MGI human-mouse homology data
# ---------------------------------------------------
human_mouse_homology = fread("~/DataProjects/cross-tissue/MGI/HOM_MouseHumanSequence.rpt")

# Separate table into human and mouse
human_homology = human_mouse_homology[
	human_mouse_homology[["Common Organism Name"]] == "human", ]

mouse_homology = human_mouse_homology[
	human_mouse_homology[["Common Organism Name"]] == "mouse, laboratory", ]

mouse2human = function(mouse_symbols) {
	# Get homologene IDs
	homologene_ids = mouse_homology[["HomoloGene ID"]][
		match(mouse_symbols, mouse_homology$Symbol)
	]

	human_symbols = human_homology$Symbol[
		match(homologene_ids, human_homology[["HomoloGene ID"]])
	]

	# out = human_symbols

	# Default to input if no human homologue is found
	# missing = is.na(out) | out == ""
	# out[missing] = mouse_symbols[missing]

	return(human_symbols)
}

# Fit DESeq 2 model
# --------------------------------------
include = apply(count_mat, 1, max) > 5  # filter transcripts based on minimum counts
sum(include)

dds = DESeqDataSetFromMatrix(countData=count_mat[include, ],
	colData=samples,
	# design=~condition
	design=~day + condition
)

dds = DESeq(dds)

# endocrines = c("EPDR1", "FCN2", "LBP")
endocrines = c("EPDR1", "FCN2", "LBP", "FSTL3")
results = lapply(endocrines, function(pert) {
	res = results(dds, contrast=c("condition", pert, "Veh"))

	res$mgi_symbol = getGeneSymbol(rownames(res))
	res$hgnc_symbol = mouse2human(res$mgi_symbol)  # human gene symbols

	res = res[order(res$pvalue), ]
	return(res)
})
names(results) = endocrines


# Write log2 FC rankings
for (i in 1:length(results)) {
	res = results[[i]]

	res = res[order(res$log2FoldChange, decreasing=TRUE), ]

	# Remove emtpy gene symbols
	res = res[!is.na(res$hgnc_symbol) & res$hgnc_symbol != "", ]

	write.table(res[, c("hgnc_symbol", "log2FoldChange")],
		paste0("mice-inject/DEG_ranking_all/log2FC_", names(results)[i], ".rnk"),
		quote=FALSE,
		row.names=FALSE,
		sep="\t"
	)
}

head(data.frame(results[[4]]), 20)
head(data.frame(results[[4]][results[[4]]$log2FoldChange > 0, ]), 20)



# Get gene lists of differentially expressed genes
fdr = 0.05
min_log2fc = 1

sig_genes_human = lapply(results, function(res) {
	# sig = res$padj < fdr
	sig = res$padj < fdr & abs(res$log2FoldChange) > min_log2fc
	sig_symbols = res$hgnc_symbol[sig]  # homologous human gene symbols

	sig_symbols = na.omit(sig_symbols)
	sig_symbols = unique(sig_symbols)
	sig_symbols = sig_symbols[sig_symbols != ""]
	return(sig_symbols)
})

sig_genes_mouse = lapply(results, function(res) {
	# sig = res$padj < fdr 
	sig = res$padj < fdr & abs(res$log2FoldChange) > min_log2fc
	sig_symbols = res$mgi_symbol[sig]  # homologous human gene symbols

	sig_symbols = na.omit(sig_symbols)
	sig_symbols = unique(sig_symbols)
	sig_symbols = sig_symbols[sig_symbols != ""]
	return(sig_symbols)
})

# Some counts
sapply(sig_genes_human, length)
sapply(sig_genes_mouse, length)

write(
	paste(intersect(sig_genes_mouse[[2]], sig_genes_mouse[[3]]), collapse=", "),
	"mice-inject/DEG_all/mouse_FCN2_LBP.txt")

# Write gene lists
for (cond in names(sig_genes_mouse)) {
	write(sig_genes_mouse[[cond]], paste0("mice-inject/DEG_all/mouse_", cond, ".txt"))
	write(sig_genes_human[[cond]], paste0("mice-inject/DEG_all/human_", cond, ".txt"))
}

# Print specific resutls
gene = "DHCR7"
gene = "FAM213A"
gene = "RDH11"
gene = "DHCR24"
gene = "FAM213A"

gene = "PCSK9"
for (i in 1:4) {
	print(results[[i]][which(results[[i]]$hgnc_symbol == gene), ])
}


# Venn counts comparing LBP and FCN2 signatures
# -------------------------------------
# Create indicator matrix
gene_sets = sig_genes_mouse[c("LBP", "FCN2")]
ind_mat = sapply(gene_sets, function(genes) {
	table(factor(genes, levels=unique(unlist(gene_sets))))
})

pdf("mice-inject/plots/venn_LBP_FCN2.pdf")
vennDiagram(ind_mat, circle.col=brewer.pal(9, "Set1")[3:2])
dev.off()



# Pathway enrichment, from Enrichr
# -------------------------------------
kegg_fcn2 = fread("mice-inject/DEG_all/enrichr/KEGG/FCN2_KEGG_2016_table.txt")
kegg_lbp = fread("mice-inject/DEG_all/enrichr/KEGG/LBP_KEGG_2016_table.txt")
kegg_fstl3 = fread("mice-inject/DEG_all/enrichr/KEGG/FSTL3_KEGG_2016_table.txt")


path_fdr = 0.05
terms = c(
	kegg_fcn2$Term[kegg_fcn2[["Adjusted P-value"]] < path_fdr],
	kegg_lbp$Term[kegg_lbp[["Adjusted P-value"]] < path_fdr],
	kegg_fstl3$Term[kegg_fstl3[["Adjusted P-value"]] < path_fdr]
)
terms = unique(terms)

default = function(value) {
	if (length(value) == 1) {
		return(value)
	} else if (length(value) > 1) {
		return(value[1])
	} else {
		return(NA)
	}
}

kegg_pvals = sapply(terms, function(x) {
	fcn2_pval = kegg_fcn2[["Adjusted P-value"]][which(kegg_fcn2$Term == x)]
	lbp_pval = kegg_lbp[["Adjusted P-value"]][which(kegg_lbp$Term == x)]
	fstl3_pval = kegg_fstl3[["Adjusted P-value"]][which(kegg_fstl3$Term == x)]

	pvals = c(default(fcn2_pval), default(lbp_pval), default(fstl3_pval))
	names(pvals) = c("Fcn2", "Lbp", "Fslt3")
	return(pvals)
})

colnames(kegg_pvals) = sapply(strsplit(colnames(kegg_pvals), "_"), function(x) x[1])


pdf("mice-inject/plots/KEGG_pathway_barplot.pdf", width=5, height=5)
idx = order(apply(kegg_pvals, 2, min, na.rm=TRUE), decreasing=TRUE)

# idx = order(apply(-log10(kegg_pvals), 2, sum, na.rm=TRUE), decreasing=TRUE)
kegg_pvals = kegg_pvals[, idx]

par(mar=c(4, 16, 3, 3))
barplot(-log10(kegg_pvals),
	beside=TRUE,
	horiz=TRUE,
	las=2,
	# col=brewer.pal(9, "Set1")[c(3, 2, 4)],
	col=brewer.pal(9, "Set1")[c(2, 3, 4)],
	legend.text=TRUE,
	ylab="KEGG pathway",
	xlab="-log10 p (BH)"
)

abline(v=0)

abline(v=-log10(0.1), col="grey", lty=2)
dev.off()


# kegg_fcn2 = fread("mice-inject/DEG/enrichr/mouse_FCN2_KEGG_2016_table.txt")
# kegg_lbp = fread("mice-inject/DEG/enrichr/mouse_LBP_KEGG_2016_table.txt")




# terms = c(
# 	"Non-alcoholic fatty liver disease (NAFLD)_Homo sapiens_hsa04932",
# 	"Protein processing in endoplasmic reticulum_Homo sapiens_hsa04141"
# )

# kegg_pvals = sapply(terms, function(x) {
# 	fcn2_pval = kegg_fcn2[["Adjusted P-value"]][kegg_fcn2$Term == x]
# 	lbp_pval = kegg_lbp[["Adjusted P-value"]][kegg_lbp$Term == x]

# 	pvals = c(lbp_pval, fcn2_pval)
# 	names(pvals) = c("Lbp", "Fcn2")
# 	return(pvals)
# })


# pdf("mice-inject/plots/KEGG_enrichrment_LBP_FCN2.pdf", width=2, height=5.0)
# par(mar=c(16, 4, 2, 2))
# barplot(-log10(kegg_pvals),
# 	beside=TRUE,
# 	las=2,
# 	ylab="-log10 p (BH)",
# 	legend.text=TRUE,
# 	col=brewer.pal(9, "Set1")[3:2])

# abline(h=-log10(0.05), lty=2, col="grey")
# abline(h=0)
# dev.off()


# Normalization based on size factors
# ---------------------------------------
# sizeFactors(dds)

library(edgeR)

count_norm = sweep(count_mat[include, ], 2, sizeFactors(dds), "/")

rownames(count_norm) = getGeneSymbol(rownames(count_norm))

# cpm_norm = cpm(count_norm)

cpm_norm = cpm(count_mat[include, ])
rownames(cpm_norm) = getGeneSymbol(rownames(cpm_norm))


# Module 98 targeted analysis
# ------------------------------------------
library(metap)

mod_tab = fread("co-expression/tables/modules.csv")

genes98 = mod_tab$gene_symbol[mod_tab$clust == 98]

results_mod98 = lapply(results, function(res) {
	res[res$hgnc_symbol %in% genes98, ]
})


lapply(results_mod98, as.data.frame)

lapply(results_mod98, function(res) {
	# sumz(res$pvalue)  # Stouffer metap
	sumlog(res$pvalue)  # Fisher's
})


par(mfrow=c(4, 1))
for (i in 1:4) {
	n = nrow(results_mod98[[i]])
	pvals = results_mod98[[i]]$pvalue
	plot(
		-log10(pvals),
		pch=16,
		bty="n"
	)
	abline(h=-log10(0.05), col="grey", lty=2)

	idx = pvals < 0.05
	text(-log10(pvals[idx]),
		xpd=TRUE,
		labels=results_mod98[[i]]$mgi_symbol[idx],
		pos=4
	)
}




# UpSet plot of overlaps with DEGs
pdf("mice-inject/plots/DEG_upset.pdf", height=3.0, width=4)
# pdf("mice-inject/plots/DEG_upset.pdf")
gene_sets = c(sig_genes_human, list("Module 98"=genes98))
upset(fromList(gene_sets),
	sets.bar.color=brewer.pal(9, "Set1")[1],
	mainbar.y.label="Gene overlap",
	sets.x.label="DEGs (FDR < 5%)",
	mb.ratio=c(0.5, 0.5)
	# empty.intersections="on"
)
dev.off()

sort(table(unlist(gene_sets)))


# Plots of differentially expressed genes
# -----------------------------------------

plotExpr = function(gene, mat, conditions=c("Veh", "EPDR1", "FCN2", "LBP", "FSTL3")) {
	values = lapply(conditions, function(cond) {
		i = match(gene, rownames(mat))
		expr_values = mat[i, samples$condition == cond]
		return(expr_values)
	})
	names(values) = conditions

	df = data.frame(
		conditions=factor(samples$condition, levels=conditions),
		gene=mat[match(gene, rownames(mat)), ]
	)

	colors = brewer.pal(9, "Set1")[c(9, 1:(length(conditions) - 1))]

	plotmeans(gene~conditions, data=df,
		ylab=paste(gene, "(CPM)"),
		xlab="",
		frame=FALSE,
		connect=FALSE,
		minbar=0,
		pch="-",
		cex=1.5,
		# pch=7,
		# barcol="black",
		barcol=colors,
		col=colors,
		las=2,
		n.label=FALSE
	)

	for (i in 1:length(values)) {
		points(
			jitter(rep(i, length(values[[i]])), amount=0.2),
			values[[i]],
			cex=1,
			# pch=16,
			# col=colors[i])
			col="black",
			pch=21,
			bg=colors[i])
	}

	# Print results table
	for (i in 1:(length(results))) {
		print(results[[i]][which(results[[i]]$mgi_symbol == gene), ])
	}

}


mod98_overlap = lapply(sig_genes_human, function(genes) {
	intersect(genes, genes98)
})
mod98_overlap = unlist(mod98_overlap)


# pdf("mice-inject/plots/gene_expr_DEG_mod98.pdf", height=2.1, width=5.5)
# pdf("mice-inject/plots/gene_expr_DEG_mod98.pdf", height=2.1, width=6.0)
pdf("mice-inject/plots/gene_expr_DEG_mod98.pdf", height=2.1, width=6.0)
# pdf("mice-inject/plots/gene_expr_DEG_mod98.pdf", height=2.1, width=7.0)
par(mfrow=c(1, 4))
plotExpr("Dhcr7", cpm_norm)
plotExpr("Dhcr24", cpm_norm)
plotExpr("Lss", cpm_norm)
# plotExpr("Fam213a", cpm_norm)
plotExpr("Hmgcs1", cpm_norm)
dev.off()


par(mfrow=c(1, 4))
plotExpr("Dhcr7", count_norm)
plotExpr("Fam213a", count_norm)
plotExpr("Rdh11", count_norm)
plotExpr("Dhcr24", count_norm)



plotExpr("Rxra", cpm_norm)
plotExpr("Spink1", cpm_norm)
plotExpr("Gpx3", cpm_norm)
plotExpr("Cela1", cpm_norm)


plotExpr("Pcsk9", cpm_norm)

plotExpr("Apoa4", cpm_norm)
plotExpr("Adrb1", cpm_norm)

plotExpr("Bco1", cpm_norm)

plotExpr("Ggt1", cpm_norm)

plotExpr("Cndp1", cpm_norm)

plotExpr("Abcc3", cpm_norm)
plotExpr("Derl3", cpm_norm)

plotExpr("Etnppl", cpm_norm)

plotExpr("Acat3", cpm_norm)
plotExpr("Tnnt2", cpm_norm)
plotExpr("Gja1", cpm_norm)

plotExpr("Tpm1", cpm_norm)

plotExpr("Hmgcs1", cpm_norm)

plotExpr("Taf15", cpm_norm)

plotExpr("Npc1", cpm_norm)

plotExpr("Hmgcs1", cpm_norm)

plotExpr("Pcsk9", cpm_norm)
plotExpr("Hmgcr", cpm_norm)

plotExpr("Ldlr", cpm_norm)


# Principal component analysis (PCA)
# ----------------------------------------

# pca = prcomp(t(log(count_mat + 1)), scale=TRUE)
pca = prcomp(t(log(count_norm + 1)), scale=TRUE)

colors = c(brewer.pal(9, "Set1")[1:3], "grey")

pt_col = colors[as.integer(samples$condition)]

# pt_col = colors[as.integer(samples$day)]

i = 1
j = 2
plot(pca$x[, i], pca$x[, j], col=pt_col, pch=16)
text(pca$x[, i], pca$x[, j], labels=colnames(count_mat),
	cex=0.8,
	pos=1,
	col="grey"
)


# Correlation matrices
# ------------------------------------------

library(pwr)
pwr.r.test(n=4, sig.level=0.1, power=0.8)

genes98 = mod_tab$gene_symbol[mod_tab$clust == 98]

# sum(genes98 %in% toupper(rownames(count_norm)))
# sum(genes98 %in% mouse2human(rownames(count_norm)))


idx = which(mouse2human(rownames(count_norm)) %in% genes98)

cmats = lapply(levels(samples$condition), function(cond) {
	cmat = cor(t(log10(count_norm[idx, samples$condition == cond] + 1)))
	cmat[is.na(cmat)] = 0
	return(cmat)
})
names(cmats) = levels(samples$condition)


cmats_rand = lapply(levels(samples$condition), function(cond) {
	rand_idx = sample(nrow(count_norm), size=length(idx))
	cmat = cor(t(log10(count_norm[rand_idx, samples$condition == cond] + 1)))
	cmat[is.na(cmat)] = 0
	return(cmat)
})
names(cmats_rand) = levels(samples$condition)




# i = 1
i = 2
# i = 3
# i = 4

cond = names(cmats)[i]
pdf(paste0("mice-inject/plots/", cond, "_cor_mat_heatmap.pdf"))
# pdf(paste0("mice-inject/plots/", cond, "_cor_mat_heatmap_rand.pdf"))
heatmap.2(
	cmats[[i]],
	# cmats_rand[[i]],
	trace="none",
	col=colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
	main=names(cmats)[i],
	key.title="",
	key.xlab="Pearson cor.",
	cexRow=0.7,
	cexCol=0.7
)
dev.off()


# cmat = cor(t(count_norm[idx, ]))
# cmat_rand = cor(t(
# 	count_norm[sample(nrow(count_norm), size=length(idx)), ]
# ))


# Log transform cmat
cmat = cor(t(log10(count_norm[idx, ] + 1)))

rand_idx = sample(nrow(count_norm), size=length(idx))
cmat_rand = cor(t(
	log10(count_norm[rand_idx, ] + 1)
))


# pdf("mice-inject/plots/cor_mat_heatmap.pdf")
pdf("mice-inject/plots/cor_mat_heatmap_rand.pdf")
heatmap.2(
	# cmat,
	cmat_rand,
	trace="none",
	col=colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
	main=names(cmats)[i],
	key.title="",
	key.xlab="Pearson cor.",
	cexRow=0.7,
	cexCol=0.7
)
dev.off()


# Correlation matrix permutation tests
# -----------------------------------------
library(sva)

min_expr = 4  # log2 norm counts adjusted, 
m = 10000  # number of permutations

genes98 = mod_tab$gene_symbol[mod_tab$clust == 98]
genes78 = mod_tab$gene_symbol[mod_tab$clust == 78]

# Translate mouse->human gene symbols
count_norm_human = count_norm
# rownames(count_norm_human) = mouse2human(rownames(cpm_norm))

# Pseudo-log transform
count_norm_human = log2(count_norm_human + 1)

# Batch correction
count_norm_human = ComBat(count_norm_human, samples$day)

# Filter out low expression 
count_norm_human = count_norm_human[apply(count_norm_human, 1, median) > min_expr, ]

# Get submatrices
# mat98 = count_norm_human[rownames(count_norm_human) %in% genes98, ]
# mat78 = count_norm_human[rownames(count_norm_human) %in% genes78, ]

mat98 = count_norm_human[mouse2human(rownames(count_norm_human)) %in% genes98, ]
mat78 = count_norm_human[mouse2human(rownames(count_norm_human)) %in% genes78, ]


# Construct samples indicies for running permutation tests
sample_idx = lapply(unique(samples$condition), function(cond) {
	samples$condition == cond
})
names(sample_idx) = unique(samples$condition)

sample_idx$all = !is.na(samples$condition)


# Run permutation tests per co-expression module
cor_tests98 = lapply(sample_idx, function(idx) {
	corPermuteTest(mat98[, idx], count_norm_human[, idx], m=m)
})

cor_tests78 = lapply(sample_idx, function(idx) {
	corPermuteTest(mat78[, idx], count_norm_human[, idx], m=m)
})
save(cor_tests98, cor_tests78, file="mice-inject/data/permuteTest/cor_tests.RData")


plotPermTest = function(test, ...) {
	hist(test$cor_mean_null,
		breaks=20,
		xlim=range(c(test$cor_mean_null, test$cor_mean)),
		ylab="",
		# xlab="Mean |r|",
		...)

	legend("topright",
		bty="n",
		legend=paste0("P = ", format(test$p_val, digits=4))
	)
	abline(v=test$cor_mean, lwd=1.5)
}


col_idx = c(1:4, 9, 7)
colors = brewer.pal(9, "Set1")[col_idx]
colors_light = brewer.pal(9, "Pastel1")[col_idx]

pdf("mice-inject/plots/mod78_98_mean_cor_permutation_test.pdf", height=4, width=10)
par(mfrow=c(2, 6))
for (i in 1:length(cor_tests78)) {
	plotPermTest(cor_tests78[[i]], main=names(cor_tests78)[i],
		xlab="Module 78 mean |r|",
		border=colors[i], col=colors_light[i])
}

for (i in 1:length(cor_tests98)) {
	plotPermTest(cor_tests98[[i]], main=names(cor_tests98)[i],
		xlab="Module 98 mean |r|",
		border=colors[i], col=colors_light[i])
}
dev.off()


# Correlation analysis for module 98 genes
# -----------------------------------------------
cmat = corAndPvalue(t(mat98))

# cmat = corAndPvalue(t(mat98[, sample_idx$FSTL3]))
# cmat = corAndPvalue(t(mat98[, sample_idx$Veh]))
# cmat = corAndPvalue(t(mat98[, sample_idx$LBP]))
# cmat = corAndPvalue(t(mat98[, sample_idx$FCN2]))
# cmat = corAndPvalue(t(mat98[, sample_idx$EPDR1]))

# Bonferoni with effective parameters
cmat$padj = matrix(cmat$p * sum(lower.tri(cmat$p)), nrow=nrow(cmat$p))
diag(cmat$padj) = NA

sum(cmat$padj[lower.tri(cmat$padj)] < 0.05, na.rm=TRUE)

pdf("mice-inject/plots/mod98_cor_heatmap.pdf", height=5)
note = matrix("", nrow=nrow(cmat$cor), ncol=ncol(cmat$cor))
note[cmat$padj < 0.05] = "*"

heatmap.2(cmat$cor,
	mar=c(4, 13),
	xlab="Module 98 genes",
	ylab="Module 98 genes",
	col=colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
	key.title="",
	key.xlab="Pearson cor.",
	cellnote=note,
	notecol="black",
	cexRow=0.5, cexCol=0.5,
	trace="none",
	tracecol="black"
)
dev.off()


# Module 98 gene PCA
# Uses normalization defined above
# ----------------------------------------------
pca98 = prcomp(t(mat98), scale=TRUE)

pdf("mice-inject/plots/mod98_pca.pdf", width=3.7, height=4)
colors = c(brewer.pal(9, "Set1")[1:4], "white")
conditions = factor(samples$condition, levels=c("EPDR1", "FCN2", "LBP", "FSTL3", "Veh"))
pt_col = colors[as.integer(conditions)]


plot(pca98$x[, 1], pca98$x[, 2],
	main=paste0("Module 98 genes (g=", nrow(mat98), ")"),
	xlab=paste0("PC1 (", summary(pca98)$importance[2, 1] * 100, "%)"),
	ylab=paste0("PC2 (", summary(pca98)$importance[2, 2] * 100, "%)"),
	bg=pt_col, pch=21)

legend("topleft", legend=levels(conditions), pch=21, pt.bg=colors)
dev.off()

# Individual correlations
# ---------------------------------------------

gene2 = "Dhcr7"
gene1 = "Pcsk9"

# gene1 = "Hmgcs1"

pdf(paste0("mice-inject/plots/scatter_", gene1, gene2, ".pdf"),
	width=3.7, height=4.0)
i = rownames(count_norm) == gene1
j = rownames(count_norm) == gene2


colors = brewer.pal(9, "Set1")[c(1:4, 9)]
levels(samples$condition)
pts_col = colors[as.integer(samples$condition)]

x = count_norm[i, ]
y = count_norm[j, ]

cor_test = cor.test(x, y)

plot(x, y,
	main=paste0(
		"P=", format(cor_test$p.value, digits=3),
		", r=", format(cor_test$estimate, digits=3)
	),
	xlab=paste0(gene1, " (adj. counts)"),
	ylab=paste0(gene2, " (adj. counts)"),
	bty="n",
	bg=pts_col,
	pch=21)

legend("bottomright",
	legend=levels(samples$condition),
	pch=21,
	pt.bg=colors)
dev.off()



# Bayesian network edge correlations
# ------------------------------------------
library(plyr)

# Load bayesian network
bnet = fread("co-expression/annotate/bayesNet/all.tsv")
bnet_nodes = fread("co-expression/annotate/bayesNet/nodes.tsv")

# Subnetwork 98
bnet98 = bnet[
	bnet$TAIL %in% bnet_nodes$NODE[bnet_nodes$MODULE == 98] |
	bnet$HEAD %in% bnet_nodes$NODE[bnet_nodes$MODULE == 98],
]

bnet98$from = sapply(strsplit(bnet98$HEAD, "_"), function(x) x[2])
bnet98$to = sapply(strsplit(bnet98$TAIL, "_"), function(x) x[2])

# Copy expression matrix, for human rownames (gene symbols)
cpm_norm_human = cpm_norm
rownames(cpm_norm_human) = mouse2human(rownames(cpm_norm))


cpm_norm_human = log2(cpm_norm_human + 1)

netw_idx = cbind(
	match(bnet98$from, rownames(cpm_norm_human)),
	match(bnet98$to, rownames(cpm_norm_human))
)

complete = apply(netw_idx, 1, function(idx) !any(is.na(idx)))
netw_idx = netw_idx[complete, ]



# Veh + all perturbations
sample_idx = samples$condition != "Veh"
cor_tests_pert = alply(netw_idx, 1, function(row_pairs) {
	i = row_pairs[1]
	j = row_pairs[2]
	cor_test = cor.test(cpm_norm_human[i, sample_idx], cpm_norm_human[j, sample_idx])
	cor_test$id = paste(rownames(cpm_norm_human)[i], rownames(cpm_norm_human)[j], sep="_")
	return(cor_test)
})


# Veh only
sample_idx = samples$condition == "Veh"
cor_tests_veh = alply(netw_idx, 1, function(row_pairs) {
	i = row_pairs[1]
	j = row_pairs[2]
	cor_test = cor.test(cpm_norm_human[i, sample_idx], cpm_norm_human[j, sample_idx])
	cor_test$id = paste(rownames(cpm_norm_human)[i], rownames(cpm_norm_human)[j], sep="_")
	return(cor_test)
})


pdf("mice-inject/plots/regulatory_interactions_98_boxplot.pdf", width=2.5, height=4)
# x = sapply(cor_tests_veh, function(x) x$p.value)
# x = -log10(x)

x = sapply(cor_tests_veh, function(x) x$estimate)
x = abs(x)

# y = sapply(cor_tests_veh, function(x) x$p.value)
# y = -log10(x)

y = sapply(cor_tests_pert, function(x) x$estimate)
y = abs(y)

t_test = t.test(x, y, paired=TRUE)
colors = brewer.pal(9, "Set1")
boxplot(list(Veh=x, "FSTL3+LBP+EPDR1+FCN2"=y),
	# ylab="|r|",
	las=2,
	frame=FALSE,
	# col=colors,
	ylab=paste0("Module 98 regulatory interactions (|r|, m=", nrow(netw_idx), ")"),
	main=paste0("P=", format(t_test$p.value, digits=4))
)

points(rep(1, length(x)), x, col=colors[1], cex=0.5)
points(rep(2, length(y)), y, col=colors[2], cex=0.5)

segments(1, x, 2, y, col=rgb(0, 0, 0, 0.3))
dev.off()