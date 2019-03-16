rm(list=ls())

library(data.table)
library(RColorBrewer)
library(DESeq2)
library(biomaRt)
library(limma)  # for vennDiagram
library(gplots)


setwd("~/GoogleDrive/projects/STARNET/cross-tissue")

counts = fread("mice-inject/data/feature_counts/gene_counts.txt")

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

endocrines = c("EPDR1", "FCN2", "LBP")
results = lapply(endocrines, function(pert) {
	res = results(dds, contrast=c("condition", pert, "Veh"))

	res$mgi_symbol = getGeneSymbol(rownames(res))
	res$hgnc_symbol = mouse2human(res$mgi_symbol)  # human gene symbols

	res = res[order(res$pvalue), ]
	return(res)
})
names(results) = endocrines


# Get gene lists of differentially expressed genes
fdr = 0.05
sig_genes_human = lapply(results, function(res) {
	sig = res$padj < fdr
	sig_symbols = res$hgnc_symbol[sig]  # homologous human gene symbols

	sig_symbols = na.omit(sig_symbols)
	sig_symbols = unique(sig_symbols)
	sig_symbols = sig_symbols[sig_symbols != ""]
	return(sig_symbols)
})

sig_genes_mouse = lapply(results, function(res) {
	sig = res$padj < fdr
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
	"mice-inject/DEG/mouse_FCN2_LBP.txt")

# Write gene lists
for (cond in names(sig_genes)) {
	write(sig_genes_mouse[[cond]], paste0("mice-inject/DEG/mouse_", cond, ".txt"))
	write(sig_genes_human[[cond]], paste0("mice-inject/DEG/human_", cond, ".txt"))
}

# Print specific resutls
gene = "DHCR7"
gene = "FAM213A"
gene = "RDH11"
gene = "DHCR24"
for (i in 1:3) {
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
kegg_fcn2 = fread("mice-inject/DEG/enrichr/mouse_FCN2_KEGG_2016_table.txt")
kegg_lbp = fread("mice-inject/DEG/enrichr/mouse_LBP_KEGG_2016_table.txt")

terms = c(
	"Non-alcoholic fatty liver disease (NAFLD)_Homo sapiens_hsa04932",
	"Protein processing in endoplasmic reticulum_Homo sapiens_hsa04141"
)

kegg_pvals = sapply(terms, function(x) {
	fcn2_pval = kegg_fcn2[["Adjusted P-value"]][kegg_fcn2$Term == x]
	lbp_pval = kegg_lbp[["Adjusted P-value"]][kegg_lbp$Term == x]

	pvals = c(lbp_pval, fcn2_pval)
	names(pvals) = c("Lbp", "Fcn2")
	return(pvals)
})


pdf("mice-inject/plots/KEGG_enrichrment_LBP_FCN2.pdf", width=2, height=5.0)
par(mar=c(16, 4, 2, 2))
barplot(-log10(kegg_pvals),
	beside=TRUE,
	las=2,
	ylab="-log10 p (BH)",
	legend.text=TRUE,
	col=brewer.pal(9, "Set1")[3:2])

abline(h=-log10(0.05), lty=2, col="grey")
abline(h=0)
dev.off()


# Normalization based on size factors
# ---------------------------------------
sizeFactors(dds)

library(edgeR)

count_norm = sweep(count_mat[include, ], 2, sizeFactors(dds), "/")

rownames(count_norm) = getGeneSymbol(rownames(count_norm))

# cpm_norm = cpm(count_norm)

cpm_norm = cpm(count_mat[include, ])
rownames(cpm_norm) = getGeneSymbol(rownames(cpm_norm))



# Plots of differentially expressed genes
# -----------------------------------------

plotExpr = function(gene, mat, conditions=c("Veh", "EPDR1", "FCN2", "LBP")) {
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

	colors = brewer.pal(9, "Set1")[c(9, 1:3)]

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
	for (i in 1:3) {
		print(results[[i]][which(results[[i]]$mgi_symbol == gene), ])
	}

}


pdf("mice-inject/plots/gene_expr_DEG_mod98.pdf", height=2.1, width=5.5)
par(mfrow=c(1, 4))
plotExpr("Dhcr7", cpm_norm)
plotExpr("Fam213a", cpm_norm)
plotExpr("Rdh11", cpm_norm)
plotExpr("Dhcr24", cpm_norm)
dev.off()


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
	

# values = lapply(values, log2)

boxplot(values, frame=FALSE, las=2)

boxplot(values, las=2)

plot(sapply(values, mean))


plotmeans(values)


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



# Module 98 targeted
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


par(mfrow=c(3, 1))
for (i in 1:3) {
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

# Individual correlations
# ---------------------------------------------

gene2 = "Dhcr7"
gene1 = "Pcsk9"

pdf(paste0("mice-inject/plots/scatter_", gene1, gene2, ".pdf"),
	width=3.5, height=4.0)
i = rownames(count_norm) == gene1
j = rownames(count_norm) == gene2


colors = brewer.pal(9, "Set1")[c(1:3, 9)]
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





