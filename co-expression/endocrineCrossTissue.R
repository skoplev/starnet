# 
rm(list=ls())

library(data.table)
library(WGCNA)
library(gplots)
library(RColorBrewer)

library(devtools)  # for installing heatmap.3
# Load heatmap.3
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")


data_dir = "/Users/sk/DataProjects/cross-tissue"  # root of data directory
setwd("/Users/sk/Google Drive/projects/STARNET/cross-tissue")
source("src/parse.R")

addAlpha = function(col, alpha=1){
	if (missing(col)) stop("Please provide a vector of colours.")
	apply(sapply(col, col2rgb)/255, 2, 
		function(x) 
		rgb(x[1], x[2], x[3], alpha=alpha)
	)  
}

# Load cross-tissue modules
# ------------------------------
between = new.env()
load(file.path(data_dir, "modules/between_within-cross-tissue.RData"),
	between,
	verbose=TRUE)

# Parse module data
between = parseModuleData(between)

# Main module table
mod_tab = fread("co-expression/tables/module_tab.csv")

secreted_proteins = fread(file.path(data_dir, "Uniprot/uniprot_human_secreted_proteins.tab"))

sec_prots = secreted_proteins[["Gene names  (primary )"]]
sec_prots = strsplit(paste(sec_prots, collapse=";"), ";")[[1]]
sec_prots = trimws(sec_prots)
sec_prots = unique(sec_prots)
sec_prots = sec_prots[sec_prots != ""]


# Load expression data
# -----------------------------------
emat_file = "STARNET/gene_exp_norm_reshape/expr_recast.RData"

load(file.path(data_dir, emat_file), verbose=TRUE)

emat = expr_recast[, 3:ncol(expr_recast)]
meta_genes = expr_recast[, 1:2]
meta_genes = as.data.frame(meta_genes)
meta_genes$gene_symbol = sapply(strsplit(as.character(meta_genes$transcript_id), "_"), function(x) x[1])

rm(expr_recast)

# Test if the gene metadata is the same as for the cross-tissue modules
if (!all(meta_genes$transcript_id == between$meta_genes$transcript_id)) {
	stop("Transcript mismatch")
}


# Secretory proteins in cross-tissue modules
# --------------------------------------
# sum(sec_prots %in% between$meta_genes$gene_symbol)  # total number of mapped secretory protein symbols

# Find secretory proteins in cross-tissue modules
# length(sec_prots)

# sec_idx = between$meta_genes$gene_symbol %in% sec_prots
# sum(sec_idx)

cross_tissue_frac_defn = 0.05  # definition of cross-tissue modules
min_mod_size = 5000

cross_tissue_modules = which(mod_tab$purity < 1 - cross_tissue_frac_defn & mod_tab$mod_size < min_mod_size)
tissue_specific_modules = which(mod_tab$purity >= 1 - cross_tissue_frac_defn & mod_tab$mod_size < min_mod_size)

CT_idx = between$meta_genes$gene_symbol %in% sec_prots & between$clust %in% cross_tissue_modules
TS_idx = between$meta_genes$gene_symbol %in% sec_prots & between$clust %in% tissue_specific_modules




# table(between$meta_genes$tissue[CT_idx])

# library("Hmisc")

# Endocrine-eigengene correlations
cmat = cor(t(emat[CT_idx, ]), between$bwnet$eigengenes,
	use="pairwise.complete.obs")



cmat_ts = cor(t(emat[TS_idx, ]), between$bwnet$eigengenes,
	use="pairwise.complete.obs")


# Table for secreted proteins
sec_info = data.frame(
	tissue=between$meta_genes$tissue[CT_idx],
	gene_symbol=between$meta_genes$gene_symbol[CT_idx],
	clust=between$clust[CT_idx])
sec_info$id = paste0(sec_info$tissue, ":", sec_info$gene_symbol)

sec_info$cross_tissue = sec_info$clust %in% cross_tissue_modules

table(sec_info$tissue)
length(unique(sec_info$gene_symbol))



CT_self = cbind(1:sum(CT_idx), between$clust[CT_idx])
cmat_no_self = cmat
cmat_no_self[CT_self] = NA

tissue_count = mod_tab[, c("AOR", "BLOOD", "LIV", "MAM", "SKLM", "SF", "VAF")]

# Primary tissue of each module
primary_tissue = apply(tissue_count, 1, function(row) {
	names(row)[which.max(row)]
})

TS_self = cbind(1:sum(TS_idx), between$clust[TS_idx])
cmat_ts_no_self = cmat_ts
cmat_ts_no_self[TS_self] = NA


qqplotAnnot = function(x, y,
	probs=c(0.001, 0.01, 0.1, 0.9, 0.99, 0.999),
	...)
{
	qqplot(
		x,
		y,
		pch=16, cex=0.5,
		...
	)
	abline(0, 1, col="red")

	q_y = quantile(y, probs=probs, na.rm=TRUE)
	q_x = quantile(x, probs=probs, na.rm=TRUE)

	points(q_x, q_y, col="grey", cex=0.6)
	text(q_x, q_y, label=paste0(probs*100, "%"), pos=1, cex=0.8, col="grey")
}


# Are CT modules have more distinct associations with endocrine factors than TS modules?
pdf("co-expression/plots/endocrine/eigengene_cor_same_module.pdf")
par(mfrow=c(2, 1))
hist(cmat[CT_self], xlim=c(-1, 1), breaks=50, col=brewer.pal(9, "Set1")[1],
	main=paste0("Cross-tissue endocrines (n=", nrow(CT_self), ")"), xlab="Endocrine-eigengene cor., self-module")
hist(cmat_ts[TS_self], xlim=c(-1, 1), breaks=50, col=brewer.pal(9, "Set1")[2],
	main=paste0("Tissue-specific endocrines (n=", nrow(TS_self), ")"), xlab="Endocrine-eigengene cor., self-module")
dev.off()


# # Matches the cmat index for cross-tissue endocrines
# bool_mat_CT = sapply(as.character(meta_genes$tissue[CT_idx]), function(tissue) {
# 	return(primary_tissue != tissue)
# })
# bool_mat_CT = t(bool_mat_CT)

# # array indices of endocrine-module associations where the endocrine is not the
# # from the primary (most prevalent) tissue.
# arr_ind_CT = which(bool_mat_CT, arr.ind=TRUE)

# # Only cross-tissue modules
# arr_ind_CT = arr_ind_CT[arr_ind_CT[, 2] %in% tissue_specific_modules, ]


# bool_mat_TS = sapply(as.character(meta_genes$tissue[TS_idx]), function(tissue) {
# 	return(primary_tissue != tissue)
# })
# bool_mat_TS = t(bool_mat_TS)

# # array indices of endocrine-module associations where the endocrine is not the
# # from the primary (most prevalent) tissue.
# arr_ind_TS = which(bool_mat_TS, arr.ind=TRUE)

# # Only cross-tissue modules
# arr_ind_TS = arr_ind_TS[arr_ind_TS[, 2] %in% tissue_specific_modules, ]

# par(mfrow=c(2, 1))
# hist(cmat[arr_ind_CT], breaks=50, xlim=c(-1, 1), main="Cross-tissue endocrines", xlab="Heterogenous, tissue-specific eigengene cor.")
# hist(cmat_ts[arr_ind_TS], breaks=50, xlim=c(-1, 1), main="Tissue-specific endocrines")


pdf("co-expression/plots/endocrine/eigengene_qqplots.pdf", width=10, height=5.5)
par(mfrow=c(1, 2))
qqplotAnnot(
	cmat_ts[TS_self],
	cmat[CT_self],
	main="QQ-plot, self-module eigengene cor.",
	xlab="Tissue-specific endocrines",
	ylab="Cross-tissue endocrines"
)
abline(0, 1, col="red")

qqplotAnnot(
	cmat_ts[arr_ind_TS],
	cmat[arr_ind_CT],
	main="QQ-plot, TS eigengene cor., different tissue",
	xlab="Tissue-specific endocrines",
	ylab="Cross-tissue endocrines"
)
dev.off()

wilcox.test(
	abs(cmat_ts[TS_self]),
	abs(cmat[CT_self])
)

wilcox.test(
	abs(cmat_ts[arr_ind_TS]),
	abs(cmat[arr_ind_CT])
)


# Overlap with Bayesian network
library(igraph)

# Load directed network
netw = new.env()
load("~/DataProjects/cross-tissue/R_workspaces/BayesNet2.RData", netw, verbose=TRUE)

idx = CT_idx | TS_idx
cmat_all = cor(t(emat[idx, ]), between$bwnet$eigengenes,
	use="pairwise.complete.obs")

rownames(cmat_all) = paste0(between$meta_genes$tissue, ":", between$meta_genes$gene_symbol)[idx]


# Find module of each endocrine factor
endocrine_module = between$clust[idx]

g = graph_from_graphnel(netw$bn$graphNEL)
edge_list = as_edgelist(g)

# Get endocrine-module associations above threshold
endocrine_edges = apply(edge_list, 1, function(edge) {
	from = as.numeric(edge[1])
	to = as.numeric(edge[2])

	cmat_sub = cmat_all[endocrine_module == from, to, drop=FALSE]
	cmat_sub = as.vector(cmat_sub)
	names(cmat_sub) = rownames(cmat_all)[endocrine_module == from]
	# rownames(cmat_sub) = rownames(cmat_all)
	cmat_sub = cmat_sub[abs(cmat_sub) > 0.2, drop=FALSE]
	return(cmat_sub)
})

# example
edge = c(78, 98)
# edge = c(98, 78)
# edge = c(145, 33)

k = which(edge_list[, 1] == edge[1] & edge_list[, 2] == edge[2])

edge_list[k, ]

pdf(paste0("co-expression/plots/endocrine/netwOverlap/", paste(edge, collapse="_"), ".pdf"), width=10, height=5.5)
par(mar=c(8, 4.1, 4.1, 2.1))
barplot(sort(endocrine_edges[k][[1]], decreasing=TRUE), las=2, 
	ylab="Eigengene cor.",
	main=paste(edge, collapse="->"))
dev.off()


# origin = 78
# origin = 28
# origin = 150
# origin = 122
# origin = 106
origin = 33

pdf(paste0("co-expression/plots/endocrine/netwOverlap/origin_", origin, ".pdf"), width=20, height=12)
edges = which(edge_list[, 1] == origin)
par(mfrow=c(4, 4))

for (k in edges) {
	edge = edge_list[k, ]
	try({
		barplot(sort(endocrine_edges[k][[1]], decreasing=TRUE), las=2, 
			ylab="Eigengene cor.",
			main=paste(edge, collapse="->"),
			cex.names=0.5
		)
	})
}
dev.off()


# Intersection 



# qqplot(
# 	cmat_ts_no_self[, tissue_specific_modules],
# 	cmat_no_self[, tissue_specific_modules],
# 	main="QQ-plot, TS eigengene cor. (no self)",
# 	xlab=paste0("Tissue-specific endocrines (n=", nrow(cmat_ts), ")"),
# 	ylab=paste0("Cross-tissue endocrines (n=", nrow(cmat), ")"),
# 	pch=16,
# 	cex=0.5,
# 	xlim=c(-1, 1),
# 	ylim=c(-1, 1)
# )
# abline(0, 1, col="red")

# wilcox.test(
# 	cmat_no_self[, tissue_specific_modules],
# 	cmat_ts_no_self[, tissue_specific_modules]
# )

# # probs = seq(0, 1, 0.05)
# probs = c(0.001, 0.01, 0.1, 0.9, 0.99, 0.999)
# q_CT = quantile(cmat_no_self[, tissue_specific_modules], probs=probs, na.rm=TRUE)
# q_TS = quantile(cmat_ts_no_self[, tissue_specific_modules], probs=probs, na.rm=TRUE)

# points(q_TS, q_CT, col="grey", cex=0.6)
# text(q_TS, q_CT, label=paste0(probs*100, "%"), pos=1, cex=0.8, col="grey")

# sum(cmat_no_self > 0.7, na.rm=TRUE)
# sum(cmat_no_self > 0.7, na.rm=TRUE)



sum(cmat_ts_no_self > 0.7, na.rm=TRUE)


# Example correlation
tissue = "MAM"
gene = "PLXNB1"
mod = 1

tissue = "LIV"
gene = "LTBP4"
mod = 171


cmat_no_self[635, 103]



# x = data.matrix(t(emat[CT_idx, ]))
# y = data.matrix(between$bwnet$eigengenes)

# cmat_test = rcorr(x, y)   # correlations with p-values, concatenates x and y, which must be separated

# n = ncol(cmat_test$P)  # total features (endocrine + eigengenes)
# k = ncol(x)  # endocrine features, 

# pmat = cmat_test$P[1:k, (k + 1):n]  # matrix of endocrine-eigengene correlation p-values
# cmat_qval = qvalue(pmat)

# sum(pmat < 0.05)
# sum(cmat_qval$qvalue < 0.05)
# sum(cmat_qval$qvalue < 0.05) / length(cmat_qval$qvalue)


# pmat[, tissue_specific_modules]

# sum(cmat_qval$qvalue[, tissue_specific_modules] < 0.05)
# sum(cmat_qval$qvalue[, tissue_specific_modules] < 0.05 & abs(cmat[, tissue_specific_modules]) > 0.3)

# sum(abs(cmat[, tissue_specific_modules]) > 0.3)

rownames(cmat) = paste0(between$meta_genes$tissue, ":", between$meta_genes$gene_symbol)[sec_idx]





# cmat[sec_info$cross_tissue, tissue_specific_modules]
# cmat[!sec_info$cross_tissue, tissue_specific_modules]

# plot(density(cmat))

# sum(cmat > 0.8)

# include = apply(abs(cmat), 1, max) > 0.7

# # sum(apply(abs(cmat), 1, max) > 0.5)
# mat = cmat[include, ]
# mat = cmat[include, mod_tab$CAD_qvalue < 0.1]

# mat = cmat[include, mod_tab$CAD_qvalue < 0.1 & mod_tab$purity < 0.99]  # Cross-tissue
# mat = cmat[include, mod_tab$CAD_qvalue < 0.1 & mod_tab$purity >= 0.99]  # Tissue-specific

# mod = 78
# mod = 145
# mod = 28
# mod = 116
# mod = 150
# mod = 189

# mod = 188
# mod = 33
# mod = 106
# mod = 122
# mod = 162
mod = 203
mat = cmat[sec_info$clust == mod, tissue_specific_modules]
mat = mat[, apply(mat, 2, function(x) max(abs(x))) > 0.3]


tissues = unique(between$meta_genes$tissue)
tissue_col = brewer.pal(9, "Set1")[-6]

# Tissues count matrix for selected modules
tissue_counts = as.data.frame(mod_tab)[as.integer(colnames(mat)), as.character(tissues)]

tissue_frac = tissue_counts / mod_tab$mod_size[as.integer(colnames(mat))]

rlab = matrix("white", ncol=ncol(tissue_counts), nrow=nrow(tissue_counts))
colnames(rlab) = tissues

for (i in 1:ncol(rlab)) {
	for (j in 1:nrow(rlab)) {
		rlab[j, i] = addAlpha(tissue_col[i], tissue_frac[j, i])
		# rlab[j, i] = addAlpha(tissue_col[i], tissue_frac[j, i] + 0.1)
	}
}

library(gplots)

pdf(paste0("co-expression/plots/endocrine/eigengeneCor/tissue_specific/tissue_specific_eigencor_mod", mod, ".pdf"),
	width=8)
heatmap.3(mat,
	trace="none",
	col=colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
	ColSideColors=rlab,
	ColSideColorsSize=2.5,
	breaks=seq(-1.0, 1.0, length.out=101),
	cexRow=0.3,
	xlab="Module", ylab=paste0("Endocrines, mod=", mod),
	KeyValueName="Pearson's r"
)
dev.off()



# 
tissues = unique(between$meta_genes$tissue)
tissue_col = brewer.pal(9, "Set1")[-6]

# mod_tab$LIV < 6

# cor_thresh = 0.25
cor_thresh = 0.20
# cor_thresh = 0.3

# cor_thresh = 0.5

# cor_thresh = 0.7

# t1 = "LIV"
# t1 = "VAF"
# t1 = "AOR"
# t1 = "MAM"
# t1 = "SKLM"
# t1 = "SF"
# t1 = "BLOOD"

for (t1 in tissues) {
	pdf(paste0("co-expression/plots/endocrine/eigengeneCor/", t1, ".pdf"), width=8)
	include_mod = mod_tab$CAD_qvalue < 0.1 | 
		mod_tab$pval_syntax_score < 0.05 |
		mod_tab$pval_ndv < 0.05 |
		mod_tab$pval_lesions < 0.05 |
		mod_tab$pval_DUKE < 0.05

	mod_idx = include_mod & mod_tab[[t1]] < 2
	# mod_idx = include_mod
	# mod_idx = rep(TRUE, length(include_mod))  # all
	# mod_idx = mod_tab[[t1]] < 2

	# gene_idx = between$meta_genes$tissue[sec_idx] == t1  # relative to
	gene_idx = sec_info$tissue == t1

	mat = cmat[gene_idx, mod_idx]  # temp

	mat = mat[apply(abs(mat), 1, max) > cor_thresh, apply(abs(mat), 2, max) > cor_thresh,
		drop=FALSE
	]
	# mat = as.matrix(mat)

	# Look up if each gene is part of module...
	sec_info$clust[match(rownames(mat), sec_info$id)]

	sec_in_module = sapply(colnames(mat), function(k) {
		sec_info$clust[match(rownames(mat), sec_info$id)] == k
	})
	sec_in_module = as.matrix(sec_in_module)

	sec_in_module_note = matrix("", ncol=ncol(sec_in_module), nrow=nrow(sec_in_module))
	sec_in_module_note[sec_in_module] = "."

	# List module statistics from
	mod_tab[as.integer(colnames(mat)), 1:10]

	# Tissues count matrix for selected modules
	tissue_counts = as.data.frame(mod_tab)[as.integer(colnames(mat)), as.character(tissues)]

	tissue_frac = tissue_counts / mod_tab$mod_size[as.integer(colnames(mat))]

	rlab = matrix("white", ncol=ncol(tissue_counts), nrow=nrow(tissue_counts))
	colnames(rlab) = tissues

	for (i in 1:ncol(rlab)) {
		for (j in 1:nrow(rlab)) {
			rlab[j, i] = addAlpha(tissue_col[i], tissue_frac[j, i])
			# rlab[j, i] = addAlpha(tissue_col[i], tissue_frac[j, i] + 0.1)
		}
	}

	clab = t(matrix(rep(tissue_col[match(t1, tissues)], nrow(mat))))

	# Endocrine name only
	rownames(mat) = sapply(strsplit(rownames(mat), ":"), function(x) x[2])

	try({
		heatmap.3(mat,
			trace="none",
			# mar=c(8, 24),
			main=t1,
			mar=c(8, 32),
			col=colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
			ColSideColors=rlab,
			ColSideColorsSize=2.5,
			RowSideColorsSize=0.7,
			RowSideColors=clab,
			# breaks=seq(-0.5, 0.5, length.out=101),
			breaks=seq(-1.0, 1.0, length.out=101),
			# cexRow=0.075,
			cexRow=0.5,
			xlab="Module", ylab="Endocrine",
			KeyValueName="Pearson's r",
			cellnote=sec_in_module_note, notecol="black"
		)
	})
	dev.off()
}



