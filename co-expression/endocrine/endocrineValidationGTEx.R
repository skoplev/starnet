rm(list=ls())

library(data.table)
library(RColorBrewer)
library(WGCNA)
library(compiler)
enableJIT(3)

setwd("~/GoogleDrive/projects/STARNET/cross-tissue")

source("src/permuteTest.R")
source("src/parse.R")
source("src/base.R")


# Load STARNET modules
modules = fread("co-expression/tables/modules.csv")
modules$gene_biotype = parseEnsemblBiotype(modules$ensembl)
modules$tissue_biotype = paste(modules$tissue, modules$gene_biotype, sep="_")
modules$ensembl_base = sapply(strsplit(modules$ensembl, "[.]"), function(x) x[1])
modules$tissue_transcript_id = paste(modules$tissue, modules$ensembl_base, sep="_")


mod_tab = fread("co-expression/tables/module_tab.csv")


# Load GTEx
gtex = list()
gtex$mat = fread("~/DataBases/GTEx/RNA-seq/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct")

gtex$annot = fread("~/DataBases/GTEx/RNA-seq/GTEx_v7_Annotations_SampleAttributesDS.txt")

# Match annotation table to column names
# note that first two columns of mat are dummy IDs
gtex$annot = gtex$annot[match(colnames(gtex$mat), gtex$annot$SAMPID), ]

gtex$annot$gtex.ID = sapply(strsplit(gtex$annot$SAMPID, "-"), function(x) paste(x[1], x[2], sep="-"))


table(gtex$annot$SMTSD)

# Map to STARNET tissues
gtex$annot$tissue[gtex$annot$SMTSD == "Adipose - Subcutaneous"] = "SF"
gtex$annot$tissue[gtex$annot$SMTSD == "Adipose - Visceral (Omentum)"] = "VAF"
gtex$annot$tissue[gtex$annot$SMTSD == "Liver"] = "LIV"
gtex$annot$tissue[gtex$annot$SMTSD == "Artery - Aorta"] = "AOR"
gtex$annot$tissue[gtex$annot$SMTSD == "Muscle - Skeletal"] = "SKLM"
gtex$annot$tissue[gtex$annot$SMTSD == "Whole Blood"] = "BLOOD"


#...

# tissues = c("SF", "VAF", "LIV")
tissues = c("SF", "VAF", "LIV", "SKLM", "BLOOD", "AOR")

gtex$annot$tissue %in% tissues
sum(gtex$annot$tissue %in% tissues)



# GTEx IDs with at 
avail = gtex$annot$tissue %in% tissues
# gtex_ids = names(which(
# 	table(gtex$annot$gtex.ID[avail]) > 1))

gtex_ids = gtex$annot$gtex.ID[avail]

length(gtex_ids)



# Make long table tissue mRNA x patient ID
mats = lapply(tissues, function(tis) {
	message(tis)
	# Tissue samples
	idx = which(gtex$annot$tissue == tis)

	# Match column idx to patients
	idx = idx[match(gtex_ids, gtex$annot$gtex.ID[idx])]

	gtex$annot$tissue[idx]

	gtex$annot$gtex.ID[idx] == gtex_ids

	meta_row = gtex$mat[, 1:2]
	meta_row$tissue = tis

	mat = data.matrix(gtex$mat)[, idx]
	colnames(mat) = gtex_ids

	# colnames(mat) = 
	stopifnot(all(na.omit(colnames(mat) == gtex_ids)))

	tab = cbind(meta_row, data.table(mat))
	return(tab)
})

expr_mat = rbindlist(mats)

meta_row = expr_mat[, 1:3]
meta_row$ensembl_base = sapply(strsplit(meta_row$Name, "[.]"), function(x) x[1])
meta_row$tissue_transcript_id = paste(meta_row$tissue, meta_row$ensembl_base, sep="_")
meta_row$gene_biotype = parseEnsemblBiotype(meta_row$ensembl)
meta_row$tissue_gene_biotype = paste(meta_row$tissue, meta_row$gene_biotype, sep="_")

expr_mat = expr_mat[, -1:-3]
expr_mat = data.matrix(expr_mat)
rownames(expr_mat) = paste(meta_row$tissue, meta_row$ensembl_base, sep="_")

# Remove zero variance genes
# idx = apply(expr_mat, 1, sd, na.rm=TRUE) > 0.001
idx = apply(expr_mat, 1, sd, na.rm=TRUE) > 0
expr_mat = expr_mat[idx, ]
meta_row = meta_row[idx, ]

expr_mat = t(expr_mat)

# sum(apply(expr_mat, 2, sd, na.rm=TRUE) < 0.000001, na.rm=TRUE)
# sum(is.na(apply(expr_mat, 2, sd, na.rm=TRUE)))
	

# Permutation tests
# ----------------------------------------
m = 1000
perm_test = list()

k = 78
mod_idx = modules$clust == k
genes = modules$tissue_transcript_id[mod_idx]
group = modules$tissue_biotype[mod_idx]

perm_test[[k]] = corPermTestExprMat(
	expr_mat=expr_mat,
	genes=genes,
	m=m,
	mat_group=meta_row$tissue_gene_biotype,
	genes_group=group
)



k = 98
mod_idx = modules$clust == k
genes = modules$tissue_transcript_id[mod_idx]
group = modules$tissue_biotype[mod_idx]

perm_test[[k]] = corPermTestExprMat(
	expr_mat=expr_mat,
	genes=genes,
	m=m,
	mat_group=meta_row$tissue_gene_biotype,
	genes_group=group
)

pdf("co-expression/plots/endocrine/validationGTEx/perm_tests/perm_test_78_98.pdf", width=3, height=5)
par(mfrow=c(2, 1))
plotPermuteTest(perm_test[[78]], main="GTEx module 78")
plotPermuteTest(perm_test[[98]], main="GTEx module 98")
dev.off()



# Corelation matrix tests
# --------------------------------------------------

# expr_mat = t(expr_mat)

gtex$mat = t(expr_mat)

# Load STARNET expression data
data_dir = "~/DataProjects/cross-tissue"  # root of data directory

load(file.path(data_dir, "STARNET/gene_exp_norm_reshape/expr_recast.RData"),
	verbose=TRUE)

starnet = parseExprTable(expr_recast)
rm(expr_recast)

# Rename matrix rows
rownames(starnet$mat) = paste(starnet$meta_row$tissue, starnet$meta_row$ensembl_base, sep="_")


k = 98


mod_k = 1:max(modules$clust)
# mod_k = c(78, 98)
# 1:max(modules$clust)
cor_tests = lapply(mod_k, function(k) {
	message(k)

	mod_idx = modules$clust == k
	genes = modules$tissue_transcript_id[mod_idx]
	tissue = modules$tissue[mod_idx]

	out = list()

	# All
	out$all = corModuleTest(mat1=starnet$mat, mat2=gtex$mat,
		module_genes=genes)

	if (length(table(tissue)) > 1) {
		# Cross-tissue
		out$ct = corModuleTest(mat1=starnet$mat, mat2=gtex$mat,
			module_genes=genes,
			method="CT")

		for (tis in names(which(table(tissue) > 20))) {
			message(tis)
			genes = modules$tissue_transcript_id[mod_idx & modules$tissue == tis]
			out[[tis]] = corModuleTest(mat1=starnet$mat, mat2=gtex$mat, module_genes=genes)
		}
	}
	return(out)
})
names(cor_tests) = mod_k

save(cor_tests, file="~/DataProjects/STARNET/moduleVali/valGTEx.RData")


xlab = "STARNET cor."
ylab = "GTEx cor."
pdf("co-expression/moduleValidation/GTEx/gtex_98_LIV.pdf")
plotCorTest(cor_tests[["98"]]$LIV, xlab=xlab, ylab=ylab, main="Module 98 LIV")
dev.off()

pdf("co-expression/moduleValidation/GTEx/gtex_78_CT.pdf")
plotCorTest(cor_tests[["78"]]$ct, xlab=xlab, ylab=ylab, main="Module 78 cross-tissue")
dev.off()

pdf("co-expression/moduleValidation/GTEx/gtex_78_VAF.pdf")
plotCorTest(cor_tests[["78"]]$VAF, xlab=xlab, ylab=ylab, main="Module 78 VAF")
dev.off()

pdf("co-expression/moduleValidation/GTEx/gtex_78_LIV.pdf")
plotCorTest(cor_tests[["78"]]$LIV, xlab=xlab, ylab=ylab, main="Module 78 LIV")
dev.off()

pdf("co-expression/moduleValidation/GTEx/gtex_78_SF.pdf")
plotCorTest(cor_tests[["78"]]$SF, xlab=xlab, ylab=ylab, main="Module 78 SF")
dev.off()


r_val_all = sapply(cor_tests, function(x) {
	result = NA
	try({result = x$all$test$estimate}, silent=TRUE)
	return(result)
})

p_val_all = sapply(cor_tests, function(x) {
	result = NA
	try({result = x$all$test$p.value}, silent=TRUE)
	return(result)
})


n_genes = sapply(cor_tests, function(x) {
	result = NA
	try({result = x$all$ngenes}, silent=TRUE)
	return(result)
})

r_val_ct = sapply(cor_tests, function(x) {
	out = NA  # default
	try({out = x$ct$test$estimate}, silent=TRUE)
	if (is.null(out)) {
		out = NA
	}
	return(out)
})

p_val_ct = sapply(cor_tests, function(x) {
	out = NA  # default
	try({out = x$ct$test$p.value}, silent=TRUE)
	if (is.null(out)) {
		out = NA
	}
	return(out)
})





mod_val = data.frame(
	clust=mod_k,
	mod_size=mod_tab$mod_size,
	r=r_val_all,
	R2=r_val_all^2,
	p=p_val_all,
	r_ct=r_val_ct,
	R2_ct=r_val_ct^2,
	p_ct=p_val_ct,
	n_genes=n_genes
)
write.csv(mod_val, file="~/DataProjects/STARNET/moduleVali/valGTEx_mod_val.csv", row.names=FALSE)

stopifnot(all(mod_tab[, 1] == mod_val$clust))

# mod_tab
tissue_cols = c("AOR", "MAM", "LIV", "VAF", "SF", "SKLM", "BLOOD")
colors = brewer.pal(9, "Set1")[c(1, 5, 7, 4, 8, 3, 2)]
colors_faint = brewer.pal(9, "Pastel1")[c(1, 5, 7, 4, 8, 3, 2)]
# colors_faint = brewer.pal(9, "Set1")

mod_val$primary_tissue = tissue_cols[apply(data.frame(mod_tab)[, tissue_cols], 1, which.max)]

cross_tissue_idx = mod_tab$purity < 0.95

# Tissue-specific
mod_val_ts = mod_val[!cross_tissue_idx, ]

mod_val_ts = mod_val_ts[order(mod_val_ts$R2, decreasing=TRUE), ]



# Exclude MAM tissues
mod_val_ts = mod_val_ts[mod_val_ts$primary_tissue != "MAM", ]

# Exclude tissues with few genes
mod_val_ts = mod_val_ts[which(mod_val_ts$n_genes > 20), ]

# Group by primary tissue
mod_val_ts = mod_val_ts[order(mod_val_ts$primary_tissue), ]

mod_val_ts = mod_val_ts[order(factor(mod_val_ts$primary_tissue, levels=tissue_cols)), ]


# Cross-tissue
mod_val_ct = mod_val[cross_tissue_idx, ]
mod_val_ct = mod_val_ct[order(mod_val_ct$R2_ct, decreasing=TRUE), ]

mod_val_ct = mod_val_ct[mod_val_ct$n_genes > 20, ]




# Is the difference in R2 due to module size?

pdf("co-expression/moduleValidation/GTEx/module_size_test.pdf", width=3.5, height=3.5)

plot(
	log10(mod_val_ts$mod_size),
	mod_val_ts$R2,
	pch=21,
	# col="black",
	# bg="grey",
	bg=brewer.pal(9, "Pastel1")[2],
	xlim=range(log10(mod_val$mod_size)),
	ylim=range(c(mod_val$R2, mod_val$R2_ct), na.rm=TRUE),
	xlab="Module size (log10 n)",
	ylab="R2"
)

points(log10(mod_val_ct$mod_size), mod_val_ct$R2_ct,
	pch=22,
	# bg="red"
	bg=brewer.pal(9, "Pastel1")[1],
)

fit_ts = lm(mod_val_ts$R2 ~ log10(mod_val_ts$mod_size))
abline(fit_ts, col=brewer.pal(9, "Set1")[2], lwd=2)
summary(fit_ts)

fit_ct = lm(mod_val_ct$R2 ~ log10(mod_val_ct$mod_size))
abline(fit_ct, col=brewer.pal(9, "Set1")[1], lwd=2)

legend("topright", legend=c(
	paste0(
		"p=",
		format(
			cor.test(log10(mod_val_ts$mod_size), mod_val_ts$R2)$p.value,
			digits=3
		)
	),
	paste0(
		"p=",
		format(
			cor.test(log10(mod_val_ct$mod_size), mod_val_ct$R2_ct)$p.value,
			digits=3
		)
	)
))
dev.off()

library(plotrix)


pdf("co-expression/moduleValidation/GTEx/all.pdf", width=8, height=9)

par(mfcol=c(4, 1))

alpha = 0.05   # for Bonferoni correction
# Number of used comparisons for Bonferoni correction
model_size = 224

# Tissue-specific
# -----------------------------
bar_cols = colors[as.integer(factor(mod_val_ts$primary_tissue, levels=tissue_cols))]

# Less strong color below R2 threshold
idx_faint = mod_val_ts$R2 < 0.2
bar_cols[idx_faint] = addAlpha(bar_cols[idx_faint], 0.5)

barplot(
	mod_val_ts$R2,
	names.arg=mod_val_ts$clust,
	las=2,
	col=bar_cols,
	cex.names=0.8,
	xlab=paste0("STARNET tissue-specific co-expression modules (n=", nrow(mod_val_ts), ")"),
	ylab=expression("GTEx variance explained R"^2),
	space=0
)


sig_col = rep("white", nrow(mod_val_ts))
sig_idx = mod_val_ts$p < (0.05 / model_size)
sig_col[sig_idx] = "black"

mean(sig_idx) * 100


plot(1:nrow(mod_val_ts), rep(0, nrow(mod_val_ts)),
	col=sig_col,
	pch=16
)

# Cross-tissue
# -------------------------------
grey_col = brewer.pal(9, "Set1")[9]
bar_cols = rep(grey_col, nrow(mod_val_ct))
idx_faint = which(mod_val_ct$R2_ct < 0.1)

bar_cols[idx_faint] = addAlpha(bar_cols[idx_faint], 0.5)

# barplot(
bp = gap.barplot(
	mod_val_ct$R2_ct,
	gap=c(0.32, 0.61),
	xaxt='n',  # suppress default x-axis for gap.barplot
	ytics=c(0,0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8),
	# ylim=c(0, 0.5),

	# # for base barplot()
	# names.arg=mod_val_ct$clust,
	# las=2,
	# cex.names=0.8,

	col=bar_cols,
	xlab=paste0("STARNET cross-tissue co-expression modules (n=", nrow(mod_val_ct), ")"),
	ylab=expression("GTEx variance explained R"^2),
	space=0
)

# Custom barplot axis for use with gap.barplot()
axis(1, at=bp, lab=mod_val_ct$clust,
	las=2,  # vertical labels
	cex.axis=0.8,
	tck=0  # no ticks
)

sig_col = rep("white", nrow(mod_val_ct))
sig_idx = mod_val_ct$p_ct < (0.05 / model_size)
sig_col[sig_idx] = "black"

mean(sig_idx, na.rm=TRUE) * 100

plot(1:nrow(mod_val_ct), rep(0, nrow(mod_val_ct)),
	col=sig_col,
	pch=16
)

dev.off()



# Co-expression module 78
# validation of cross-tissue correlations. Only adipose to liver axis
# --------------------------------------------------------
k = 78

mod_idx = modules$clust == k
genes = modules$tissue_transcript_id[mod_idx]
tissue = modules$tissue[mod_idx]


# Get transcript IDs found in both GTEx and STARNET
genes_found = genes[
	genes %in% rownames(starnet$mat) &
	genes %in% rownames(gtex$mat)
]

tissue_genes_found = sapply(strsplit(genes_found, "_"), function(x) x[1])

adipose_genes_found =  genes_found[tissue_genes_found %in% c("SF", "VAF")]
liver_genes_found = genes_found[tissue_genes_found == "LIV"]


# adipose-liver cross correlations
cmat_starnet = corAndPvalue(
	t(starnet$mat[rownames(starnet$mat) %in% adipose_genes_found, ]),
	t(starnet$mat[rownames(starnet$mat) %in% liver_genes_found, ]),
	use="pairwise.complete.obs"
)

# matching adipose-liver cross correlations for GTEx

cmat_gtex = corAndPvalue(
	t(gtex$mat[rownames(gtex$mat) %in% adipose_genes_found, ]),
	t(gtex$mat[rownames(gtex$mat) %in% liver_genes_found, ]),
	use="pairwise.complete.obs"
)

# r_cutoff = 0.3

# cmat_starnet$cor > r_cutoff

# sum(
# 	cmat_starnet$cor > r_cutoff & cmat_gtex$cor > r_cutoff
# )

# sum(cmat_starnet$cor > r_cutoff & cmat_gtex$cor > r_cutoff)

alpha = 0.05
# r_cutoff = 0.4
r_cutoff = 0.4

mat_bool = cmat_starnet$p < (alpha / length(cmat_starnet$p)) &
	cmat_gtex$p < (alpha / length(cmat_gtex$p)) &
	abs(cmat_starnet$cor) > r_cutoff & abs(cmat_gtex$cor) > r_cutoff &
	# same direction
	sign(cmat_starnet$cor) == sign(cmat_gtex$cor)

sum(mat_bool)

mat_index = which(mat_bool, arr.ind=TRUE)

edges = as.data.frame(mat_index, row.names="")

edges$from_id = rownames(cmat_starnet$cor)[mat_index[, 1]]
edges$to_id = colnames(cmat_starnet$cor)[mat_index[, 2]]

edges$starnet_cor = cmat_starnet$cor[mat_index]
edges$gtex_cor = cmat_gtex$cor[mat_index]

edges$from_tissue = sapply(strsplit(edges$from_id, "_"), function(x) x[1])
edges$to_tissue = sapply(strsplit(edges$to_id, "_"), function(x) x[1])

edges$from_gene_symbol = sapply(edges$from_id, function(id) {
	ensembl = sapply(strsplit(id, "_"), function(x) x[2])
	gene_symbol = starnet$meta_row$gene_symbol[match(ensembl, starnet$meta_row$ensembl_base)]
	return(gene_symbol)
})

edges$to_gene_symbol = sapply(edges$to_id, function(id) {
	ensembl = sapply(strsplit(id, "_"), function(x) x[2])
	gene_symbol = starnet$meta_row$gene_symbol[match(ensembl, starnet$meta_row$ensembl_base)]
	return(gene_symbol)
})


library(igraph)


g = graph_from_edgelist(
	as.matrix(edges[, c("from_id", "to_id")]),
	directed=FALSE)

# E(g)$cor = edges$

E(g)$color = brewer.pal(9, "Set1")[as.integer(edges$starnet_cor > 0)]

V(g)$tissue = sapply(
	strsplit(as_ids(V(g)), "_"), function(x) x[1])

V(g)$type[V(g)$tissue == "LIV"] = "liver"
V(g)$type[V(g)$tissue %in% c("SF", "VAF")] = "adipose"

V(g)$gene_symbol = sapply(as_ids(V(g)), function(id) {
	ensembl = sapply(strsplit(id, "_"), function(x) x[2])
	gene_symbol = starnet$meta_row$gene_symbol[match(ensembl, starnet$meta_row$ensembl_base)]
	return(gene_symbol)
})

V(g)$color = brewer.pal(9, "Set1")[c(4, 8, 7)][match(V(g)$tissue, c("VAF", "SF", "LIV"))]
# V(g)


pdf("co-expression/moduleValidation/GTEx/module_78_bipartite_graph.pdf",
	width=2.0, height=5
)
lay = layout_as_bipartite(g,
	types=(V(g)$type == "adipose")
)
# rotate
lay = lay[, 2:1]
plot(g,
	layout=lay,
	asp=0,
	rescale=FALSE,
	xlim=range(lay[, 1]),
	ylim=range(lay[, 2]),
	vertex.size=15,
	vertex.label=V(g)$gene_symbol,
	vertex.label.family="Helvetica",
	vertex.label.dist=7,
	vertex.label.degree=-pi/2,  # top
	vertex.label.color="black"
)
dev.off()



# Check correlations
i = which(rownames(starnet$mat) == "VAF_ENSG00000106772")  # PRUNE2
j = which(rownames(starnet$mat) == "LIV_ENSG00000115457")  # IGFBP2

x = starnet$mat[i, ]
y = starnet$mat[j, ]
plot(x, y)

cor.test(x, y)
