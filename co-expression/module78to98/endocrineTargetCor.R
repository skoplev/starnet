# Calculates gene-level correlations between endocrine factors
# Adopted from endocrineCrossTissueReadouts.R, which calculates global correlations between endocrine factors and target tissue genes.
# 

rm(list=ls())

library(data.table)
library(WGCNA)
library(gplots)
library(RColorBrewer)

setwd("~/GoogleDrive/projects/STARNET/cross-tissue")
data_dir = "~/DataProjects/cross-tissue"  # root of data directory

source("src/parse.R")


# Load cross-tissue modules
# ------------------------------
between = new.env()
load(file.path(data_dir, "modules/between_within-cross-tissue.RData"),
	between,
	verbose=TRUE)

# Parse module data
between = parseModuleData(between)


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


# Load key driver analysis results

# load results table
kda = read.table(
	file.path("co-expression/annotate/bayesNet", "kda", "modules.results.txt"),
	header=TRUE
)

kda$tissue = sapply(strsplit(as.character(kda$NODE), "_"), function(x) x[1])
kda$gene_symbol = sapply(strsplit(as.character(kda$NODE), "_"), function(x) x[2])

# Some remane, auxiliary, columns for merging
kda$target_tissue = kda$tissue
kda$target_clust = kda$MODULE
kda$target_gene_symbol = kda$gene_symbol
kda$key_driver_FDR = kda$FDR


# Load validated endocrine table
endo = fread("co-expression/tables/CT_endocrines_TS_interactions.csv")
endo = endo[endo$clust == 78 & endo$target_clust == 98, ]

endo_target_cor = lapply(1:nrow(endo), function(k) {
	# Find endocrine factor index in gene expression matrix
	endo_idx = which(
		meta_genes$tissue == endo$tissue[k] &
		meta_genes$gene_symbol == endo$gene_symbol[k])

	target_clust_idx = which(between$clust == endo$target_clust[k])

	# Calculate correlation statistics
	cors = corAndPvalue(t(emat[endo_idx, ]), t(emat[target_clust_idx, ]))

	# Make table
	df = data.frame(
		endocrine=endo$gene_symbol[k],
		source_tissue=endo$tissue[k],
		source_clust=endo$clust[k],
		target_tissue=meta_genes$tissue[target_clust_idx],
		target_gene_symbol=meta_genes$gene_symbol[target_clust_idx],
		target_clust=between$clust[target_clust_idx],
		cor=t(cors$cor),
		p=t(cors$p)
	)


	# Is the target transcript a key driver?
	df = merge(
		df,
		kda[c("target_tissue", "target_clust", "target_gene_symbol", "key_driver_FDR")],
		by=c("target_tissue", "target_clust", "target_gene_symbol"),
		all.x=TRUE
	)

	# Order table by correlation strength
	# df = df[order(df$p), ]

	return(df)
})


# Test if target gene symbols agree
# so that the correlations can be collected in matrix
stopifnot(
	all(
		# test if all rows are the same
		apply(
			# matrix of target gene symbol
			sapply(endo_target_cor, function(x) x$target_gene_symbol),
			1,
			function(x) length(unique(x)) == 1  # only one entry
		)
	)
)


# Write table of sorted target genes
endo_target_cor_sorted = lapply(endo_target_cor, function(df) {
	df = df[order(df$p), ]
})

endo_target_cor_sorted_all = Reduce(rbind, endo_target_cor_sorted)
write.csv(endo_target_cor_sorted_all,
	"co-expression/module78to98/tables/endocrine_targets_78to98.csv",
	row.names=FALSE)


# Heatmap
cmat = sapply(endo_target_cor, function(x) x$cor)
colnames(cmat) = endo$gene_symbol
rownames(cmat) = endo_target_cor[[1]]$target_gene_symbol  # assumes that all targets are in the same order

key_driver_col = rep("white", nrow(cmat))

# key_driver_col[endo_target_cor[[1]]$key_driver_FDR < 0.05] = "grey"
key_driver_col[endo_target_cor[[1]]$key_driver_FDR < 0.05] = brewer.pal(9, "Set1")[3]

pdf("co-expression/module78to98/plots/target_gene_heatmap.pdf",
	height=5,
	width=8
)

idx = apply(cmat, 1, function(x) max(abs(x))) > 0.2
heatmap.2(
	# t(cmat),
	# ColSideColors=key_driver_col,
	ColSideColors=key_driver_col[idx],
	t(cmat[idx, ]),
	trace="none",
	Rowv=FALSE,
	col=colorRampPalette(rev(brewer.pal(9, "PuOr")))(64),
	breaks=seq(-0.3, 0.3, length.out=64 + 1),  # cap of coloring 
	xlab="Target genes in module 98",
	ylab="Endocrine factors",
	key.title="",
	key.xlab="cor.",
	density.info="none",
	mar=c(8, 6)
)
dev.off()


# Grouped barplot
m = 20  # number of top genes to show
row_idx = order(apply(cmat, 1, function(x) max(abs(x))), decreasing=TRUE)[1:m]
col_idx = match(c("FSTL3", "LBP", "STC2", "LEP"), colnames(cmat))
# col_idx = match(c("FSTL3", "LBP", "STC2", "FCN2", "LEP"), colnames(cmat))

cmat_sub = cmat[row_idx, col_idx]
cmat_sub = cmat_sub[order(apply(cmat_sub, 1, max), decreasing=TRUE), ]

key_drivers  = endo_target_cor[[1]]$target_gene_symbol[endo_target_cor[[1]]$key_driver_FDR < 0.05]
key_drivers = na.omit(key_drivers)

# Add key driver asterix annotation
rownames(cmat_sub) = sapply(rownames(cmat_sub), function(symbol) {
	if (symbol %in% key_drivers) {
		return(paste0(symbol, "*"))
	} else {
		return(symbol)
	}
})

pdf("co-expression/module78to98/plots/target_gene_grouped_barplot.pdf",
	height=3.5)

bar_cols = c(brewer.pal(9, "Set1")[1:3], "white")
barplot(
	t(cmat_sub),
	beside=TRUE,
	las=2,
	col=bar_cols,
	legend=colnames(cmat_sub),
	ylab="Cor."
)
dev.off()



endo_target_cor_all = Reduce(rbind, endo_target_cor)
endo_target_cor_all = endo_target_cor_all[order(endo_target_cor_all$p), ]

