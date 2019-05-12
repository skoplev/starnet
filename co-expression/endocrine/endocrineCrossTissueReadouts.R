rm(list=ls())


library(data.table)
library(WGCNA)

setwd("/Users/sk/Google Drive/projects/STARNET/cross-tissue")
data_dir = "/Users/sk/DataProjects/cross-tissue"  # root of data directory


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
endo = fread("co-expression/tables/CT_endocrines_TS_interactions_mouse.csv")

k = 1

# endo$gene_symbol[k]
# endo$target_clust[k]
# endo$tissue[k]


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
	df = df[order(df$p), ]

	return(df)
})


endo_target_cor = Reduce(rbind, endo_target_cor)

write.csv(endo_target_cor, "co-expression/tables/endocrines_target_correlations_mouse_validated.csv",
	row.names=FALSE)