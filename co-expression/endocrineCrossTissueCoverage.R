rm(list=ls())

library(WGCNA)
library(data.table)
library(RColorBrewer)
enableWGCNAThreads(nThreads=8)

data_dir = "/Users/sk/DataProjects/cross-tissue"  # root of data directory
setwd("/Users/sk/Google Drive/projects/cross-tissue")
source("src/parse.R")



# endocrine = fread("~/Google Drive/projects/STARNET-endocrine/data/endo_scores_large.csv")

endocrine = fread("~/Google Drive/projects/STARNET-endocrine/data/endo_scores.csv")







# Load cross-tissue modules
# ------------------------------
between = new.env()
load(file.path(data_dir, "modules/between_within-cross-tissue.RData"),
	between,
	verbose=TRUE)

# Parse module data
between = parseModuleData(between)


# Load module annotations
module_tab = fread("co-expression/tables/module_tab.csv")


# Find modules of each endocrine factor
# ----------------------------------
endocrine_ids = paste(endocrine$from_tissue, endocrine$endocrine_factor, sep=":")
transcript_ids = paste(between$meta_genes$tissue, between$meta_genes$gene_symbol, sep=":")

endocrine$clust = between$clust[match(endocrine_ids, transcript_ids)]

# Add module purity to endocrine table
endocrine$mod_purity = module_tab$purity[endocrine$clust]

# Cross-tissue?
endocrine$cross_tissue_10 = endocrine$mod_purity < 0.9
endocrine$cross_tissue_5 = endocrine$mod_purity < 0.95
endocrine$cross_tissue_1 = endocrine$mod_purity < 0.99

# Write example table for particular module
write.csv(endocrine[endocrine$clust == 98, ], "co-expression/plots/endocrines_mod98.csv")

sec_score_thresh = 2
# sec_score_thresh = 1.5

pdf("co-expression/plots/cross_tissue_endocrine.pdf", height=3.5, width=7)
par(mfcol=c(1, 2),
	mar=c(5.1, 4.1, 4.1, 2.1)  # default
)
hist(endocrine$sec_score,
	breaks=100,
	xlab="Endocrine secretion score",
	main="",
	col="black")
abline(v=2, col=brewer.pal(9, "Set1")[1], lty=3)
# sum(endocrine$sec_score > sec_score_thresh)

unique(endocrine$endocrine_factor)
unique(endocrine$endocrine_factor[endocrine$sec_score > sec_score_thresh])

mat = cbind(
	table(endocrine$cross_tissue_1[endocrine$sec_score > sec_score_thresh]),
	table(endocrine$cross_tissue_5[endocrine$sec_score > sec_score_thresh]),
	table(endocrine$cross_tissue_10[endocrine$sec_score > sec_score_thresh])
)
colnames(mat) = c("1", "5", "10")
mat = mat[2:1, ]  # flip rows


# Calculate cross-tissue percentages
mat[1,] / apply(mat, 2, sum)

par(mar=c(5.1, 4.1, 4.1, 8.1))
barplot(
	mat,
	ylab="Endocrines",
	xlab="Cross-tissue transcripts (%)"
)
dev.off()



# Load expression data
# -----------------------------------
emat_file = "STARNET/gene_exp_norm_reshape/expr_recast.RData"

load(file.path(data_dir, emat_file), verbose=TRUE)

emat = expr_recast[, 3:ncol(expr_recast)]
meta_genes = expr_recast[, 1:2]
meta_genes = as.data.frame(meta_genes)

rm(expr_recast)

# Test if the gene metadata is the same as for the cross-tissue modules
if (!all(meta_genes$transcript_id == between$meta_genes$transcript_id)) {
	stop("Transcript mismatch")
}


t1 = "VAF"
t2 = "LIV"

cmat = cor(t(emat[meta_genes$tissue == t1, ]), t(emat[meta_genes$tissue == t2, ]),
	use="pairwise.complete.obs")

