# Positive control for the permutation test

rm(list=ls())

library(data.table)
library(RColorBrewer)
library(WGCNA)
library(compiler)
enableJIT(3)

setwd("~/Google Drive/projects/STARNET/cross-tissue")

source("src/permuteTest.R")
source("src/parse.R")


data_dir = "/Users/sk/DataProjects/cross-tissue"  # root of data directory


# Load STARNET modules
modules = fread("co-expression/tables/modules.csv")
modules$tissue_transcript_id = paste(modules$tissue, modules$transcript_id, sep="_")
modules$gene_biotype = parseEnsemblBiotype(modules$ensembl)
modules$tissue_biotype = paste(modules$tissue, modules$gene_biotype, sep="_")

# Load STARNET expression data
load(file.path(data_dir, "STARNET/gene_exp_norm_reshape/expr_recast.RData"),
	verbose=TRUE)

starnet = parseExprTable(expr_recast)
rm(expr_recast)


# Cross-tissue correlation matrix for STARNET. VAF-LIV
# ------------------------------------------------------------------------


m = 1000  # permutation iterations

# Cross-tissue module 78
# Get adipose-liver data
liv_adipose = list()
idx = starnet$meta_row$tissue %in% c("VAF", "SF", "LIV")
liv_adipose$mat = starnet$mat[idx, ]
liv_adipose$mat = t(liv_adipose$mat)
liv_adipose$meta_col = starnet$meta_row[idx, ]

# Module genes
k = 78
mod_idx = modules$clust == k
genes = modules$tissue_transcript_id[mod_idx]
# group = modules$tissue[mod_idx]
group = modules$tissue_biotype[mod_idx]

perm_test = list()

perm_test[[k]] = corPermTestExprMat(
	expr_mat=liv_adipose$mat,
	genes=genes,
	m=m,
	# mat_group=liv_adipose$meta_col$tissue,
	mat_group=liv_adipose$meta_col$tissue_biotype,
	genes_group=group
)



# Module 98
liv = list()
idx = starnet$meta_row$tissue == "LIV"
liv$mat = starnet$mat[idx, ]
liv$mat = t(liv$mat)
liv$meta_col = starnet$meta_row[idx, ]

k = 98
mod_idx = modules$clust == k
genes = modules$tissue_transcript_id[mod_idx]
group = modules$tissue_biotype[mod_idx]

perm_test[[k]] = corPermTestExprMat(
	expr_mat=liv$mat,
	genes=genes,
	m=m,
	mat_group=liv$meta_col$tissue_biotype,
	genes_group=group
)

pdf("co-expression/plots/endocrine/validationSTARNET/perm_tests/perm_test_78_98.pdf", width=3, height=5)
par(mfrow=c(2, 1))
plotPermuteTest(perm_test[[78]], main="STARNET module 78")
plotPermuteTest(perm_test[[98]], main="STARNET module 98")
dev.off()

# perm_test = corPermTestExprMat(expr_mat=liv_adipose_mat, genes, m)
