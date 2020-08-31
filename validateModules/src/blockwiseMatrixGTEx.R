rm(list=ls())

library(data.table)

library(compiler)
enableJIT(3)

library(WGCNA)
# enableWGCNAThreads(nThreads=2)  # assuming 2 threads per core


setwd("~/GoogleDrive/projects/STARNET/cross-tissue")

source("src/parse.R")
source("validateModules/lib/core.R")



# Load STARNET modules
modules = fread("co-expression/tables/modules.csv")
modules$gene_biotype = parseEnsemblBiotype(modules$ensembl)
modules$tissue_biotype = paste(modules$tissue, modules$gene_biotype, sep="_")
modules$ensembl_base = sapply(strsplit(modules$ensembl, "[.]"), function(x) x[1])
modules$tissue_transcript_id = paste(modules$tissue, modules$ensembl_base, sep="_")


# Load GTEx data
# -------------------------------------
gtex = loadGTEx()
gc()

# Match by GTEx transcripts by tissue transcript ID
idx = match(modules$tissue_transcript_id, gtex$meta_row$tissue_transcript_id)

sum(!is.na(idx))

gtex$meta_row = gtex$meta_row[idx, ]
gtex$mat = gtex$mat[idx, ]


# Standardize expression data
gtex$mat = scale(t(gtex$mat))
gtex$mat[is.na(gtex$mat)] = 0.0  # impute to average
gc()


# Load gene blocks used for WGCNA
# ------------------------------------
between_within = new.env()
load("/Users/sk/DataProjects/cross-tissue/modules/between_within-cross-tissue.RData",
	between_within,
	verbose=TRUE)

# Test if genes are identifical
stopifnot(all(between_within$meta_genes$transcript_id == modules$transcript_id))

gene_blocks = between_within$bwnet$gene_blocks


beta_mat = betaMat()

stopifnot(all(unique(na.omit(gtex$meta_row$tissue)) %in% colnames(beta_mat)))


gc()

computeBlockwiseCmatTOM(gtex,
	gene_blocks,
	beta_mat,
	out_dir="/Users/sk/DataProjects/cross-tissue/modules/matrices_GTEx")
