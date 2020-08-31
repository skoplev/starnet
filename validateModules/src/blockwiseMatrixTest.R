# For increased memory heapsize open R with:
# open -n /Applications/R.app

rm(list=ls())

library(data.table)

library(compiler)
enableJIT(3)

library(WGCNA)
# enableWGCNAThreads(nThreads=2)  # assuming 2 threads per core


setwd("~/GoogleDrive/projects/STARNET/cross-tissue")

source("src/parse.R")
source("validateModules/lib/core.R")


# Configuration
# ---------------------------
opts = list()


opts$data_dir = "/Users/sk/DataProjects/cross-tissue/STARNET"  # root of data directory
opts$emat_file = "gene_exp_norm_reshape/expr_recast.RData"

opts$beta_mat = betaMat()


# Load and parse expression data
# ---------------------------------
load(file.path(opts$data_dir, opts$emat_file), verbose=TRUE)

expr = parseExprTable(expr_recast)
names(expr)

rm(expr_recast)
gc()


# Standardize expression data
expr$mat = scale(t(expr$mat))
expr$mat[is.na(expr$mat)] = 0.0  # impute to average



# Load co-expression module data
# ------------------------------------
between_within = new.env()
load("/Users/sk/DataProjects/cross-tissue/modules/between_within-cross-tissue.RData",
	between_within,
	verbose=TRUE)

# Test if genes are identifical
stopifnot(all(between_within$meta_genes$transcript_id == expr$meta_row$transcript_id))

gene_blocks = between_within$bwnet$gene_blocks



gene_sel = 1:1000

table(gene_blocks[gene_sel])

expr_sub = list()
expr_sub$mat = expr$mat[, gene_sel]
expr_sub$meta_row = expr$meta_row[gene_sel, ]

test_dir = "/Users/sk/DataProjects/cross-tissue/modules/matrices_STARNET_test"
computeBlockwiseCmatTOM(
	expr_sub,
	gene_blocks[gene_sel],
	opts$beta_mat,
	out_dir=test_dir)


# Investigate output
block = 4

tom = readRDS(paste0(test_dir, "/tom_block", block, ".rds"))
cmat = readRDS(paste0(test_dir, "/cmat_block", block, ".rds"))

cmat[10:20, 10:20]
tom[10:20, 10:20]

