#!/usr/bin/env Rscript

# Preparation scripts, invoked from local machine
# Transfer commands to Minerva:

# STARNET network matrices
# rsync -avz --progress ~/DataProjects/cross-tissue/modules/matrices_STARNET koples01@minerva.hpc.mssm.edu:/sc/arion/projects/STARNET/koples01/module-validation/data/
# GTEx network matrices
# rsync -avz --progress ~/DataProjects/cross-tissue/modules/matrices_GTEx koples01@minerva.hpc.mssm.edu:/sc/arion/projects/STARNET/koples01/module-validation/data/

# Co-expression module data
# rsync -avz --progress ~/GoogleDrive/projects/STARNET/cross-tissue/co-expression/tables/modules.csv minerva:/sc/arion/projects/STARNET/koples01/module-validation/data/
# rsync -avz --progress ~/DataProjects/cross-tissue/modules/between_within-cross-tissue.RData minerva:/sc/arion/projects/STARNET/koples01/module-validation/data/


# Expression data

# STARNET
# rsync -avz --progress ~/DataProjects/cross-tissue/STARNET/gene_exp_norm_reshape/expr_recast.RData minerva:/sc/arion/projects/STARNET/koples01/module-validation/data/


# GTEx
# rsync -avz --progress ~/DataBases/GTEx/RNA-seq/ minerva:/sc/arion/projects/STARNET/koples01/module-validation/data/


# Transfer scripts
# Libraries
# rsync -avz --progress ~/GoogleDrive/projects/STARNET/cross-tissue/src/parse.R minerva:/sc/arion/projects/STARNET/koples01/module-validation/lib/
# rsync -avz --progress ~/GoogleDrive/projects/STARNET/cross-tissue/validateModules/lib/core.R minerva:/sc/arion/projects/STARNET/koples01/module-validation/lib/


# Transfer main scripts
# rsync -avz --progress ~/GoogleDrive/projects/STARNET/cross-tissue/validateModules/src minerva:/sc/arion/projects/STARNET/koples01/module-validation/

# Transfer results back from Minerva
# rsync -avz --progress minerva:/sc/arion/projects/STARNET/koples01/module-validation/results ~/GoogleDrive/projects/STARNET/cross-tissue/validateModules/

rm(list=ls())

library(NetRep)
library(data.table)


# Check command line argument input
# ------------------------------
args = commandArgs(trailingOnly=TRUE)

block = args[1]

if (length(args) != 1) {
	stop("Please provide block as command-line argument. USAGE: networkReplication.R <block_n>")
}


# Config
# -----------------------------------------

nPerm = 1000  # number of permutations of NetRep 
# nPerm = 50  # testing
nThreads = 4  # Number of threads for NetRep, 2 threads per core


# # Local setup
# # ------------------------------------
# wd = "~/GoogleDrive/projects/STARNET/cross-tissue"

# parse_script = "src/parse.R"
# core_script = "validateModules/lib/core.R"


# # Directories for network matrices
# starnet_dir = "/Users/sk/DataProjects/cross-tissue/modules/matrices_STARNET"
# gtex_dir = "/Users/sk/DataProjects/cross-tissue/modules/matrices_GTEx"
# temp_dir = "/Users/sk/DataProjects/cross-tissue/modules/matrices_temp"

# modules_path = "co-expression/tables/modules.csv"

# network_results_path = "/Users/sk/DataProjects/cross-tissue/modules/between_within-cross-tissue.RData"

# expr_path = "/Users/sk/DataProjects/cross-tissue/STARNET/gene_exp_norm_reshape/expr_recast.RData"

# gtex_expr_path = "~/DataBases/GTEx/RNA-seq"


# Minerva setup
# --------------------------------

wd = "/sc/arion/projects/STARNET/koples01/module-validation"

parse_script = "lib/parse.R"
core_script = "lib/core.R"


# Directories for network matrices
starnet_dir = "data/matrices_STARNET"
gtex_dir = "data/matrices_GTEx"
temp_dir = "tmp"

modules_path = "data/modules.csv"

network_results_path = "data/between_within-cross-tissue.RData"

expr_path = "data/expr_recast.RData"

gtex_expr_path = "data"



setwd(wd)
source(parse_script)
source(core_script)

# Convert list of matrices to disk matrices
diskMatrix = function(mat_list, dir, file_base) {
	dir.create(dir)
	mat_pointers = list()
	for (i in 1:length(mat_list)) {
		mat_pointers[[i]] = as.disk.matrix(
			x=mat_list[[i]],
			file=file.path(dir, paste0(file_base, "_", names(mat_list)[i], ".rds")),
			serialize=TRUE
		)
	}
	names(mat_pointers) = names(mat_list)
	return(mat_pointers)
}

# Load STARNET modules
modules = fread(modules_path)
modules$gene_biotype = parseEnsemblBiotype(modules$ensembl)
modules$tissue_biotype = paste(modules$tissue, modules$gene_biotype, sep="_")
modules$ensembl_base = sapply(strsplit(modules$ensembl, "[.]"), function(x) x[1])
modules$tissue_transcript_id = paste(modules$tissue, modules$ensembl_base, sep="_")
modules$tissue_symbol_transcript_id = paste(modules$tissue, modules$transcript_id, sep="_")


between_within = new.env()
load(network_results_path,
	between_within,
	verbose=TRUE)

gene_blocks = between_within$bwnet$gene_blocks


# Load and parse STARNET expression data
# ---------------------------------
load(expr_path, verbose=TRUE)

expr = parseExprTable(expr_recast)
names(expr)

rm(expr_recast)
gc()


# Standardize expression data
expr$mat = scale(t(expr$mat))
expr$mat[is.na(expr$mat)] = 0.0  # impute to average


# Load GTEx gene expression data, math to STARNET transcripts
# -------------------------------------
gtex = loadGTEx(dir=gtex_expr_path)
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


# Specify block
# ---------------------------------


# results = lapply(unique(gene_blocks), function(block) {

message("Block: ", block)

gc()

idx = gene_blocks == block

module_assign = list()
module_assign$starnet = modules$clust[idx]
names(module_assign$starnet) = modules$tissue_symbol_transcript_id[idx]


data_list = list(starnet=expr$mat[, idx], gtex=gtex$mat[, idx])


# Remove GTEx NA names, for missing transcript IDs
# colnames(data_list$gtex)[is.na(colnames(data_list$gtex))] = ""

# Transfer transcript IDs, assumes that dataset is matched (which it is)
colnames(data_list$gtex) = colnames(data_list$starnet)


correlation_list = list(
	starnet=readRDS(paste0(starnet_dir, "/cmat_block", block, ".rds")),
	gtex=readRDS(paste0(gtex_dir, "/cmat_block", block, ".rds")))

network_list = list(
	starnet=readRDS(paste0(starnet_dir, "/tom_block", block, ".rds")),
	gtex=readRDS(paste0(gtex_dir, "/tom_block", block, ".rds")))


# Convert TOM similarity to network weights
network_list$starnet = 1 - network_list$starnet
network_list$gtex = 1 - network_list$gtex

# replace missing correlation coefficients with zeros
correlation_list$gtex[is.na(correlation_list$gtex)] = 0

# Replace missing network weights with zeros
network_list$gtex[is.na(network_list$gtex)] = 0


# Transfer gene ids from STARNET
colnames(correlation_list$gtex) = colnames(correlation_list$starnet)
rownames(correlation_list$gtex) = rownames(correlation_list$starnet)

# Use row and column names for TOM matrices from correlation matrices
rownames(network_list$starnet) = rownames(correlation_list$starnet)
colnames(network_list$starnet) = colnames(correlation_list$starnet)

# # GTEx matrix column and row names: NA -> ""
# rownames(correlation_list$gtex)[is.na(rownames(correlation_list$gtex))] = ""
# colnames(correlation_list$gtex)[is.na(colnames(correlation_list$gtex))] = ""

# 
rownames(network_list$gtex) = rownames(correlation_list$gtex)
colnames(network_list$gtex) = colnames(correlation_list$gtex)


# Store prepareted data on disk to increase memory efficiency
correlation_list = diskMatrix(correlation_list, dir=temp_dir, file_base=paste("cmat_block_", block))
gc()
network_list = diskMatrix(network_list, dir=temp_dir, file_base=paste0("tom_block", block))
gc()

preservation = modulePreservation(
	network=network_list,
	data=data_list,
	correlation=correlation_list, 
	moduleAssignments=module_assign,
	discovery="starnet", test="gtex", 
	nPerm=nPerm,
	nThreads=nThreads
)

	# return(preservation)
# })

# names(results) = unique(gene_blocks)

dir.create("results")

saveRDS(preservation, file=paste0("results/network_preservation_block", block, ".rds"))
