rm(list=ls())

library(NetRep)
# library(RColorBrewer)
# library(viridis)
# library(pheatmap)
library(data.table)

setwd("~/GoogleDrive/projects/STARNET/cross-tissue")

source("src/parse.R")
source("validateModules/lib/core.R")


# Load STARNET modules
modules = fread("co-expression/tables/modules.csv")
modules$gene_biotype = parseEnsemblBiotype(modules$ensembl)
modules$tissue_biotype = paste(modules$tissue, modules$gene_biotype, sep="_")
modules$ensembl_base = sapply(strsplit(modules$ensembl, "[.]"), function(x) x[1])
modules$tissue_transcript_id = paste(modules$tissue, modules$ensembl_base, sep="_")
modules$tissue_symbol_transcript_id = paste(modules$tissue, modules$transcript_id, sep="_")


between_within = new.env()
load("/Users/sk/DataProjects/cross-tissue/modules/between_within-cross-tissue.RData",
	between_within,
	verbose=TRUE)

gene_blocks = between_within$bwnet$gene_blocks


# Load and parse STARNET expression data
# ---------------------------------
load("/Users/sk/DataProjects/cross-tissue/STARNET/gene_exp_norm_reshape/expr_recast.RData", verbose=TRUE)

expr = parseExprTable(expr_recast)
names(expr)

rm(expr_recast)
gc()


# Standardize expression data
expr$mat = scale(t(expr$mat))
expr$mat[is.na(expr$mat)] = 0.0  # impute to average


# Load GTEx gene expression data, math to STARNET transcripts
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


# Specify block
# ---------------------------------

# gene_blocks[modules$clust == 154]
# gene_blocks[modules$clust == 98]

# block = 1
# block = 3
# block = 5


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

results = lapply(unique(gene_blocks), function(block) {
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


	# Directories for network matrices
	starnet_dir = "/Users/sk/DataProjects/cross-tissue/modules/matrices_STARNET"
	gtex_dir = "/Users/sk/DataProjects/cross-tissue/modules/matrices_GTEx"
	temp_dir = "/Users/sk/DataProjects/cross-tissue/modules/matrices_temp"

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
	correlation_list = diskMatrix(correlation_list, dir=temp_dir, file_base="cmat")
	gc()
	network_list = diskMatrix(network_list, dir=temp_dir, file_base="tom")
	gc()

	preservation = modulePreservation(
		network=network_list,
		data=data_list,
		correlation=correlation_list, 
		moduleAssignments=module_assign,
		discovery="starnet", test="gtex", 
		# nPerm=10000,
		nPerm=1000,
		# nPerm=100,  # for testing
		nThreads=6
	)

	return(preservation)

})

names(results) = unique(gene_blocks)

# discovery_tom = readRDS(paste0(starnet_dir, "/tom_block", block, ".rds"))
# discovery_cmat = readRDS(paste0(starnet_dir, "/cmat_block", block, ".rds"))


# idx = gene_blocks == block
# table(modules$clust[idx])

# # # Block 1 examples
# # module = 189
# # module = 203
# # module = 152
# # module = 62
# module = 154  # Jason's study


# # # Block 3 examples
# # module = 98

# idx_block = modules$clust[idx] == module
# sum(idx_block)


# # discovery_tom[idx_block, idx_block]

# # mat = discovery_cmat[idx_block, idx_block]
# mat = discovery_tom[idx_block, idx_block]

# mat[mat == 0] = NA
# # mat[lower.tri(mat)] = NA

# mat = 1 - mat

# modules$tissue_gene_symbol = paste0(modules$tissue, "_", modules$gene_symbol)

# transcript_ids = modules$tissue_gene_symbol[idx][idx_block]
# gene_symbols = modules$gene_symbol[idx][idx_block]

# # gene_highlight = c("PAN2", "PHACTR1", "SMG6", "THOC5", "MRAS")

# labels = rep("", nrow(mat))

# show_idx = gene_symbols %in% gene_highlight

# labels[show_idx] = transcript_ids[show_idx]  # selected genes
# # labels = transcript_ids  # all



# rownames(mat) = labels
# colnames(mat) = labels

# sort(apply(mat, 2, mean, na.rm=TRUE), decreasing=TRUE)[1:50]

# ncolor_grad = 100
# pheatmap(
# 	mat,
# 	col=rev(inferno(ncolor_grad)),
# 	fontsize_row=6,
# 	fontsize_col=6
# 	)


# # discovery_tom = attach.disk.matrix(paste0(starnet_dir, "/tom_block", block, ".rds"))
# # as.matrix(attach.disk.matrix)[1:10, 1:10]

