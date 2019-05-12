# Permutation test of STARNET liver modules in HMDP ApoE data (atherosclerosis)
# Adpoted from endocrineValidationHMDP.R

rm(list=ls())

library(data.table)
library(WGCNA)
library(RColorBrewer)
library(gplots)

setwd("~/Google Drive/projects/STARNET/cross-tissue")

source("src/permuteTest.R")

data_dir = "/Users/sk/DataProjects/cross-tissue"  # root of data directory

# Load data, returning a numeric matrix -- genes in column, mice IDs in columns
# fread warning is okay: "Starting data input on line 2 and discarding line 1 because".
parseArrayData = function(file_path) {
	# load
	d = fread(file_path)

	# drop first column, which is repeated gene_symbol
	stopifnot(all(d[, 1] == d[, 2]))
	d = d[, -1]

	# Get header of file
	header = read.table(file_path, nrows=1, as.is=TRUE)[1, ]
	header = as.character(header)

	colnames(d) = header

	# Numeric matrix
	mat = data.matrix(d[, -1])
	rownames(mat) = d$gene_symbol
	return(t(mat))
}


# Load HMDP ApoE data
hmdp = list()
hmdp$female = parseArrayData(file.path(data_dir, "HMDP/ath_liver/HMDP Female Ath Liver Expression.txt"))
hmdp$male = parseArrayData(file.path(data_dir, "HMDP/ath_liver/HMDP Male Ath Liver Expression.txt"))

# all(hmdp$female == hmdp$male)


# Load STARNET modules
modules = fread("co-expression/tables/modules.csv")


# Load mouse homology and map STARNET gene symbols to mouse gene symbols
human_mouse_homology = fread(file.path(data_dir, "MGI/HOM_MouseHumanSequence.rpt"))

# Separate table into human and mouse
human_homology = human_mouse_homology[
	human_mouse_homology[["Common Organism Name"]] == "human", ]

mouse_homology = human_mouse_homology[
	human_mouse_homology[["Common Organism Name"]] == "mouse, laboratory", ]


modules$mouse_symbol = findMouseHomologue(modules$gene_symbol, human_homology, mouse_homology)






# Permutation test
# -------------------------------------------------

# Calculate all correlation matrices, for mouse genes with homologues in STARNET
cor_mats = lapply(hmdp, matchCor, unique(modules$mouse_symbol))

# all(cor_mats[[1]] == cor_mats[[2]])  # ? same data.

m = 10000

# Module 98 genes
mod = 98

# mod = 20
# mod = 5
# mod = 28

genes = modules$mouse_symbol[modules$clust == mod]
genes = na.omit(genes)

perm_tests = lapply(cor_mats, function(cmat) {
	corPermTest(cmat, genes, m)
})

pdf(paste0("co-expression/plots/endocrine/validationHMDP/perm_tests/HMDP_ApoE_mod", mod, ".pdf"),
	width=4, height=6)
par(mfrow=c(2, 1))
plotPermuteTest(perm_tests$female, main=paste(mod, "HMDP ApoE Leiden female"))
plotPermuteTest(perm_tests$male, main=paste(mod, "HMDP ApoE Leiden male"))
dev.off()
