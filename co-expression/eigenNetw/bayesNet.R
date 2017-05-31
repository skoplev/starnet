# Bayesian network model of module eigengene including CIBERSORT frequencies
# and phenotype features.

options(java.parameters = "-Xmx8g")  # Max Java memory heap size, for rcausal

library(rcausal)
library(data.table)
library(igraph)

data_dir = "/Users/sk/DataProjects/cross-tissue"  # root of data directory
setwd("/Users/sk/Google Drive/projects/cross-tissue")

source("src/parse.R")


# Load cross-tissues modules
between = new.env()
load(file.path(data_dir, "modules/between_within-cross-tissue.RData"),
	between,
	verbose=TRUE)

between = parseModuleData(between)


# Load phenotype data
# ---------------------------
# STARNET phenotype data
pheno = fread(file.path(
	"/Volumes/SANDY/phenotype_data",
	"STARNET_main_phenotype_table.cases.Feb_29_2016.tbl"
))

# Match phenotype data tables to
pheno_matched = pheno[match(between$patient_ids, pheno$starnet.ID), ]


# Load module overview table
module_tab = read.table("co-expression/tables/module_tab.csv",
	sep=",",
	header=TRUE
)


# CIBERSORT frequencies, association with eigengenes
# ---------------------------------------------------------------
# Load all files
freq_files = list.files(file.path(data_dir, "CIBERSORT/out_freq"))

exclude_freq_files = c(
	"exp.GRCh38.GENCODE_r24.COR.mat.tsv",
	"exp.GRCh38.GENCODE_r24.FOC.mat.tsv",
	"exp.GRCh38.GENCODE_r24.MAC.mat.tsv")

freq_files = freq_files[!freq_files %in% exclude_freq_files]

# Parse files returning only matrices
ciber_freq = parseCibersortFiles(freq_files, data_dir)

# Align CIBERSORT frequency data to patient_ids
ciber_freq = lapply(ciber_freq, function(freq_mat) {
	idx = match(
		between$patient_ids,
		sapply(strsplit(rownames(freq_mat), "_"), function(x) x[2])
		)

	freq_mat = freq_mat[idx, ]

	return(freq_mat)
})

# Test if all rownames agrees with patient IDs (patient_ids), throws error if not.
void = sapply(ciber_freq, function(freq_mat) {
	freq_ids = sapply(strsplit(rownames(freq_mat), "_"), function(x) x[2])
	if(!all(freq_ids == between$patient_ids, na.rm=TRUE)) {
		stop("Patient ID mismatch")
	}
})

# Combine frequencies across tissue into single matrix
ciber_freq_mat = Reduce(cbind, ciber_freq)
rownames(ciber_freq_mat) = between$patient_ids


# Collect features for Bayesian network
sub_ciber_freq_mat = ciber_freq_mat[,
	c(
		"AOR:blood:macrophage",
		"AOR:blood:dendritic cell, myeloid, immature",
		"AOR:blood:dendritic cell, plasmacytoid",
		"AOR:blood:monocyte",
		"AOR:blood:natural killer cell",
		"AOR:blood:T cell, CD4+",
		"AOR:blood:B cell",
		"AOR:blood:neutrophil",
		"AOR:heart:preadipocyte",
		"AOR:adipose tissue:preadipocyte",
		"AOR:aorta:fibroblast",
		"AOR:omentum:preadipocyte"
	), drop=FALSE]


phenotypes = c("syntax_score", "ndv", "lesions", "DUKE", "BMI", "HbA1c", "LDL", "HDL", "CRP", "p_chol", "bl.glucose")
# phenotypes = c("DUKE", "BMI", "HbA1c", "LDL", "HDL", "CRP", "p_chol", "bl.glucose")
sub_pheno_matched = as.data.frame(pheno_matched)[, phenotypes]

# Select modules
include_modules = which(module_tab$cad_pval < 0.05 | module_tab$pval_DUKE < 0.05)
length(include_modules)

variables = cbind(
	between$bwnet$eigengenes[, include_modules],
	sub_ciber_freq_mat,
	sub_pheno_matched
)

bn = fges(variables, maxDegree=100)  # Fast-greedy equivalence search

g = graph_from_graphnel(bn$graphNEL)

V(g)$log_P_DUKE = -log10(1.0)
V(g)$log_P_DUKE[1:length(include_modules)] = -log10(module_tab$pval_DUKE[include_modules])

V(g)$log_P_cad = -log10(1.0)
V(g)$log_P_cad[1:length(include_modules)] = -log10(module_tab$cad_pval[include_modules])

V(g)$type = c(
	rep("module", length(include_modules)),
	rep("cell_type_freq", ncol(sub_ciber_freq_mat)),
	rep("phenotype", ncol(sub_pheno_matched))
)	

write_graph(g,
	file="co-expression/eigenNetw/v1/bayes_net2.gml",
	format="gml")
