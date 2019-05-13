# adopted from endocrineValidationHMPD.R

rm(list=ls())

library(data.table)
library(WGCNA)
library(RColorBrewer)
library(gplots)

setwd("~/GoogleDrive/projects/STARNET/cross-tissue")

source("src/permuteTest.R")
source("src/parse.R")

data_dir = "/Users/sk/DataProjects/cross-tissue"  # root of data directory


# Load and parse HMDP gene expression data
hmdp = lapply(list.files(file.path(data_dir, "HMDP/chow_HF")),
	function(file_name) {
	file_path = file.path(data_dir, "HMDP/chow_HF", file_name)
	d = fread(file_path)
	d = data.frame(d)


	rownames(d) = d[, 1]  # mouse strain names
	d = data.matrix(d[, -1])

	genes = read.table(file_path, nrows=1)
	genes = t(genes)[, 1]
	genes = unname(genes)

	colnames(d) = genes

	return(d)
})
names(hmdp) = list.files(file.path(data_dir, "HMDP/chow_HF"))

# sapply(hmdp, dim)


# Select all of FDR threshold endocrine candidates
# endocrine_sel = "all"
endocrine_sel = "FDR"

# Load endocrine table
if (endocrine_sel == "FDR") {
	endo = fread("co-expression/tables/CT_endocrines_TS_interactions.csv")
} else if (endocrine_sel == "all") {
	endo = fread("co-expression/tables/CT_endocrines_TS_interactions_all.csv")
} else {
	stop("Wrong selector: ", endocrine_sel)
}


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

endo$mouse_symbol = findMouseHomologue(endo$gene_symbol, human_homology, mouse_homology)


endo[endo$clust == 78 & endo$target_clust == 98, ]


# Endocrine-eigenegene correlations
# ----------------------------------------------------

# mod_tab is a dataframe of module (clust) assignments to gene symbols, used for computing eigengenes for target_mat.
endoEigenCor = function(target_mat, mod_tab, source_mat, endocrines) {

	stopifnot(c("gene_symbol", "clust") %in% colnames(mod_tab))
	stopifnot(rownames(target_mat) == rownames(source_mat))

	# Get assigned clusters of target expression matrix
	idx = match(colnames(target_mat), mod_tab$gene_symbol)
	clust = mod_tab$clust[idx]

	# Calculate eigengenes
	# idx = match(colnames(target_mat), symbols)
	message("Calculating eigengenes...")
	eig = moduleEigengenes(target_mat, clust)
	colnames(eig$eigengenes) = substring(colnames(eig$eigengenes), 3)  # strip ME from eigengene names

	message("Estimating correlations...")
	# endocrine correlation tests and annotation of endocrine table.
	endo_eigen_cor = corAndPvalue(
		source_mat[, match(endocrines, colnames(source_mat))],
		eig$eigengenes)

	return(endo_eigen_cor)
}


# Adds correlation statistics to endocrine table
annotateEndoEigenCor = function(endo, cor_stats, prefix) {
	# Annotate endocrine table
	row_idx = match(endo$mouse_symbol, rownames(cor_stats$p))
	col_idx = match(endo$target_clust, colnames(cor_stats$p))

	# endo$HMDP_chow_ts_endo_cor = cor_stats$cor[cbind(row_idx, col_idx)]
	endo[[paste0(prefix, "_cor")]] = cor_stats$cor[cbind(row_idx, col_idx)]
	endo[[paste0(prefix, "_p")]] = cor_stats$p[cbind(row_idx, col_idx)]

	return(endo)
}


# Calculate eigengenes for all liver modules
idx = modules$tissue == "LIV"
modules_LIV_mouse = data.frame(
	gene_symbol=modules$mouse_symbol[idx],
	clust=modules$clust[idx])

endo78_98 = endo[endo$clust == 78 & endo$target_clust == 98, ]
endocrines = endo78_98$mouse_symbol
# # adipose -> liver
# adipose_liver_endocrines = endo$mouse_symbol[(endo$tissue == "SF" | endo$tissue == "VAF") & endo$target_tissue == "LIV"]
# adipose_liver_endocrines = as.character(na.omit(adipose_liver_endocrines))
# adipose_liver_endocrines = unique(adipose_liver_endocrines)


# idx = match(colnames(hmdp$Liver_HF_male), modules_LIV$mouse_symbol)

lapply(hmdp, dim)

adipose_liver_HF = endoEigenCor(
	hmdp$Liver_HF_male,
	modules_LIV_mouse,
	hmdp$Adipose_HF_male,
	# adipose_liver_endocrines
	endocrines
)

endo78_98 = annotateEndoEigenCor(endo78_98, adipose_liver_HF, "HMDP_HF")

adipose_liver_chow = endoEigenCor(
	target_mat=hmdp$Liver_chow_male,
	mod_tab=modules_LIV_mouse,
	source_mat=hmdp$Adipose_chow_male,
	endocrines=endocrines
)

endo78_98 = annotateEndoEigenCor(endo78_98, adipose_liver_chow, "HMDP_chow")


# Write table
endo[idx_78_98, ]
write.csv(
	endo[idx_78_98,
		c("clust", "target_clust", "tissue", "gene_symbol", "mouse_symbol",
			"ts_endo_cor", "ts_endo_cor_p", "ts_endo_cor_p_adj",
			"HMDP_HF_ts_endo_cor", "HMDP_HF_ts_endo_cor_p",
			"HMDP_chow_ts_endo_cor", "HMDP_chow_ts_endo_cor_p")
	],
	paste0("co-expression/plots/endocrine/endocrines_module_78_98_validation_heatmap_", endocrine_sel, ".csv"),
	row.names=FALSE
)



library(devtools)  # for installing heatmap.3
# Load heatmap.3
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")


# Correlation matrix of endocrine candidates from modules 78->98

pdf(paste0("co-expression/plots/endocrine/endocrines_module_78_98_validation_heatmap_", endocrine_sel, ".pdf"))
# pdf(paste0("co-expression/plots/endocrine/endocrines_module_78_98_validation_heatmap_", endocrine_sel, ".pdf"))
# idx_78_98 = endo$clust == 78 & endo$target_clust == 98

# endocrine_sel
# cmat_78_98 = endo[idx_78_98, c("ts_endo_cor", "HMDP_HF_ts_endo_cor", "HMDP_chow_ts_endo_cor")]
cmat_78_98 = endo78_98[, c("ts_endo_cor", "HMDP_HF_cor", "HMDP_chow_cor")]
rownames(cmat_78_98) = endo$id[idx_78_98]

colnames(cmat_78_98) = c("STARNET", "HMDP HF", "HMDP chow")


# colors = brewer.pal(9, "Set1")[c(8, 4)]
colors = brewer.pal(9, "Set1")[c(5, 8, 4)]

rlab = colors[as.integer(factor(endo$tissue[idx_78_98]))]
rlab = as.matrix(rlab)
colnames(rlab) = "Tissue of origin"

heatmap.3(data.matrix(cmat_78_98), trace="none",
	col=colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
	Rowv=FALSE,
	Colv=FALSE,
	mar=c(18, 34),
	# cexRow=1.0,
	cexRow=0.4,
	cexCol=1.0,
	RowSideColors=t(rlab),
	RowSideColorsSize = 0.7,
	KeyValueName="Module 98 eigengene cor."
)
dev.off()
