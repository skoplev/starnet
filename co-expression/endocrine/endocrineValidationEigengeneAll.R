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

sapply(hmdp, dim)


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

# Calculate eigengenes for all liver modules
tissue = "LIV"
# idx = modules$tissue == tissue

modules_LIV = modules[modules$tissue == "LIV", ]

liver_modules_all = list()
for (k in unique(modules$clust)) {
	liver_modules_all[[k]] = modules$mouse_symbol[modules$clust == k & modules$tissue == tissue]
}


liver_modules_all = lapply(liver_modules_all, function(x) as.character(na.omit(x)))





# Calculate eigengenes
idx = match(colnames(hmdp$Liver_HF_male), modules_LIV$mouse_symbol)
eigen_liver_HF = moduleEigengenes(hmdp$Liver_HF_male, modules_LIV$clust[idx])
colnames(eigen_liver_HF$eigengenes) = substring(colnames(eigen_liver_HF$eigengenes), 3)  # strip ME from eigengene names


idx = match(colnames(hmdp$Liver_chow_male), modules_LIV$mouse_symbol)
eigen_liver_chow = moduleEigengenes(hmdp$Liver_chow_male, modules_LIV$clust[idx])
colnames(eigen_liver_chow$eigengenes) = substring(colnames(eigen_liver_chow$eigengenes), 3)  # strip ME from eigengene names



# Check if rows are aligned, for correlations between HMDP tissues
if (!all(rownames(hmdp$Liver_HF_male) == rownames(hmdp$Adipose_HF_male))) {
	stop("Liver_HF_male is not matched to Adipose_HF_male")
}

if (!all(rownames(hmdp$Liver_chow_male) == rownames(hmdp$Adipose_chow_male))) {
	stop("Liver_chow_male is not matched to Adipose_chow_male")
}


adipose_liver_endocrines = endo$mouse_symbol[(endo$tissue == "SF" | endo$tissue == "VAF") & endo$target_tissue == "LIV"]
adipose_liver_endocrines = as.character(na.omit(adipose_liver_endocrines))
adipose_liver_endocrines = unique(adipose_liver_endocrines)

# Annotate endocrine table: VAF/SF -> LIV
idx = (endo$tissue == "VAF" | endo$tissue == "SF") & endo$target_tissue == "LIV"


# HF, endocrine correlation tests and annotation of endocrine table.
endocrine_eigen_cor_adipose_liv_HF = corAndPvalue(
	hmdp$Adipose_HF_male[, match(adipose_liver_endocrines, colnames(hmdp$Adipose_HF_male))],
	eigen_liver_HF$eigengenes)

row_idx = match(endo$mouse_symbol[idx], rownames(endocrine_eigen_cor_adipose_liv_HF$p))
col_idx = match(endo$target_clust[idx], colnames(endocrine_eigen_cor_adipose_liv_HF$p))

endo$HMDP_HF_ts_endo_cor[idx] = endocrine_eigen_cor_adipose_liv_HF$cor[cbind(row_idx, col_idx)]
endo$HMDP_HF_ts_endo_cor_p[idx] = endocrine_eigen_cor_adipose_liv_HF$p[cbind(row_idx, col_idx)]


# chow, endocrine correlation
endocrine_eigen_cor_adipose_liv_chow = corAndPvalue(
	hmdp$Adipose_chow_male[, match(adipose_liver_endocrines, colnames(hmdp$Adipose_chow_male))],
	eigen_liver_chow$eigengenes)

# Annotate endocrine table
row_idx = match(endo$mouse_symbol[idx], rownames(endocrine_eigen_cor_adipose_liv_chow$p))
col_idx = match(endo$target_clust[idx], colnames(endocrine_eigen_cor_adipose_liv_chow$p))

endo$HMDP_chow_ts_endo_cor[idx] = endocrine_eigen_cor_adipose_liv_chow$cor[cbind(row_idx, col_idx)]
endo$HMDP_chow_ts_endo_cor_p[idx] = endocrine_eigen_cor_adipose_liv_chow$p[cbind(row_idx, col_idx)]


# Some counts
# sum(idx)
# sum(!is.na(endo$mouse_symbol[idx]))
# sum(!is.na(endo$HMDP_chow_ts_endo_cor_p[idx]))
# sum(!is.na(endo$HMDP_HF_ts_endo_cor_p[idx]))


# endo[endo$clust == 14, ]
# endo[endo$clust == 122, ]

# endo[endo$clust == 218, ]
# endo[endo$clust == 28, ]
# endo[endo$clust == 106, ]
# endo[endo$clust == 167, ]
# endo[endo$clust == 135, ]
# endo[endo$clust == 134, ]

# endo$HMDP_chow_ts_endo_cor_p < 0.05 | 

# sum(endo$HMDP_chow_ts_endo_cor_p < 0.05 | endo$HMDP_HF_ts_endo_cor_p < 0.05, na.rm=TRUE)



library(devtools)  # for installing heatmap.3
# Load heatmap.3
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")


# Correlation matrix of endocrine candidates from modules 78->98

pdf(paste0("co-expression/plots/endocrine/endocrines_module_78_98_validation_heatmap_", endocrine_sel, ".pdf"))
# pdf(paste0("co-expression/plots/endocrine/endocrines_module_78_98_validation_heatmap_", endocrine_sel, ".pdf"))
idx_78_98 = endo$clust == 78 & endo$target_clust == 98

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

# endocrine_sel

cmat_78_98 = endo[idx_78_98, c("ts_endo_cor", "HMDP_HF_ts_endo_cor", "HMDP_chow_ts_endo_cor")]
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
