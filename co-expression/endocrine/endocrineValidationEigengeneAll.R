# adopted from endocrineValidationHMPD.R

rm(list=ls())

library(data.table)
library(WGCNA)
library(RColorBrewer)
library(gplots)
library(GEOquery)

setwd("~/GoogleDrive/projects/STARNET/cross-tissue")

source("src/permuteTest.R")
source("src/parse.R")

data_dir = "/Users/sk/DataProjects/cross-tissue"  # root of data directory


# Load gene expression data
# ---------------------------------------------

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


# Load GTEx
gtex = list()
gtex$mat = fread("~/DataBases/GTEx/RNA-seq/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct")
gtex$genes = gtex$mat[, 1:2]
gtex$mat = data.matrix(gtex$mat[, c(-1, -2)])
rownames(gtex$mat) = make.names(gtex$genes$Description, unique=TRUE)

gtex$annot = fread("~/DataBases/GTEx/RNA-seq/GTEx_v7_Annotations_SampleAttributesDS.txt")

# Match annotation table to column names
# note that first two columns of mat are dummy IDs
gtex$annot = gtex$annot[match(colnames(gtex$mat), gtex$annot$SAMPID), ]

gtex$annot$gtex.ID = sapply(strsplit(gtex$annot$SAMPID, "-"), function(x) paste(x[1], x[2], sep="-"))


# Extract and match gene tissue-specific gene expression matrices
getGtex = function(idx) {
	mat = gtex$mat[, idx]
	colnames(mat) = gtex$annot$gtex.ID[idx]
	return(mat)
}

gtex$LIV_mat = getGtex(gtex$annot$SMTSD == "Liver")
gtex$SF_mat = getGtex(gtex$annot$SMTSD == "Adipose - Subcutaneous")
gtex$VAF_mat = getGtex(gtex$annot$SMTSD == "Adipose - Visceral (Omentum)")

# Match to LIV
gtex$SF_mat = gtex$SF_mat[, match(colnames(gtex$LIV_mat), colnames(gtex$SF_mat))]
gtex$VAF_mat = gtex$VAF_mat[, match(colnames(gtex$LIV_mat), colnames(gtex$VAF_mat))]


# Morbid obesity gene expression
# ------------------------------------------------------
# morbid = loadMorbidObesity()

# save(morbid, file="~/DataBases/MorbidObesity/morbid.RData")
load(file="~/DataBases/MorbidObesity/morbid.RData")

rownames(morbid$emat) = morbid$meta_row$hgnc_symbol
morbid$emat = t(morbid$emat)
rownames(morbid$emat) = morbid$sample_annot$MGH.ID


# Get submatrices
morbid$emat_LIV = morbid$emat[morbid$sample_annot$tissue_name == "liver", ]
morbid$emat_SF = morbid$emat[morbid$sample_annot$tissue_name == "subcutaneous adipose", ]
morbid$emat_VAF = morbid$emat[morbid$sample_annot$tissue_name == "omental adipose", ]


# Match sample IDs to those of liver
morbid$emat_SF = morbid$emat_SF[match(rownames(morbid$emat_LIV), rownames(morbid$emat_SF)), ]
morbid$emat_VAF = morbid$emat_VAF[match(rownames(morbid$emat_LIV), rownames(morbid$emat_VAF)), ]

stopifnot(rownames(morbid$emat_LIV) == rownames(morbid$emat_SF))
stopifnot(rownames(morbid$emat_LIV) == rownames(morbid$emat_VAF))


# Select all of FDR threshold endocrine candidates
# ------------------------------------------------------
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
# ------------------------------------------------
human_mouse_homology = fread(file.path(data_dir, "MGI/HOM_MouseHumanSequence.rpt"))

# Separate table into human and mouse
human_homology = human_mouse_homology[
	human_mouse_homology[["Common Organism Name"]] == "human", ]

mouse_homology = human_mouse_homology[
	human_mouse_homology[["Common Organism Name"]] == "mouse, laboratory", ]


modules$mouse_symbol = findMouseHomologue(modules$gene_symbol, human_homology, mouse_homology)

endo$mouse_symbol = findMouseHomologue(endo$gene_symbol, human_homology, mouse_homology)

# Based on old gene symbol (@Marcus)
endo$mouse_symbol[endo$mouse_symbol == "Fcn2"] = "Fcna"


endo[endo$clust == 78 & endo$target_clust == 98, ]


# Endocrine-eigenegene correlations
# ----------------------------------------------------

# mod_tab is a dataframe of module (clust) assignments to gene symbols, used for computing eigengenes for target_mat.
endoEigenCor = function(target_mat, mod_tab, source_mat, endocrines) {

	stopifnot(c("gene_symbol", "clust") %in% colnames(mod_tab))
	stopifnot(na.exclude(rownames(target_mat) == rownames(source_mat)))

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
annotateEndoEigenCor = function(endo, cor_stats, prefix, gene_column) {
	# Annotate endocrine table
	row_idx = match(endo[[gene_column]], rownames(cor_stats$p))
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
endocrines_mouse = endo78_98$mouse_symbol
# endocrines_mouse = endo78_98$mouse_symbol

endocrines_human = endo78_98$gene_symbol

idx = modules$tissue == "LIV"
modules_LIV_human = data.frame(
	gene_symbol=modules$gene_symbol[idx],
	clust=modules$clust[idx])


# HMDP
# --------------------------------------------

# # adipose -> liver
# adipose_liver_endocrines = endo$mouse_symbol[(endo$tissue == "SF" | endo$tissue == "VAF") & endo$target_tissue == "LIV"]
# adipose_liver_endocrines = as.character(na.omit(adipose_liver_endocrines))
# adipose_liver_endocrines = unique(adipose_liver_endocrines)

adipose_liver_HF = endoEigenCor(
	hmdp$Liver_HF_male,
	modules_LIV_mouse,
	hmdp$Adipose_HF_male,
	endocrines_mouse
)

endo78_98 = annotateEndoEigenCor(endo78_98, adipose_liver_HF, "HMDP_HF", "mouse_symbol")

adipose_liver_chow = endoEigenCor(
	target_mat=hmdp$Liver_chow_male,
	mod_tab=modules_LIV_mouse,
	source_mat=hmdp$Adipose_chow_male,
	endocrines=endocrines_mouse
)

endo78_98 = annotateEndoEigenCor(endo78_98, adipose_liver_chow, "HMDP_chow", "mouse_symbol")


# GTEx
# ------------------------------------------

VAF_liver_gtex = endoEigenCor(
	# Use subset of gene expression matrix
	target_mat=t(gtex$LIV_mat[rownames(gtex$LIV_mat) %in% modules_LIV_human$gene_symbol, ]),
	mod_tab=modules_LIV_human,
	source_mat=t(gtex$VAF_mat),
	endocrines=endocrines_human
)
endo78_98 = annotateEndoEigenCor(endo78_98, VAF_liver_gtex, "GTEx_VAF", "gene_symbol")

SF_liver_gtex = endoEigenCor(
	# Use subset of gene expression matrix
	target_mat=t(gtex$LIV_mat[rownames(gtex$LIV_mat) %in% modules_LIV_human$gene_symbol, ]),
	mod_tab=modules_LIV_human,
	source_mat=t(gtex$SF_mat),
	endocrines=endocrines_human
)

endo78_98 = annotateEndoEigenCor(endo78_98, SF_liver_gtex, "GTEx_SF", "gene_symbol")


# Morbid obesity
# ---------------------------

VAF_liver_morbid = endoEigenCor(
	target_mat=morbid$emat_LIV[, colnames(morbid$emat_LIV) %in% modules_LIV_human$gene_symbol],
	mod_tab=modules_LIV_human,
	source_mat=morbid$emat_VAF,
	endocrines=endocrines_human
)
endo78_98 = annotateEndoEigenCor(endo78_98, VAF_liver_morbid, "Obese_VAF", "gene_symbol")

SF_liver_morbid = endoEigenCor(
	target_mat=morbid$emat_LIV[, colnames(morbid$emat_LIV) %in% modules_LIV_human$gene_symbol],
	mod_tab=modules_LIV_human,
	source_mat=morbid$emat_SF,
	endocrines=endocrines_human
)
endo78_98 = annotateEndoEigenCor(endo78_98, SF_liver_morbid, "Obese_SF", "gene_symbol")



# # Write table
# endo[idx_78_98, ]
# write.csv(
# 	endo[idx_78_98,
# 		c("clust", "target_clust", "tissue", "gene_symbol", "mouse_symbol",
# 			"ts_endo_cor", "ts_endo_cor_p", "ts_endo_cor_p_adj",
# 			"HMDP_HF_ts_endo_cor", "HMDP_HF_ts_endo_cor_p",
# 			"HMDP_chow_ts_endo_cor", "HMDP_chow_ts_endo_cor_p")
# 	],
# 	paste0("co-expression/plots/endocrine/endocrines_module_78_98_validation_heatmap_", endocrine_sel, ".csv"),
# 	row.names=FALSE
# )


# STARNET endocrine candidate correlation with clinical traits
# ------------------------------------------------

library(devtools)  # for installing heatmap.3
# Load heatmap.3
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")


# Correlation matrix of endocrine candidates from modules 78->98

pdf(paste0("co-expression/plots/endocrine/endocrines_module_78_98_validation_heatmap_", endocrine_sel, ".pdf"))
# pdf(paste0("co-expression/plots/endocrine/endocrines_module_78_98_validation_heatmap_", endocrine_sel, ".pdf"))
# idx_78_98 = endo$clust == 78 & endo$target_clust == 98

# endocrine_sel
features = c("ts_endo", "HMDP_HF", "HMDP_chow", "GTEx_VAF", "GTEx_SF", "Obese_VAF", "Obese_SF")

# for systematically extracting p-values
colnames(endo78_98)[colnames(endo78_98) == "ts_endo_cor_p"] = "ts_endo_p"

cmat_78_98 = data.frame(endo78_98)[, paste0(features, "_cor")]
pmat_78_98 = data.frame(endo78_98)[, paste0(features, "_p")]

# cmat_78_98 = endo78_98[, c("ts_endo_cor", "HMDP_HF_cor", "HMDP_chow_cor", "GTEx_VAF_cor", "GTEx_SF_cor", "Obese_VAF_cor", "Obese_SF_cor")]
rownames(cmat_78_98) = endo78_98$id

# colnames(cmat_78_98) = c("STARNET", "HMDP HF", "HMDP chow")

tissues = c("AOR", "BLOOD", "SKLM", "VAF", "MAM", "LIV",  "SF")
tissue_cols = brewer.pal(9, "Set1")[-6]

rlab = tissue_cols[as.integer(factor(endo78_98$tissue, levels=tissues))]
rlab = as.matrix(rlab)
colnames(rlab) = "Tissue of origin"

heatmap.3(data.matrix(cmat_78_98),
	trace="none",
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
