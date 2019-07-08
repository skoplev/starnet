# adopted from endocrineValidationHMPD.R
# Open R session from terminal
# $ open /Applications/R.app
# with max heap size set above default, ie. from .bash_profile
# $ export R_MAX_VSIZE=32000000000

rm(list=ls())

library(data.table)
library(WGCNA)
library(RColorBrewer)
library(gplots)
library(GEOquery)
library(stringr)
library(plyr)
library(biomaRt)


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

# Load HMDP apoe-leiden data. Liver (microarray) and adipose (RNA-seq).
# ----------------------------------------------

# Load biomaRt database for annotating genes
ensembl = useEnsembl(biomart="ensembl",
	dataset="mmusculus_gene_ensembl")
	# version=89)

gene_map = getBM(
	attributes=c("ensembl_gene_id", "mgi_symbol", "gene_biotype"),
	mart=ensembl)


apoe = list()

# Adipose gene expression
apoe$adipose_mat = fread("data/apoe-leiden-mice/adipose_RNAseq/gene_counts.txt")
apoe$adipose_genes = apoe$adipose_mat[, 1:6]

# Add gene symbols
apoe$adipose_genes = cbind(apoe$adipose_genes,
	gene_map[match(apoe$adipose_genes$Geneid, gene_map$ensembl_gene_id), c("mgi_symbol", "gene_biotype")]
)

apoe$adipose_mat = data.matrix(apoe$adipose_mat[, -1:-6])
rownames(apoe$adipose_mat) = apoe$adipose_genes$mgi_symbol

# Pseudo-log transform
apoe$adipose_mat = log2(apoe$adipose_mat + 1)

# Filter out zero variance genes
apoe$adipose_mat = apoe$adipose_mat[apply(apoe$adipose_mat, 1, sd) > 0, ]

# Sample annotation table
samples = data.frame(
	str_match(
		colnames(apoe$adipose_mat),
		"^ms([0-9]+)"),
	stringsAsFactors=FALSE)
colnames(samples) = c("id_str", "id")

apoe$adipose_annot = fread("data/apoe-leiden-mice/adipose_RNAseq/adipose mouse IDs and labels.txt")
colnames(apoe$adipose_annot)[colnames(apoe$adipose_annot) == "Sample #"] = "id"

apoe$adipose_samples = join(samples, apoe$adipose_annot, by="id")

colnames(apoe$adipose_mat) = apoe$adipose_samples$mouse_number  # for matching samples
	

# liver gene expression
apoe$liver_female_mat = fread("data/apoe-leiden-mice/ath_liver/HMDP Female Ath Liver Expression.txt")
apoe$liver_male_mat = fread("data/apoe-leiden-mice/ath_liver/HMDP Male Ath Liver Expression.txt")

stopifnot(apoe$liver_female_mat$gene_symbol == apoe$liver_male_mat$gene_symbol)

# Combine
apoe$liver_mat = cbind(apoe$liver_female_mat, apoe$liver_male_mat[, c(-1, -2)])
gene_symbol = apoe$liver_mat$gene_symbol
apoe$liver_mat = data.matrix(apoe$liver_mat[, c(-1, -2)])
rownames(apoe$liver_mat) = gene_symbol


# Match to available adipose samples 
apoe$liver_mat_match = apoe$liver_mat[, match(colnames(apoe$adipose_mat), colnames(apoe$liver_mat))]
stopifnot(colnames(apoe$adipose_mat) == colnames(apoe$liver_mat_match))



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


# Load STARNET data 
# ------------------------------------
modules = fread("co-expression/tables/modules.csv")


# Expression matrix
load("~/DataProjects/cross-tissue/STARNET/gene_exp_norm_reshape/expr_recast.RData",
	verbose=TRUE)

starnet = parseExprTable(expr_recast)

rm(expr_recast)

starnet$pheno = fread(
	"~/GoogleDrive/projects/STARNET/phenotype/data/current/STARNET_main_phenotype_table.2017_12_03.tsv"
)

# Match phenotype data
starnet$pheno_match = starnet$pheno[match(colnames(starnet$mat), starnet$pheno$starnet.ID), ]


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

endo_sel = endo[endo$clust == 78 & endo$target_clust == 98, ]  # 78 -> 98 endocrines

endocrines_mouse = endo_sel$mouse_symbol
# endocrines_mouse = endo_sel$mouse_symbol

endocrines_human = endo_sel$gene_symbol

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

endo_sel = annotateEndoEigenCor(endo_sel, adipose_liver_HF, "HMDP_HF", "mouse_symbol")

adipose_liver_chow = endoEigenCor(
	target_mat=hmdp$Liver_chow_male,
	mod_tab=modules_LIV_mouse,
	source_mat=hmdp$Adipose_chow_male,
	endocrines=endocrines_mouse
)

endo_sel = annotateEndoEigenCor(endo_sel, adipose_liver_chow, "HMDP_chow", "mouse_symbol")



# HMDP ApoE-Leiden
# ------------------------------------------

adipose_liver_apoe = endoEigenCor(
	target_mat=t(apoe$liver_mat_match),
	mod_tab=modules_LIV_mouse,
	source_mat=t(apoe$adipose_mat),
	endocrines=endocrines_mouse
)

endo_sel = annotateEndoEigenCor(endo_sel, adipose_liver_apoe, "HMDP_apoe", "mouse_symbol")


# GTEx
# ------------------------------------------
VAF_liver_gtex = endoEigenCor(
	# Use subset of gene expression matrix
	target_mat=t(gtex$LIV_mat[rownames(gtex$LIV_mat) %in% modules_LIV_human$gene_symbol, ]),
	mod_tab=modules_LIV_human,
	source_mat=t(gtex$VAF_mat),
	endocrines=endocrines_human
)
endo_sel = annotateEndoEigenCor(endo_sel, VAF_liver_gtex, "GTEx_VAF", "gene_symbol")

SF_liver_gtex = endoEigenCor(
	# Use subset of gene expression matrix
	target_mat=t(gtex$LIV_mat[rownames(gtex$LIV_mat) %in% modules_LIV_human$gene_symbol, ]),
	mod_tab=modules_LIV_human,
	source_mat=t(gtex$SF_mat),
	endocrines=endocrines_human
)

endo_sel = annotateEndoEigenCor(endo_sel, SF_liver_gtex, "GTEx_SF", "gene_symbol")


# Morbid obesity
# ---------------------------
VAF_liver_morbid = endoEigenCor(
	target_mat=morbid$emat_LIV[, colnames(morbid$emat_LIV) %in% modules_LIV_human$gene_symbol],
	mod_tab=modules_LIV_human,
	source_mat=morbid$emat_VAF,
	endocrines=endocrines_human
)
endo_sel = annotateEndoEigenCor(endo_sel, VAF_liver_morbid, "Obese_VAF", "gene_symbol")

SF_liver_morbid = endoEigenCor(
	target_mat=morbid$emat_LIV[, colnames(morbid$emat_LIV) %in% modules_LIV_human$gene_symbol],
	mod_tab=modules_LIV_human,
	source_mat=morbid$emat_SF,
	endocrines=endocrines_human
)
endo_sel = annotateEndoEigenCor(endo_sel, SF_liver_morbid, "Obese_SF", "gene_symbol")



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

# TODO...

# starnet$pheno_match

features = c(
	"syntax_score",
	"DUKE",
	# "ndv",
	# "lesions",
	"BMI(kg/m2)",
	"CRP(mg/l)",
	"HbA1c(%)",
	# "Waist/Hip",
	"P-Chol(mmol/l)",
	"fP-LDL-Chol(mmol/l)",
	"fP-HDL-Chol(mmol/l)",
	"fP-TG(mmol/l)"
)

# features %in% colnames(starnet$pheno_match)
# starnet$pheno_match[, features]
# data.frame(starnet$pheno_match, check.names=FALSE)[, features]

# Match transcripts by id
starnet$meta_row$tissue_gene_symbol = paste(starnet$meta_row$tissue,
	starnet$meta_row$gene_symbol,
	sep="_")

endo_sel$tissue_gene_symbol = paste(endo_sel$tissue, endo_sel$gene_symbol, sep="_")

# Correlation statistics
idx = match(
	endo_sel$tissue_gene_symbol,
	starnet$meta_row$tissue_gene_symbol)

starnet_pheno = corAndPvalue(
	t(starnet$mat[idx, ]),
	data.frame(starnet$pheno_match, check.names=FALSE)[, features]
)


# Visualization
# ------------------------------------------------

# library(devtools)
# install_github("jokergoo/ComplexHeatmap")

library(ComplexHeatmap)
library(circlize)


pdf(paste0("co-expression/plots/endocrine/endocrines_module_78_98_validation_heatmap_v3_", endocrine_sel, ".pdf"),
	height=5,
	width=7.5)

# endocrine_sel
features = c("ts_endo", "HMDP_HF", "HMDP_chow", "GTEx_VAF", "GTEx_SF", "Obese_VAF", "Obese_SF", "HMDP_apoe")

# for systematically extracting p-values
colnames(endo_sel)[colnames(endo_sel) == "ts_endo_cor_p"] = "ts_endo_p"

paste0(features, "_cor") %in% colnames(endo_sel)

cmat = data.frame(endo_sel)[, paste0(features, "_cor")]
cmat = data.matrix(cmat)
rownames(cmat) = endo_sel$gene_symbol

pmat = data.frame(endo_sel)[, paste0(features, "_p")]

# Cell notes based on significance

sigNotes = function(pmat) {
	cellnote = matrix("", ncol=ncol(pmat), nrow=nrow(pmat))
	cellnote[pmat < 0.1] = "."
	cellnote[pmat < 0.05] = "*"
	cellnote[pmat < 0.01] = "**"
	cellnote[pmat < 0.001] = "***"
	return(cellnote)
}

cellnote1 = sigNotes(pmat)
cellnote2 = sigNotes(starnet_pheno$p)


# Coloring
tissues = c("AOR", "BLOOD", "SKLM", "VAF", "MAM", "LIV",  "SF")
tissue_cols = brewer.pal(8, "Set1")[-6]
names(tissue_cols) = tissues
na_col = rgb(220, 220, 220, maxColorValue=255)

ha = rowAnnotation(
# ha = HeatmapAnnotation(
	Tissue=endo_sel$tissue,
	# col=brewer.pal(9, "Set1")[-6]
	# col=as.list(tissue_cols)
	# annotation_width=2.01,
	# annotation_height=2.01,
	# simple_anno_size=0.5,
	border=TRUE,
	# width=0.2,
	col=list(Tissue=tissue_cols)
)

ht1 = Heatmap(
	cmat,
	name="Eigengene 98 cor. (LIV)",
	row_title="Adipose endocrine candidate",
	col=colorRamp2(
		seq(-0.3, 0.3, length.out=9),
		rev(brewer.pal(9, "RdBu"))),
	na_col=na_col,
	left_annotation=ha,
	row_names_side="left",
	border=TRUE,
	cell_fun=function(j, i, x, y, w, h, col) {
	    grid.text(cellnote1[i, j], x, y)
	}
)

ht2 = Heatmap(
	starnet_pheno$cor,
	name="Clinical cor.",
	col=colorRamp2(
		seq(-0.5, 0.5, length.out=9),
		rev(brewer.pal(9, "PuOr"))),
	na_col=na_col,
	border=TRUE,
	show_row_names=FALSE,
	cell_fun=function(j, i, x, y, w, h, col) {
	    grid.text(cellnote2[i, j], x, y)
	}
)

# Plot
ht1 + ht2

dev.off()



# OLD
# ---------------------------

library(devtools)  # for installing heatmap.3
# Load heatmap.3
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")






# Correlation matrix of endocrine candidates from modules 78->98

pdf(paste0("co-expression/plots/endocrine/endocrines_module_78_98_validation_heatmap_", endocrine_sel, ".pdf"))
# pdf(paste0("co-expression/plots/endocrine/endocrines_module_78_98_validation_heatmap_", endocrine_sel, ".pdf"))

# endocrine_sel
features = c("ts_endo", "HMDP_HF", "HMDP_chow", "GTEx_VAF", "GTEx_SF", "Obese_VAF", "Obese_SF")

# for systematically extracting p-values
colnames(endo_sel)[colnames(endo_sel) == "ts_endo_cor_p"] = "ts_endo_p"

cmat = data.frame(endo_sel)[, paste0(features, "_cor")]
pmat = data.frame(endo_sel)[, paste0(features, "_p")]


cellnote = matrix("", ncol=ncol(pmat), nrow=nrow(pmat))
cellnote[pmat < 0.1] = "."
cellnote[pmat < 0.05] = "*"
cellnote[pmat < 0.01] = "**"
cellnote[pmat < 0.001] = "***"

rownames(cmat) = endo_sel$id

tissues = c("AOR", "BLOOD", "SKLM", "VAF", "MAM", "LIV",  "SF")
tissue_cols = brewer.pal(9, "Set1")[-6]

rlab = tissue_cols[as.integer(factor(endo_sel$tissue, levels=tissues))]
rlab = as.matrix(rlab)
colnames(rlab) = "Tissue of origin"

heatmap.3(data.matrix(cmat),
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
	KeyValueName="Module 98 eigengene cor.",
	cellnote=cellnote,
	notecol="black"
)
dev.off()
