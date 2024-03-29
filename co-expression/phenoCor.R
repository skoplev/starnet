rm(list=ls())

# Analysis of modules
library(reshape2)
library(RColorBrewer)
library(gplots)
library(magicaxis)
library(devtools)
library(qvalue)
library(VennDiagram)
library(data.table)

# heatmap.3
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")

data_dir = "/Users/sk/DataProjects/cross-tissue"  # root of data directory
setwd("/Users/sk/GoogleDrive/projects/STARNET/cross-tissue")

source("src/base.R")
source("src/parse/io.R")
source("src/models/regr.R")  # regression models




# tissue_col = brewer.pal(9, "Pastel1")  # tissue colors
tissue_col = brewer.pal(9, "Set1")  # tissue colors

# False discovery rates for linear regression coefficients
fdr = 0.2
alpha = 0.001  # nomical minimum p-value to consider clinical feature and eigengene


# Load data
# -------------------------------------------------------------
# Load gwnet and row_meta tables
load(file.path(data_dir, "modules/cross-tissue.RData"))
modules = as.integer(factor(bwnet$colors))  # the module assignment for each gene-tissue pair
# Rename modules to numbers
# colnames(bwnet$MEs) = 1:ncol(bwnet$MEs)


# Loads expr_recast data frame with tissue-specific expression
load(file.path(data_dir, "STARNET/gene_exp_norm_reshape/expr_recast.RData"))

# Get expression matrix from recasted data frame with all tissue 
mat = expr_recast[, 3:ncol(expr_recast)]

# Missing data rate
mean(is.na(mat))


# Count tissue composition of each module
module_tissue_comp = lapply(1:length(unique(modules)), function(i) {
	return(table(row_meta$tissue[modules == i]))
})

module_tissue = countMat(module_tissue_comp)

# Calculate tissue frequencies
module_tissue_freq = sweep(module_tissue, 2, apply(module_tissue, 2, sum), "/")

# Load STARNET phenotype data
pheno = fread(file.path(
	"/Volumes/SANDY/phenotype_data",
	"STARNET_main_phenotype_table.cases.Feb_29_2016.tbl"
))

# Load Brainshake phenotype data
brainshake = fread(file.path(
	"/Volumes/SANDY/phenotype_data",
	"tz.mat"
))

# Load CIBERSORT frequency data
cibersort_freq = loadCibersortFreq(file.path(data_dir, "CIBERSORT/out_freq"))

# Match phenotype to sample order
pheno_matched = pheno[match(patient_ids, pheno$starnet.ID), ]
brainshake_matched = brainshake[match(patient_ids, brainshake$id), ]

cibersort_freq_matched = matchTrimCibersortFreq(cibersort_freq, patient_ids)


# Fit linear models, eigengenes->phenotype
# --------------------------------------------------

fits_pheno = fitLinearEigenPheno(pheno_matched, bwnet$MEs)
# summary(fits_pheno[["Age"]])


syntax_coef = summary(fits_pheno$syntax_score)$coef
# rownames(syntax_coef)
syntax_coef = syntax_coef[syntax_coef[,4] < 0.1, ]


# Get matrix of coefficient p values
pmat = getPmat(fits_pheno)
rownames(pmat) = 1:nrow(pmat)  # module# rename

# q-values for control of FDR
qmat = qvalue(pmat)$lfdr

# mat = filterMatRowColMin(qmat, fdr)
alpha = 0.05
mat = filterMatRowColMin(pmat, alpha)

# Make colors for tissue frequencies
include_module_pheno = as.integer(rownames(mat))
rlab = freqTissueColor(module_tissue_freq[,include_module_pheno], tissue_col)
rownames(rlab) = rownames(module_tissue_freq)

# Plot heatmap of significance
pdf("co-expression/plots/pheno-eigengene.pdf")
heatmap.3(
	-log10(t(mat)),
	trace="none",
	col=colorRampPalette(brewer.pal(9, "Blues"))(100),
	breaks=seq(0, 4, length.out=101),  # cap of coloring 
	ColSideColors=t(rlab),
	ColSideColorsSize=5,
	margins=c(6, 12),
	keysize=0.9,
	KeyValueName=expression("-log"[10] * " p"),
	xlab="Module eigengene"
)
dev.off()


# Brainshake regression model
# -----------------------------------------------------------------
fits_brainshake = fitLinearEigenPheno(brainshake_matched, bwnet$MEs)

pmat_brainshake = getPmat(fits_brainshake)
rownames(pmat_brainshake) = 1:nrow(pmat_brainshake)  # rename modules to numbers

qmat_brainshake = qvalue(pmat_brainshake)$lfdr

mat = filterMatRowColMin(pmat_brainshake, alpha)
# mat = filterMatRowColMin(qmat_brainshake, fdr)

# Make tissue frequency colors
include_module_brainshake = as.integer(rownames(mat))
rlab = freqTissueColor(module_tissue_freq[, include_module_brainshake], tissue_col)
rownames(rlab) = rownames(module_tissue_freq)

pdf("co-expression/plots/brainshake-eigengene.pdf", height=10)
heatmap.3(
	-log10(t(mat)),
	trace="none",
	col=colorRampPalette(brewer.pal(9, "Blues"))(100),
	breaks=seq(0, 4, length.out=101),  # cap of coloring 
	ColSideColors=t(rlab),
	ColSideColorsSize=3,
	margins=c(6, 12),
	keysize=0.9,
	KeyValueName=expression("-log"[10] * " p"),
	xlab="Module eigengene"
)
dev.off()


# Venn diagram of module associations
# ----------------------------------------------------------
# Disable log for Venn diagram plots
futile.logger::flog.threshold(
	futile.logger::ERROR,
	name = "VennDiagramLogger")

# Make Venn diagram of annotated modules
venn = venn.diagram(
	list(
		Phenotype=include_module_pheno,
		Brainshake=include_module_brainshake),
	fill=brewer.pal(9, "Set1")[1:2],
	alpha=0.15,
	lwd=0.5,
	cex=2,
	main="Module associations",
	filename=NULL)

pdf("co-expression/plots/module_annot_venn.pdf")
grid.draw(venn)
dev.off()


# Correlation analysis
# ------------------------------------------------------------------


# calculate all eigengene-phenotype correlations
pheno_cor = cor(bwnet$MEs, pheno_matched, use="pairwise.complete.obs")
pheno_cor = as.data.frame(pheno_cor)

# exclude all NA columns
pheno_cor = pheno_cor[, 
	apply(pheno_cor, 2, function(col) any(!is.na(col)))]

# Exclude particular columns
exclude_pheno_cols = c("starnet.ID")
pheno_cor = pheno_cor[, !colnames(pheno_cor) %in% exclude_pheno_cols]

# Visualize phenotype eigengene correlations
alpha=0.5
heatmap.2(
	t(as.matrix(pheno_cor)),
	trace="none",
	col=colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
	breaks=seq(-alpha, alpha, length.out=101)  # cap of coloring 
)


# BRAINSHAKE correlations
brainshake_cor = cor(bwnet$MEs, brainshake_matched, use="pairwise.complete.obs")
brainshake_cor = as.data.frame(brainshake_cor)
brainshake_cor = brainshake_cor[, 
	apply(brainshake_cor, 2, function(col) any(!is.na(col)))]

alpha=0.5
heatmap.2(
	t(as.matrix(brainshake_cor)),
	trace="none",
	col=colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
	breaks=seq(-alpha, alpha, length.out=101),  # cap of coloring 
	cexRow=0.3
)



# Sanity checks
# ----------------------------------------------------------------
mod_i = which.min(pheno_cor$HDL)

plot(bwnet$MEs[,mod_i], pheno_matched$HDL)

mod_idx = which(abs(pheno_cor$HDL) > 0.3)
which(abs(pheno_cor$BMI) > 0.3)

which(abs(pheno_cor) > 0.3)

module_tissue[, mod_idx]

barplot(module_tissue[, mod_idx], col=tissue_col,
	)
legend("topright", legend=rownames(module_tissue),
	pch=22,
	pt.bg=tissue_col
	# cex=0.7
)


row_meta[modules == mod_i,]

x = as.matrix(mat[modules == mod_i,])
x = scale(t(x))
x[is.na(x)] = 0.0

alpha = 2
heatmap.2(x,
	trace="none",
	col=colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
	breaks=seq(-alpha, alpha, length.out=101)  # cap of coloring 
)



# Test if samples of eigen genes and correlation matrix agrees.
# Runtime ~5 min
testSampleCor = function() {
	# Loads expr_recast data frame with tissue-specific expression
	load(file.path(data_dir, "STARNET/gene_exp_norm_reshape/expr_recast.RData"))

	eigen_genes = bwnet$MEs
	rownames(eigen_genes) = patient_ids

	eigen_cor = cor(t(eigen_genes))

	# Loads expr_recast data frame with tissue-specific expression
	load(file.path(data_dir, "STARNET/gene_exp_norm_reshape/expr_recast.RData"))

	# Get expression matrix from recasted data frame with all tissue 
	mat = expr_recast[, 3:ncol(expr_recast)]
	row_meta = expr_recast[, 1:2]

	colnames(mat) = colnames(expr_recast)[3:ncol(expr_recast)]  # STARNET patient name IDs

	expr_cor = cor(mat, use="pairwise.complete.obs")

	plot(expr_cor[lower.tri(expr_cor)], eigen_cor[lower.tri(eigen_cor)])

	cor_test = cor.test(expr_cor[lower.tri(expr_cor)], eigen_cor[lower.tri(eigen_cor)],
		use="pairwise.complete.obs")
	return(cor_test)
}

