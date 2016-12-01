rm(list=ls())

# Analysis of modules
library(reshape2)
library(RColorBrewer)
library(gplots)
library(magicaxis)
library(devtools)
library(qvalue)
library(VennDiagram)

# heatmap.3
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")

data_dir = "/Users/sk/DataProjects/cross-tissue"  # root of data directory
setwd("/Users/sk/Google Drive/projects/cross-tissue")

source("src/base.R")
source("src/parse/io.R")


# Fits multivariate linear model for each phenotype based the eigengenes
# Input matrices are required to be sample matched.
# Samples in rows, features in columns
fitLinearEigenPheno = function(
	phenotype,
	module_eigengenes,
	exclude_phenotypes=c("starnet.ID", "id"))
{
	phenotype = as.data.frame(phenotype)  # if data.table
	module_eigengenes = as.data.frame(module_eigengenes)

	stopifnot(nrow(phenotype) == nrow(module_eigengenes))

	# Exclude phenotype columns
	phenotype = phenotype[, !colnames(phenotype) %in% exclude_phenotypes]

	# Fit linear models of eigengenes vs phenotype.
	fits = list()
	for (pfeat in colnames(phenotype)) {
		message("Fitting ", pfeat)

		# Make target~input data frame
		df = cbind(phenotype[[pfeat]], module_eigengenes)
		if (length(levels(factor(df[,1]))) == 2) {
			# {-1, 1} dummy encoding for binomial categorical variables
			df[,1] = as.integer(factor(df[,1])) * 2 - 3
		}

		# exclude samples with missing target variable
		df = df[!is.na(df[,1]), ]
		df = data.frame(df)  # ensure same length? (don't know why it is necessary)

		tryCatch({
			# Linear fit formula comparing
			first_rest_form = as.formula(
				paste(colnames(df)[1], "~", 
					paste(colnames(df)[2:ncol(df)], collapse="+")
				)
			)

			# first vs rest linear fit
			fits[[pfeat]] = lm(first_rest_form, data=df)
		}, error=function(e) {
			# Warning for failed linear fits
			warning(pfeat, " excluded.")
		})
	}

	return(fits)
}

# Get pmat matrix from list of linear fits
getPmat = function(fits, remove_intercept=TRUE) {
	# number of parameters
	pars = sapply(fits, function(fit) {
		nrow(summary(fit)$coef)
	})

	# exclude fits without the correct number of coefficients
	# 
	warnings("Excluding ", names(which(pars != median(pars))))
	fits = fits[pars == median(pars)]

	pmat = sapply(fits, function(fit) {
		summary(fit)$coef[,4]
	})

	if (remove_intercept) {
		pmat = pmat[!rownames(pmat) %in% "(Intercept)", ]
	}
	return(pmat)
}


# filter matrix such that it contains at least one row or column lower than alpha
# Used for p-value matrix filtering.
filterMatRowColMin = function(mat, alpha) {
	include_row = apply(mat, 1, min) < alpha
	include_col = apply(mat, 2, min) < alpha
	return(mat[include_row, include_col])
}


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

# Get matrix of coefficient p values
pmat = getPmat(fits_pheno)
rownames(pmat) = 1:nrow(pmat)  # module# rename

# q-values for control of FDR
qmat = qvalue(pmat)$lfdr

# mat = filterMatRowColMin(qmat, fdr)
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

