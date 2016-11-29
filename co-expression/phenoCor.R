rm(list=ls())

# Analysis of modules
library(reshape2)
library(RColorBrewer)
library(gplots)
library(magicaxis)
library(devtools)

# heatmap.3
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")


data_dir = "/Users/sk/DataProjects/cross-tissue"  # root of data directory
setwd("/Users/sk/Google Drive/projects/cross-tissue")

source("co-expression/base.R")

# Load gwnet and row_meta tables
load(file.path(data_dir, "modules/cross-tissue.RData"))

# tissue_col = brewer.pal(9, "Pastel1")  # tissue colors
tissue_col = brewer.pal(9, "Set1")  # tissue colors
modules = as.integer(factor(bwnet$colors))  # the modules detected


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


# Module eigengenes, test if sample names are correct

# Load STARNET phenotype data
pheno = fread(file.path(
	"/Volumes/SANDY/phenotype_data",
	"STARNET_main_phenotype_table.cases.Feb_29_2016.tbl"
))

brainshake = fread(file.path(
	"/Volumes/SANDY/phenotype_data",
	"tz.mat"
))

# Match phenotype to sample order
pheno_matched = pheno[match(patient_ids, pheno$starnet.ID), ]

# calculate all eigengene-phenotype correlations
pheno_cor = cor(bwnet$MEs, pheno_matched, use="pairwise.complete.obs")
pheno_cor = as.data.frame(pheno_cor)

# exclude all NA columns
pheno_cor = pheno_cor[, 
	apply(pheno_cor, 2, function(col) any(!is.na(col)))]

# Exclude particular columns
exclude_pheno_cols = c("starnet.ID")
pheno_cor = pheno_cor[, !colnames(pheno_cor) %in% exclude_pheno_cols]



# Fit linear models of eigengenes vs phenotype.
fits = list()
for (pfeat in colnames(pheno_matched)[2:ncol(pheno_matched)]) {
	print(pfeat)

	# Make target~input data frame
	df = cbind(pheno_matched[[pfeat]], bwnet$MEs)
	if (length(levels(factor(df[,1]))) == 2) {
		# {-1, 1} dummy encoding for binomial categorical variables
		df[,1] = as.integer(factor(df[,1])) * 2 - 3
	}

	# exclude samples with missing target variable
	df = df[!is.na(df[,1]), ]
	df = data.frame(df)  # ensure same length? (don't know why it is necessary)

	try({
		fits[[pfeat]] = lm(
			# first vs rest
			as.formula(
				paste(colnames(df)[1], "~", 
					paste(colnames(df)[2:ncol(df)], collapse="+")
				)
			),
			data=df
		)
	})
}


# Get matrix of coefficient p values
alpha = 0.001  # minimum p-value to consider clinical feature and eigengene
# alpha = 0.01

pmat = sapply(fits, function(fit) {
	summary(fit)$coef[,4]
})

pmat = pmat[2:nrow(pmat),]  # remove intercept
rownames(pmat) = 1:nrow(pmat)  # module#

include_module = apply(pmat, 1, min) < alpha
include_phenotype = apply(pmat, 2, min) < alpha
pmat = pmat[include_module, include_phenotype]

module_tissue_freq[,include_module]


rlab = freqTissueColor(module_tissue_freq[,include_module], tissue_col)
rownames(rlab) = rownames(module_tissue_freq)

sub_freq = module_tissue_freq[,include_module]
colnames(sub_freq) = names(which(include_module))

pdf("co-expression/plots/pheno-eigengene.pdf")
heatmap.3(
	-log10(t(pmat)),
	trace="none",
	col=colorRampPalette(brewer.pal(9, "Blues"))(100),
	breaks=seq(0, 5, length.out=101),  # cap of coloring 
	ColSideColors=t(rlab),
	# labCol=rownames(module_tissue_freq),
	ColSideColorsSize=5,
	margins=c(6, 12),
	keysize=0.9,
	KeyValueName=expression("-log"[10] * " p"),
	xlab="Module eigengene"
)
dev.off()


# Visualize phenotype eigengene correlations
alpha=0.5
heatmap.2(
	t(as.matrix(pheno_cor)),
	trace="none",
	col=colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
	breaks=seq(-alpha, alpha, length.out=101)  # cap of coloring 
)


# BRAINSHAKE correlations
brainshake_matched = brainshake[match(patient_ids, brainshake$id), ]
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




# Sanity check
# -------------------------------------
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

