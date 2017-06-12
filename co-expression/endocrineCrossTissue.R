# 
rm(list=ls())

library(data.table)
library(WGCNA)
library(gplots)
library(RColorBrewer)

library(devtools)  # for installing heatmap.3
# Load heatmap.3
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")


data_dir = "/Users/sk/DataProjects/cross-tissue"  # root of data directory
setwd("/Users/sk/Google Drive/projects/cross-tissue")
source("src/parse.R")

addAlpha = function(col, alpha=1){
	if (missing(col)) stop("Please provide a vector of colours.")
	apply(sapply(col, col2rgb)/255, 2, 
		function(x) 
		rgb(x[1], x[2], x[3], alpha=alpha)
	)  
}

# Load cross-tissue modules
# ------------------------------
between = new.env()
load(file.path(data_dir, "modules/between_within-cross-tissue.RData"),
	between,
	verbose=TRUE)

# Parse module data
between = parseModuleData(between)

# Main module table
mod_tab = fread("co-expression/tables/module_tab.csv")

secreted_proteins = fread(file.path(data_dir, "Uniprot/uniprot_human_secreted_proteins.tab"))

sec_prots = secreted_proteins[["Gene names  (primary )"]]
sec_prots = strsplit(paste(sec_prots, collapse=";"), ";")[[1]]
sec_prots = trimws(sec_prots)
sec_prots = unique(sec_prots)
sec_prots = sec_prots[sec_prots != ""]


# Load expression data
# -----------------------------------
emat_file = "STARNET/gene_exp_norm_reshape/expr_recast.RData"

load(file.path(data_dir, emat_file), verbose=TRUE)

emat = expr_recast[, 3:ncol(expr_recast)]
meta_genes = expr_recast[, 1:2]
meta_genes = as.data.frame(meta_genes)

rm(expr_recast)

# Test if the gene metadata is the same as for the cross-tissue modules
if (!all(meta_genes$transcript_id == between$meta_genes$transcript_id)) {
	stop("Transcript mismatch")
}


# Secretory proteins in cross-tissue modules
# --------------------------------------
sum(sec_prots %in% between$meta_genes$gene_symbol)  # total number of mapped secretory protein symbols

# Find secretory proteins in cross-tissue modules
# length(sec_prots)

sec_idx = between$meta_genes$gene_symbol %in% sec_prots
sum(sec_idx)

cross_tissue_frac_defn = 0.05
# cross_tissue_frac_defn = 0.01

# cross_tissue_modules = which(mod_tab$purity < 0.99)
cross_tissue_modules = which(mod_tab$purity < 1 - cross_tissue_frac_defn)

sec_idx = between$meta_genes$gene_symbol %in% sec_prots & between$clust %in% cross_tissue_modules
sum(sec_idx)



# unique(between$meta_genes$gene_symbol[sec_idx])

table(between$meta_genes$tissue[sec_idx])

# Endocrine-eigengene correlations
cmat = cor(t(emat[sec_idx, ]), between$bwnet$eigengenes,
	use="pairwise.complete.obs")

rownames(cmat) = paste0(between$meta_genes$tissue, ":", between$meta_genes$gene_symbol)[sec_idx]

# Table for secreted proteins
sec_info = data.frame(
	tissue=between$meta_genes$tissue[sec_idx],
	gene_symbol=between$meta_genes$gene_symbol[sec_idx],
	clust=between$clust[sec_idx])
sec_info$id = paste0(sec_info$tissue, ":", sec_info$gene_symbol)

table(sec_info$tissue)

length(unique(sec_info$gene_symbol))


# plot(density(cmat))

# sum(cmat > 0.8)

# include = apply(abs(cmat), 1, max) > 0.7

# # sum(apply(abs(cmat), 1, max) > 0.5)
# mat = cmat[include, ]
# mat = cmat[include, mod_tab$CAD_qvalue < 0.1]

# mat = cmat[include, mod_tab$CAD_qvalue < 0.1 & mod_tab$purity < 0.99]  # Cross-tissue
# mat = cmat[include, mod_tab$CAD_qvalue < 0.1 & mod_tab$purity >= 0.99]  # Tissue-specific


tissues = unique(between$meta_genes$tissue)
tissue_col = brewer.pal(9, "Set1")[-6]

# mod_tab$LIV < 6

# cor_thresh = 0.25
cor_thresh = 0.20
# cor_thresh = 0.3

# cor_thresh = 0.5

# cor_thresh = 0.7

# t1 = "LIV"
# t1 = "VAF"
# t1 = "AOR"
# t1 = "MAM"
# t1 = "SKLM"
# t1 = "SF"
# t1 = "BLOOD"

for (t1 in tissues) {
	pdf(paste0("co-expression/plots/endocrine/eigengeneCor/", t1, ".pdf"), width=8)
	include_mod = mod_tab$CAD_qvalue < 0.1 | 
		mod_tab$pval_syntax_score < 0.05 |
		mod_tab$pval_ndv < 0.05 |
		mod_tab$pval_lesions < 0.05 |
		mod_tab$pval_DUKE < 0.05

	mod_idx = include_mod & mod_tab[[t1]] < 2
	# mod_idx = include_mod
	# mod_idx = rep(TRUE, length(include_mod))  # all
	# mod_idx = mod_tab[[t1]] < 2

	# gene_idx = between$meta_genes$tissue[sec_idx] == t1  # relative to
	gene_idx = sec_info$tissue == t1

	mat = cmat[gene_idx, mod_idx]  # temp

	mat = mat[apply(abs(mat), 1, max) > cor_thresh, apply(abs(mat), 2, max) > cor_thresh,
		drop=FALSE
	]
	# mat = as.matrix(mat)

	# Look up if each gene is part of module...
	sec_info$clust[match(rownames(mat), sec_info$id)]

	sec_in_module = sapply(colnames(mat), function(k) {
		sec_info$clust[match(rownames(mat), sec_info$id)] == k
	})
	sec_in_module = as.matrix(sec_in_module)

	sec_in_module_note = matrix("", ncol=ncol(sec_in_module), nrow=nrow(sec_in_module))
	sec_in_module_note[sec_in_module] = "."

	# List module statistics from
	mod_tab[as.integer(colnames(mat)), 1:10]

	# Tissues count matrix for selected modules
	tissue_counts = as.data.frame(mod_tab)[as.integer(colnames(mat)), as.character(tissues)]

	tissue_frac = tissue_counts / mod_tab$mod_size[as.integer(colnames(mat))]

	rlab = matrix("white", ncol=ncol(tissue_counts), nrow=nrow(tissue_counts))
	colnames(rlab) = tissues

	for (i in 1:ncol(rlab)) {
		for (j in 1:nrow(rlab)) {
			rlab[j, i] = addAlpha(tissue_col[i], tissue_frac[j, i])
			# rlab[j, i] = addAlpha(tissue_col[i], tissue_frac[j, i] + 0.1)
		}
	}

	clab = t(matrix(rep(tissue_col[match(t1, tissues)], nrow(mat))))

	# Endocrine name only
	rownames(mat) = sapply(strsplit(rownames(mat), ":"), function(x) x[2])

	try({
		heatmap.3(mat,
			trace="none",
			# mar=c(8, 24),
			main=t1,
			mar=c(8, 32),
			col=colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
			ColSideColors=rlab,
			ColSideColorsSize=2.5,
			RowSideColorsSize=0.7,
			RowSideColors=clab,
			# breaks=seq(-0.5, 0.5, length.out=101),
			breaks=seq(-1.0, 1.0, length.out=101),
			# cexRow=0.075,
			cexRow=0.5,
			xlab="Module", ylab="Endocrine",
			KeyValueName="Pearson's r",
			cellnote=sec_in_module_note, notecol="black"
		)
	})
	dev.off()
}



