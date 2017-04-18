rm(list=ls())

library(RColorBrewer)
library(gplots)
library(qvalue)
library(data.table)
library(psych)

library(devtools)  # for installing heatmap.3
# Load heatmap.3
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")

library(compiler)
enableJIT(3)

data_dir = "/Users/sk/DataProjects/cross-tissue"  # root of data directory
setwd("/Users/sk/Google Drive/projects/cross-tissue")

source("src/models/regr.R")  # regression models
source("src/models/cor.R")
source("src/models/enrichment.R")
source("src/parse.R")


# Parses module data
parseModuleData = function(mod_env)  {
	# Clusters as integers
	mod_env$clust = as.integer(factor(mod_env$bwnet$colors))
	mod_env$meta_genes = parseTranscriptId(mod_env$meta_genes)

	# Rename eigengene matrix
	eigen_gene_names = substring(colnames(mod_env$bwnet$eigengenes), 3)
	eigen_gene_n = match(eigen_gene_names, levels(factor(mod_env$bwnet$colors)))

	colnames(mod_env$bwnet$eigengenes) = eigen_gene_n

	mod_env$bwnet$eigengenes = mod_env$bwnet$eigengenes[, order(eigen_gene_n)]

	return(mod_env)
}


countModuleTissueStat = function(modules) {

	out = list()

	out$tissue_counts = sapply(1:max(modules$clust), function(cl) {
		table(modules$meta_genes$tissue[modules$clust == cl])
	})

	out$purity = apply(out$tissue_counts, 2, max) / apply(out$tissue_counts, 2, sum)

	out$n_tissues = apply(out$tissue_counts, 2, function(col) {
		sum(col > 0)
	})

	out$size = apply(out$tissue_counts, 2, sum)

	return(out)
}



# Load cross-tissue modules
# ------------------------------
between = new.env()
load(file.path(data_dir, "modules/between_within-cross-tissue.RData"),
	between,
	verbose=TRUE)


complete = new.env()
load(file.path(data_dir, "modules/complete-cross-tissue.RData"),
	complete,
	verbose=TRUE)

# Parse module data
between = parseModuleData(between)
complete = parseModuleData(complete)

between_stats = countModuleTissueStat(between)
complete_stats = countModuleTissueStat(complete)


# Load phenotype data
# -----------------------------------

# STARNET phenotype data
pheno = fread(file.path(
	"/Volumes/SANDY/phenotype_data",
	"STARNET_main_phenotype_table.cases.Feb_29_2016.tbl"
))

# Brainshake phenotype data
brainshake = fread(file.path(
	"/Volumes/SANDY/phenotype_data",
	"tz.mat"
))


# Match phenotype to sample order
if (!all(between$patient_ids == complete$patient_ids)) {
	stop("Patient ID mismatch between module detection methods.")
} # else the same patient match can be used for both methods
patient_ids = between$patient_ids

# Match phenotype data tables to
pheno_matched = pheno[match(patient_ids, pheno$starnet.ID), ]
brainshake_matched = brainshake[match(patient_ids, brainshake$id), ]


# Correlation test with selected phenotype data
# -----------------------------------------------------
phenotypes = c("syntax_score", "ndv", "lesions", "DUKE", "BMI", "HbA1c", "LDL", "HDL", "CRP", "p_chol", "bl.glucose")
between_pheno_cor = phenoCorTest(
	mat=between$bwnet$eigengenes,
	# mat=complete$bwnet$eigengenes,
	pheno_matched=pheno_matched,
	phenotypes=phenotypes
)

# Combining p-values into matrix
pheno_cor_pmat = sapply(between_pheno_cor, function(x) x$pval)
rownames(pheno_cor_pmat) = colnames(between$bwnet$eigengenes)
colnames(pheno_cor_pmat) = paste0("pval_", colnames(pheno_cor_pmat))



# # Nominally significant modules
# sig_mods = apply(pmat, 1, min) < 0.05

# pmat = pmat[sig_mods, ]

# hmap = heatmap.2(-log10(pmat), trace="none",
# 	cexRow=1.2,
# 	# cexCol=0.8,
# 	key.title="",
# 	key.xlab=expression("-log"[10] * " p"),
# 	margins=c(10, 8),
# 	col=colorRampPalette(brewer.pal(9, "YlGnBu"))(100)
# )

# between_pheno_cor = lapply(between_pheno_cor, function(tab) {
# 	tab[order(tab$pval), ]
# })

# i = 150
# i = 166
# plot(between$bwnet$eigengenes[, i], pheno_matched$DUKE)

# cor.test(between$bwnet$eigengenes[, i], pheno_matched$syntax_score)


# GO enrichment of modules
# ---------------------------------------------
between_go_enrich = enrichmentGO(between)
# complete_go_enrich = enrichmentGO(complete)

between_go_enrich_filter = filterCommonGOTerms(between_go_enrich, 50)
# complete_go_enrich_filter = filterCommonGOTerms(complete_go_enrich, 50)


# Extract top enrichment terms
top_go_terms_bp = sapply(1:max(between$clust), function(k) {
	enrichment = between_go_enrich_filter$bestPTerms$BP$enrichment
	idx = enrichment$module == k
	top_terms = enrichment$termName[idx][1:5]
	top_terms = top_terms[!is.na(top_terms)]
	return(top_terms)
})

top_go_terms_cc = sapply(1:max(between$clust), function(k) {
	enrichment = between_go_enrich_filter$bestPTerms$CC$enrichment
	idx = enrichment$module == k
	top_terms = enrichment$termName[idx][1:5]
	top_terms = top_terms[!is.na(top_terms)]
	return(top_terms)
})


# between_go_enrich$bestPTerms$BP$enrichment
go_tab = data.frame(
	top_go_bp=sapply(top_go_terms_bp, paste, collapse=";"),
	top_go_cc=sapply(top_go_terms_cc, paste, collapse=";")
)


write.table(between_go_enrich_filter$bestPTerms$BP$enrichment,
	"co-expression/tables/go_bp_tab.csv",
	sep=",",
	quote=FALSE,
	row.names=FALSE)

write.table(between_go_enrich_filter$bestPTerms$CC$enrichment,
	"co-expression/tables/go_cc_tab.csv",
	sep=",",
	quote=FALSE,
	row.names=FALSE)


# head(between_go_enrich$bestPTerms$BP$enrichment)

# i = 98  # Ariella's liver module
# idx = between_go_enrich_filter$bestPTerms$BP$enrichment$module == i
# between_go_enrich_filter$bestPTerms$BP$enrichment[idx, ]

# term_name = "immune response"
# term_name = "lipid metabolic process"
# idx = between_go_enrich_filter$bestPTerms$BP$enrichment$termName == term_name
# between_go_enrich_filter$bestPTerms$BP$enrichment[idx, ]


# # head(between_go_enrich$bestPTerms$BP$enrichment$termName, 100)


# Find CAD-associated genes in cross-tissue modules
# -----------------------------------------------------

# Load Deloukas (2013) genes associated with CAD
delou = read.table(file.path(data_dir, "Deloukas/ng.csv"),
	skip=1,
	sep=",",
	header=TRUE)

# Get CAD gene symbols from proximal loci
cad_genes = delou$Loci_Nearest_Transcript
cad_genes = paste(cad_genes, collapse="/")
cad_genes = strsplit(cad_genes, "/")[[1]]
cad_genes = unique(cad_genes)

# between$clust
# between$meta_genes

cad_gene_bool = between$meta_genes$gene_symbol %in% cad_genes
sum(cad_gene_bool)

# Aggregate CAD gene symbols by module
cad_genes_module = sapply(1:max(between$clust), function(k) {
	idx = between$clust == k & cad_gene_bool

	if (sum(idx) == 0) {
		return(c())
	} else {
		return(between$meta_genes[idx, ])
	}
})

# Count number of found CAD-associated genes
cad_tab = data.frame(n_cad_genes=sapply(cad_genes_module, function(df) max(0, nrow(df))))
cad_tab$cad_genes = sapply(cad_genes_module, function(df) paste(df$gene_symbol, collapse=";"))

# Hypergeometric test of CAD-associated transcripts in each module
m_cad_mRNA = sum(between$meta_genes$gene_symbol %in% cad_genes)
n_non_cad_mRNA = sum(! between$meta_genes$gene_symbol %in% cad_genes)

p_cad_hyper = sapply(1:max(between$clust), function(k) {
	module_size = between_stats$size[k]
	cad_module_mRNA = max(0, nrow(cad_genes_module[[k]]))

	p = 1 - phyper(cad_module_mRNA, m_cad_mRNA, n_non_cad_mRNA, module_size)
	return(p)
})

cad_tab$cad_pval = p_cad_hyper
cad_tab$cad_qvalue = qvalue(cad_tab$cad_pval)$qvalue


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
ciber_freq = lapply(freq_files, function(file_name) {
	file_path = file.path(data_dir, "CIBERSORT/out_freq", file_name)

	tissue = strsplit(file_name, "[.]")[[1]][4]

	# Load CIBERSORT data
	freq = fread(file_path)

	# Parse header separately
	header = read.table(file_path, nrows=1, sep="\t")

	header = unlist(lapply(header, as.character))
	header = c("sample", header)

	# header
	colnames(freq) = header

	# Only numerical fraction entries
	freq_mat = freq[, 2:(ncol(freq) - 3)]
	freq_mat = data.matrix(freq_mat)

	rownames(freq_mat) = freq$sample
	colnames(freq_mat) = paste(tissue, colnames(freq_mat), sep=":")

	return(freq_mat)
})

# Align CIBERSORT frequency data to patient_ids
ciber_freq = lapply(ciber_freq, function(freq_mat) {

	idx = match(
		patient_ids,
		sapply(strsplit(rownames(freq_mat), "_"), function(x) x[2])
		)

	freq_mat = freq_mat[idx, ]

	return(freq_mat)
})

# Test if all rownames agrees with patient IDs (patient_ids), throws error if not.
void = sapply(ciber_freq, function(freq_mat) {
	freq_ids = sapply(strsplit(rownames(freq_mat), "_"), function(x) x[2])
	if(!all(freq_ids == patient_ids, na.rm=TRUE)) {
		stop("Patient ID mismatch")
	}
})

# Combine frequencies across tissue into single matrix
ciber_freq_mat = Reduce(cbind, ciber_freq)
rownames(ciber_freq_mat) = patient_ids

# Remove low variance CIBERSORT features
freq_sd = apply(ciber_freq_mat, 2, sd, na.rm=TRUE)
ciber_freq_mat = ciber_freq_mat[, freq_sd > 1e-9]

# Correlation test
between_ciber_cor = phenoCorTest(
	mat=between$bwnet$eigengenes,
	as.data.frame(ciber_freq_mat)
)

# Add module and feature name to each table
between_ciber_cor = lapply(1:length(between_ciber_cor), function(i) {
	# ith feature
	ciber_name = names(between_ciber_cor)[i]

	df = between_ciber_cor[[i]]

	out_df = cbind(
		data.frame(module=rownames(df), cibersort=ciber_name),
		df)

	return(out_df)
})

# Get p-values for all correlations, in frequency x module matrix
ciber_cor_pmat = sapply(between_ciber_cor, function(x) x$pval)
rownames(ciber_cor_pmat) = colnames(between$bwnet$eigengenes)
# Reomve commas from CIBERSORT features
colnames(ciber_cor_pmat) = gsub(",", "",
	colnames(ciber_freq_mat))

# colnames(pheno_cor_pmat) = paste0("pval_", colnames(pheno_cor_pmat))

# sapply(between_ciber_cor, function(x) x$pval)
# sapp

# FDR
ciber_cor_qval = qvalue(ciber_cor_pmat)$qvalue


# Combine into single table, flat format
between_ciber_cor_all = rbindlist(between_ciber_cor)

# Overwrite q-value from module specific, same calculation as q-value matrix above
between_ciber_cor_all$qval = qvalue(between_ciber_cor_all$p)$qvalue

# Sort by significance
between_ciber_cor_all = between_ciber_cor_all[
	order(between_ciber_cor_all$pval),
]
between_ciber_cor_all$cibersort = gsub(",", "", between_ciber_cor_all$cibersort)  # remove commas from

# Write table
write.table(between_ciber_cor_all,
	"co-expression/tables/cibersort_eigengene_tab.csv",
	sep=",",
	quote=FALSE,
	row.names=FALSE)


# Make table per module
ciber_cor_qval_tissues = sapply(strsplit(colnames(ciber_cor_qval), ":"), function(x) x[1])
top_ciber_cor = sapply(1:nrow(ciber_cor_qval), function(k) {
	# Identify dominant tissue in kth module

	# Filter based on tissue composition of module required for considering CIBERSORT frequencies based on each tissue.
	dominant_tissue = names(which(prop.table(between_stats$tissue_counts[, k]) > 0.05))

	# Map between tissue codes.
	dominant_tissue = replace(dominant_tissue, dominant_tissue == "BLOOD", "BLO")
	dominant_tissue = replace(dominant_tissue, dominant_tissue == "SKLM", "SKM")
	dominant_tissue = replace(dominant_tissue, dominant_tissue == "SF", "SUF")

	# Only the q-values from the dominant tissue
	qvals = ciber_cor_qval[k, ciber_cor_qval_tissues %in% dominant_tissue]

	qvals = qvals[qvals < 0.05]  # FDR of eigengene-frequency correlation
	qvals = sort(qvals)[1:5]  # top fractions
	qvals = qvals[!is.na(qvals)]  # remove missing

	# Format string of q-value results
	if (length(qvals) == 0) {
		return("")
	} else {
		cell_type = names(qvals)

		# combine into single string, separated by ";" and including p-values in parenthesis
		cell_type = paste0(cell_type, "(q=", format(qvals, digits=3), ")")

		cell_type = paste(cell_type, collapse=";")
		return(cell_type)
	}

})

ciber_tab = data.frame(top_cibersort=top_ciber_cor)


# Endocrine factors explaining each module, cross-tissue.
# --------------------------------------------------------

# Load endocrine factors,
endocrine = fread("~/Google Drive/projects/STARNET-endocrine/data/endo_scores_large.csv")

# Get list of vectors of target genes for each endocrine factor
endocrine_targets = strsplit(endocrine$sig_symbols, ";")

# Prepend tissue creating target transcript IDs 
endocrine_targets = lapply(1:length(endocrine_targets), function(k) {
	paste0(endocrine$target_tissue[k], "_", endocrine_targets[[k]])
})

# Endocrine statistics
endocrine_targets_length = sapply(endocrine_targets, length)
n_transcripts = nrow(between$meta_genes)  # Total number of cross-tissue transcripts, used for Hypergeometric test


# Calculate enrichment for each factor in each module using Hypergeometric test
targets_endocrine_enrich = lapply(1:max(between$clust), function(k) {
	message(k)
	idx = between$clust == k

	# Gene symbols for kth module, inculding the tissue codes
	module_symbols = paste0(between$meta_genes$tissue[idx], "_", between$meta_genes$gene_symbol[idx])

	# Calculate overlap between endocrine target genes and kth module
	module_overlap = sapply(endocrine_targets, function(targets) {
		return(length(intersect(targets, module_symbols)))
	})


	# Hypergeometric test, enrichment of module transcripts in endocrine target genes.
	target_p = sapply(1:length(endocrine_targets), function(i) {
		p_out = 1 - phyper(
			module_overlap[i],
			between_stats$size[k],
			n_transcripts - between_stats$size[k],
			endocrine_targets_length[i])
		return(p_out)
	})

	target_tab = data.frame(
		module=k,
		endocrine_factor=endocrine$endocrine_factor,
		from_tissue=endocrine$from_tissue,
		target_tissue=endocrine$target_tissue,
		sec_score=endocrine$sec_score)

	# Is the endocrine factor found in cross-tissue module?
	target_tab$endocrine_in_module = paste0(target_tab$from_tissue, "_", target_tab$endocrine_factor) %in% module_symbols

	# Module overlap enrichment
	target_tab$p = target_p

	target_tab$padj = p.adjust(target_tab$p, method="BH")
	# target_tab$padj = qvalue(p)$qvalue

	# Calculate fraction of module overlap
	target_tab$module_coverage = module_overlap / between_stats$size[k]

	# Endocrine factor inclusion criterion
	# idx = target_tab$padj < 0.05 & target_tab$module_coverage > 0.1
	idx = target_tab$padj < 0.1 & target_tab$module_coverage > 0.05

	# Order significantly enriched endocrine targets by module coverage
	target_sub = target_tab[idx, ]
	target_sub = target_sub[order(target_sub$module_coverage, decreasing=TRUE), ]

	return(target_sub)
})

# targets_endocrine_enrich = Reduce(rbind, targets_endocrine_enrich)
targets_endocrine_enrich = rbindlist(targets_endocrine_enrich)

write.table(targets_endocrine_enrich, "co-expression/tables/endocrine_tab.csv", sep=",",
	# col.names=NA,
	row.names=FALSE,
	quote=FALSE)


# Make string of top-5 highest coverage
endocrine_ids = paste0(targets_endocrine_enrich$from_tissue, "_", targets_endocrine_enrich$endocrine_factor)
top_endocrine = sapply(1:max(between$clust), function(k) {
	idx = targets_endocrine_enrich$module == k
	top_endocrine = endocrine_ids[idx][1:5]  # top-5
	top_endocrine = top_endocrine[!is.na(top_endocrine)]  # remove missing
	return(paste(top_endocrine, collapse=";"))
})

endocrine_tab = data.frame(top_endocrine=top_endocrine)




# Risk enrichment, eQTLs in each cross-tissue module
# ----------------------------------------------------------

# Load eQTL data
eqtl_files = list.files(
	file.path(data_dir, "eQTL/adjusted.final"),
	pattern="fdr-0.10")

exclude_files = c(
	"STARNET.eQTLs.MatrixEQTL.FC.cis.tbl.fdr-0.10",
	"STARNET.eQTLs.MatrixEQTL.MP.cis.tbl.fdr-0.10")

eqtl_files = eqtl_files[!eqtl_files %in% exclude_files]

# Load eQTL data
eqtls = lapply(eqtl_files, function(file_name) {
	d = fread(file.path(data_dir, "eQTL/adjusted.final", file_name))
	return(d)
})

# Name each entry based on tissue
names(eqtls) = sapply(strsplit(eqtl_files, "[.]"), function(x) x[4])  # tissue
names(eqtls) = replace(names(eqtls), names(eqtls) == "Blood", "BLOOD")

# Add tissue name to eqtls tables
for (tissue in names(eqtls)) {
	eqtls[[tissue]]$tissue = tissue
}

# Collapse all eQTL to single table
eqtls = Reduce(rbind, eqtls)


# Find eQTL gene symbols

eqtls$gene_nosuf = sapply(strsplit(eqtls$gene, "[.]"), function(x) x[1])
eqtls$tissue_id = paste0(eqtls$tissue, "_", eqtls$gene_nosuf)

# For improved query performance based on ensembl ID without suffix
between$meta_genes$ensembl_nosuf = sapply(strsplit(between$meta_genes$ensembl, "[.]"), function(x) x[1])
between$meta_genes$tissue_id = paste0(between$meta_genes$tissue, "_", between$meta_genes$ensembl_nosuf)

eqtls$gene_symbol = between$meta_genes$gene_symbol[
	match(
		eqtls$gene_nosuf,
		between$meta_genes$ensembl_nosuf)
]


# Loop over each module 
module_eqtls = lapply(1:max(between$clust), function(k) {
	message(k)
	idx = between$clust == k  # module genes in kth module

	# find eQTLs associated with 
	idx_eqtls = which(eqtls$tissue_id %in% between$meta_genes$tissue_id[idx])

	sub_eqtls = eqtls[idx_eqtls, ]
	sub_eqtls = sub_eqtls[order(sub_eqtls[["p-value"]]), ]

	return(sub_eqtls)
})


# Get eQTL genes per module
module_eqtl_genes = sapply(module_eqtls, function(sub_eqtls) {
	symbols = paste0(sub_eqtls$tissue, "_", sub_eqtls$gene_symbol)
	symbols = unique(na.omit(symbols))
	# module_eqtls = unique(na.omit(sub_eqtls$gene_symbol))
	# return(module_eqtls)
	return(symbols)
})

# Output table
eqtl_tab = data.frame(
	n_eQTL_genes=sapply(module_eqtl_genes, length)
)

eqtl_tab$eQTL_gene_frac = eqtl_tab$n_eQTL_genes / between_stats$size
# eqtl_tab$eQTL_genes = sapply(module_eqtl_genes, paste, collapse=";")
eqtl_tab$top_eQTL_genes = sapply(module_eqtl_genes, function(genes) {
	paste(na.omit(genes[1:5]), collapse=";")
})



# k = 33
# k = 98
# k = 177
# k = 150
# module_eqtl_genes[[k]]

# length(module_eqtls)

# data.frame(eqtls[idx_eqtls, ])


# Combined module table
# -----------------------------------------------

# general statistics
mod_tab = data.frame(
	mod_size=between_stats$size,
	purity=between_stats$purity)

mod_tab = cbind(mod_tab, t(between_stats$tissue_counts))

# CAD-associations
mod_tab = cbind(mod_tab, cad_tab)

# eQTL
mod_tab = cbind(mod_tab, eqtl_tab)

mod_tab = cbind(mod_tab, pheno_cor_pmat)

mod_tab = cbind(mod_tab, endocrine_tab)

# Top GO
mod_tab = cbind(mod_tab, go_tab)

mod_tab = cbind(mod_tab, ciber_tab)

write.table(mod_tab, "co-expression/tables/module_tab.csv", sep=",", col.names=NA)


# mod_tab = mod_tab[order(mod_tab$cad_pval), ]

# plot(mod_tab$eQTL_gene_frac, type="l")

# plot(-log10(mod_tab$cad_pval), mod_tab$eQTL_gene_frac)

plot(-log10(mod_tab$cad_pval), -log10(mod_tab$pval_DUKE))
cor.test(-log10(mod_tab$cad_pval), -log10(mod_tab$pval_DUKE))


plot(-log10(mod_tab$cad_pval), -log10(mod_tab$pval_p_chol))
cor.test(-log10(mod_tab$cad_pval), -log10(mod_tab$pval_p_chol))


head(mod_tab[order(mod_tab$cad_pval), ], 20)
head(mod_tab[order(mod_tab$pval_DUKE), ], 20)
head(mod_tab[order(mod_tab$pval_CRP), ], 20)

plot(-log10(mod_tab$cad_pval), -log10(mod_tab$pval_DUKE))
cor.test(-log10(mod_tab$cad_pval), -log10(mod_tab$pval_DUKE))

plot(-log10(mod_tab$cad_pval), -log10(mod_tab$pval_syntax_score))
cor.test(-log10(mod_tab$cad_pval), -log10(mod_tab$pval_syntax_score))

# head(mod_tab[order(mod_tab$cad_pval), ], 30)


plot(-log10(mod_tab$cad_pval), mod_tab$purity)