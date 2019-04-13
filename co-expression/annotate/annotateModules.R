# rm(list=ls())

options(java.parameters = "-Xmx8g")  # Max Java memory heap size, for rcausal

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
setwd("	~/GoogleDrive/projects/STARNET/cross-tissue")

source("src/models/regr.R")  # regression models
source("src/models/cor.R")
source("src/models/enrichment.R")
source("src/parse.R")
# source("src/models/bayesNet.R")



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


# Load normalized gene expression data, used for evaluating Bayesian networks within cross-tissue modules
# --------------------------------------------------------------
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


olink = fread("~/Google Drive/projects/STARNET/olink/blood/STARNET_protein_uppsala/untitled.csv")
olink[olink == "NAN"] = NA


# Match phenotype to sample order
if (!all(between$patient_ids == complete$patient_ids)) {
	stop("Patient ID mismatch between module detection methods.")
} # else the same patient match can be used for both methods
patient_ids = between$patient_ids

# Match phenotype data tables to
pheno_matched = pheno[match(patient_ids, pheno$starnet.ID), ]
brainshake_matched = brainshake[match(patient_ids, brainshake$id), ]
olink_matched = olink[match(patient_ids, olink$NPX)]


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


brainshake_cor = phenoCorTest(
	mat=between$bwnet$eigengenes,
	pheno_matched=brainshake_matched
)

olink_cor = phenoCorTest(
	mat=between$bwnet$eigengenes,
	pheno_matched=olink_matched)


olink_pmat = sapply(olink_cor[-1], function(x) x$pval)
olink_qmat = sapply(olink_cor[-1], function(x) x$qval)

colnames(olink_pmat) = sapply(strsplit(colnames(olink_pmat), "_"), function(x) x[2])





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

cad_genes = getCADGenes(data_dir)

# Test for enrichment
cad_tab = hyperGeometricModuleTest(between, cad_genes)
sum(cad_tab$qvalue < 0.1)

colnames(cad_tab) = paste0("CAD_", colnames(cad_tab))


# GWAS traits in general
# ---------------------------------------------------------------------

gwas = fread("~/DataBases/GWAS/gwas_catalog_v1.0-associations_e88_r2017-05-29.tsv")


# Overview of GWAS traits, exploration of traits
gwas_counts = sort(table(gwas[["DISEASE/TRAIT"]]))
gwas_traits = names(gwas_counts)

write.csv(rev(gwas_counts), "GWAS_counts.csv", row.names=FALSE)

gwas_counts[grep("diabetes", gwas_traits, ignore.case=TRUE)]
gwas_counts[grep("lipid", gwas_traits, ignore.case=TRUE)]

gwas_counts[grep("BMI", gwas_traits, ignore.case=TRUE)]
gwas_counts[grep("body", gwas_traits, ignore.case=TRUE)]
gwas_counts[grep("waist", gwas_traits, ignore.case=TRUE)]

gwas_counts[grep("LDL", gwas_traits, ignore.case=TRUE)]
gwas_counts[grep("HDL", gwas_traits, ignore.case=TRUE)]
gwas_counts[grep("cholesterol", gwas_traits, ignore.case=TRUE)]
gwas_counts[grep("triglyceride", gwas_traits, ignore.case=TRUE)]
gwas_counts[grep("glucose", gwas_traits, ignore.case=TRUE)]
gwas_counts[grep("insulin", gwas_traits, ignore.case=TRUE)]

# gwas[gwas[["FIRST AUTHOR"]] == "Shungin D", ]
# table(gwas[gwas[["FIRST AUTHOR"]] == "Shungin D", ][["DISEASE/TRAIT"]])

trait = "Visceral fat"
trait = "Phospholipid levels"
trait = "Glucose homeostasis traits"
idx = gwas[["DISEASE/TRAIT"]] == trait
data.frame(gwas[idx, ])


# GWAS traits/phenotypes to include in analysis
traits = c(
	"Type 2 diabetes",
	"Phospholipid levels (plasma)",
	"Lipid metabolism phenotypes",
	"LDL cholesterol",
	"HDL cholesterol",
	"Cholesterol, total",
	"Waist-to-hip ratio adjusted for body mass index",
	"Glucose homeostasis traits",
	"Triglycerides",
	"Fasting glucose-related traits",
	"Lipid traits",
	"Body mass index",
	"Fasting plasma glucose",
	"Blood pressure",
	"Systolic blood pressure",
	"Diastolic blood pressure",
	"Coronary artery calcification",
	"Hypertension",
	"Myocardial infarction",
	"Visceral fat",
	"Response to statin therapy",
	"C-reactive protein"
)

gwas_genes = list()
gwas_genes = lapply(traits, function(trait) getGWAS(gwas, trait))
names(gwas_genes) = traits

gwas_tabs = lapply(gwas_genes, function(genes) hyperGeometricModuleTest(between, genes))

gwas_enrich_pval = sapply(gwas_tabs, function(x) x$pval)
gwas_enrich_pval = as.data.frame(gwas_enrich_pval)
gwas_enrich_pval$CAD = cad_tab$CAD_pval

rownames(gwas_enrich_pval) = 1:nrow(gwas_enrich_pval)

# Combine into single data frame
gwas_tabs_comb = lapply(1:length(gwas_tabs), function(i) {
	tab = gwas_tabs[[i]]
	colnames(tab) = paste0(names(gwas_tabs[i]), "_", colnames(tab))
	return(tab)
})
gwas_tabs_comb = Reduce(cbind, gwas_tabs_comb)

# Secreted proteins
# --------------------------------------------------

secreted_proteins = fread(file.path(data_dir, "Uniprot/uniprot_human_secreted_proteins.tab"))

secreted_proteins_symbols = strsplit(
	paste(
		secreted_proteins[["Gene names  (primary )"]],
		collapse=";"),
	";")[[1]]

secreted_proteins_symbols = trimws(secreted_proteins_symbols)
secreted_proteins_symbols = unique(secreted_proteins_symbols)
secreted_proteins_symbols = secreted_proteins_symbols[secreted_proteins_symbols != ""]


sec_tab = hyperGeometricModuleTest(between, secreted_proteins_symbols)

colnames(sec_tab) = paste0("secreted_protein_", colnames(sec_tab))



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
# Remove commas from CIBERSORT features
colnames(ciber_cor_pmat) = gsub(",", "",
	colnames(ciber_freq_mat))

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

targets_endocrine_enrich = rbindlist(targets_endocrine_enrich)

write.table(targets_endocrine_enrich,
	"co-expression/tables/endocrine_tab.csv",
	sep=",",
	row.names=FALSE,
	quote=FALSE)

write.table(targets_endocrine_enrich[targets_endocrine_enrich$endocrine_in_module, ],
	"co-expression/tables/endocrine_tab_in_module.csv",
	sep=",",
	row.names=FALSE)

# Make string of top-5 highest coverage
endocrine_ids = paste0(targets_endocrine_enrich$from_tissue, "_", targets_endocrine_enrich$endocrine_factor)
top_endocrine = sapply(1:max(between$clust), function(k) {
	idx = targets_endocrine_enrich$module == k
	top_endocrine = endocrine_ids[idx][1:5]  # top-5
	top_endocrine = top_endocrine[!is.na(top_endocrine)]  # remove missing
	return(paste(top_endocrine, collapse=";"))
})

# Count total number of endocrine factors found in each module
within_module_endocrines = table(
	factor(
		targets_endocrine_enrich[targets_endocrine_enrich$endocrine_in_module, ]$module,
		levels=1:max(between$clust)
	)
)

endocrine_tab = data.frame(
	top_endocrine=top_endocrine,
	n_within_mod_endocrines=as.vector(within_module_endocrines))


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

	sub_eqtls$module = k

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
eqtl_tab$top_eQTL_genes = sapply(module_eqtl_genes, function(genes) {
	paste(na.omit(genes[1:5]), collapse=";")
})

# write.table(rbindlist(module_eqtls), "co-expression/tables/eQTL_tab.csv",
# 	col.names=NA,
# 	sep=",")

# CAD genes overlap with eQTLs
eqtl_tab$eQTL_cad_genes = sapply(module_eqtl_genes, function(tissue_gene_ids) {
	gene_symbols = sapply(strsplit(as.character(tissue_gene_ids), "_"), function(x) x[2])

	eqtl_cad_genes = tissue_gene_ids[gene_symbols %in% cad_genes]
	return(paste0(eqtl_cad_genes, collapse=";"))
})

eqtl_tab$n_eQTL_cad_genes = sapply(eqtl_tab$eQTL_cad_genes, function(genes) length(strsplit(genes, ";")[[1]]))


# Key driver analysis, see bayesNet1.R
# -------------------------------------------


# load results table
kda_results = read.table(
	file.path(bayes_dir, "kda", paste0(kda_label, ".results.txt")),
	header=TRUE
)

# k = 98
# k = 20
# k = 150
# k = 33
# k = 177
# k = 74
# k = 65
# k = 117
# k = 162
# k = 82
# k = 150
# kda_results[kda_results$MODULE == k, ]

top_key_drivers = sapply(1:max(between$clust), function(k) {
	idx = kda_results$MODULE == k
	top_nodes = na.omit(kda_results$NODE[idx][1:5])
	top_nodes = as.character(top_nodes)
	top_nodes = sapply(strsplit(top_nodes, "_"), function(x) paste(x[1:2], collapse="_"))
	return(paste(top_nodes, collapse=";"))
})


# Key drivers that are also eQTLs
key_driver_eQTL_genes = sapply(1:max(between$clust), function(k) {
	idx = kda_results$MODULE == k
	nodes = as.character(kda_results$NODE[idx])
	nodes = sapply(strsplit(nodes, "_"), function(x) paste(x[1:2], collapse="_"))

	eqtl_nodes = nodes[nodes %in% module_eqtl_genes[[k]]]

	return(paste(eqtl_nodes, collapse=";"))
})

kd_tab = data.frame(top_key_drivers=top_key_drivers,
	key_driver_eQTL_genes=key_driver_eQTL_genes,
	n_key_driver_eQTL_genes=sapply(strsplit(key_driver_eQTL_genes, ";"), length)
	)


# Combine module tables
# -----------------------------------------------

# general statistics
mod_tab = data.frame(
	mod_size=between_stats$size,
	purity=between_stats$purity)

mod_tab = cbind(mod_tab, t(between_stats$tissue_counts))

# CAD-associations, enrichment
mod_tab = cbind(mod_tab, cad_tab)

# Secretion enrichment
mod_tab = cbind(mod_tab, sec_tab)

# eQTL
mod_tab = cbind(mod_tab, eqtl_tab)

# Phenotype correlations
mod_tab = cbind(mod_tab, pheno_cor_pmat)

# Endocrine factors
mod_tab = cbind(mod_tab, endocrine_tab)


# Top GO
mod_tab = cbind(mod_tab, go_tab)

# CIBERSORT correlations
mod_tab = cbind(mod_tab, ciber_tab)

mod_tab = cbind(mod_tab, kd_tab)

# Other GWAS enrichment results
mod_tab = cbind(mod_tab, gwas_tabs_comb)

write.table(mod_tab, "co-expression/tables/module_tab.csv", sep=",", col.names=NA)



# Collect pvalue for features
# --------------------------------
library(GO.db)


pval_mats = list()  # matrix of p-values for each module 
pval_mats$pheno = pheno_cor_pmat

pval_mats$GWAS = gwas_enrich_pval

include_go_terms = c(
	"GO:0006629", # lipid metabolic process
	"GO:0097006",  # regulation of plasma lipoprotein particle levels
	"GO:0006695",  # cholesterol biosynthetic process
	"GO:0006954",  # inflammatory response
	"GO:0045321",  # leukocyte activation
	"GO:0045087",  # innate immune response
	"GO:0002250",  # adaptive immune response
	"GO:0034377",  # plasma lipoprotein particle assembly
	"GO:0098609",  # cell-cell adhesion
	"GO:0001816",  # cytokine production
	"GO:0019752",  # carboxylic acid metabolic process
	"GO:0042981",  # regulation of apoptotic process
	"GO:0030334"  # regulation of cell migration
)

pval_mats$GO = between_go_enrich$enrichmentP[, colnames(between_go_enrich$enrichmentP) %in% include_go_terms]
colnames(pval_mats$GO) = Term(include_go_terms)


include_ciber = c(
	"AOR:blood:macrophage",
	"AOR:blood:monocyte",
	"MAM:blood:macrophage",
	"MAM:blood:monocyte"
)

pval_mats$CIBERSORT = ciber_cor_pmat[, colnames(ciber_cor_pmat) %in% include_ciber]


# Clamp p-values
pval_mats = lapply(pval_mats, function(pval) {
	min_pval = 10^-10
	pval[pval < min_pval] = min_pval
	return(pval)
})





# Module annotation significance plots
# ---------------------------------

# Nominally significant modules

pval_all = Reduce(cbind, pval_mats)
# sig_mods = apply(pval_all, 1, min) < 0.001

# sig_mods = apply(pval_mats$pheno, 1, min) < 0.01 &
# 	apply(pval_mats$GWAS, 1, min) < 0.01

sig_mods = apply(pval_mats$pheno, 1, min) < 0.01 &
	apply(pval_mats$GWAS, 1, min) < 0.01 &
	apply(pval_mats$GO, 1, min) < 0.01

# sig_mods = apply(pval_mats$pheno, 1, min) < 0.05 &
# 	apply(pval_mats$GWAS, 1, min) < 0.05 &
# 	apply(pval_mats$GO, 1, min) < 0.05


# sig_mods = apply(pval_mats$pheno, 1, min) < 0.01 |
# 	apply(pval_mats$GWAS, 1, min) < 0.01 |
# 	apply(pval_mats$GO, 1, min) < 0.01


sum(sig_mods)

heatmap.2(
	# -log10(t(pval_all)),
	-log10(t(pval_all[sig_mods, ])),
	trace="none",
	col=colorRampPalette(brewer.pal(9, "YlGnBu"))(100),
	breaks=seq(0, 6, length.out=101)
)

# mod_order = hclust(dist(-log(pval_all)))$order
hc = hclust(dist(-log(pval_all[sig_mods, ])))
mod_order = hc$order

pdf("co-expression/annotate/plots/module_dendrogram.pdf", width=12, height=5)
plot(hc)
dev.off()

# sig_mods = apply(pval_mats$pheno, 1, min) < 0.01


pdf("co-expression/annotate/plots/pheno_cor.pdf", height=5)
# mat = pval_mats$pheno[mod_order, ]
mat = pval_mats$pheno[sig_mods, ][mod_order, ]
hmap = heatmap.2(
	-log10(t(mat)),
	Colv=FALSE,
	trace="none",
	key.title="",
	key.xlab=expression("-log"[10] * " p"),
	col=colorRampPalette(brewer.pal(9, "YlOrRd"))(100),
	breaks=seq(0, 5, length.out=101),
	margins=c(12, 8),
	cexCol=0.3,
	cexRow=0.5,
	ylab="Phenotype"
)
dev.off()


# mat = pval_mats$GWAS
# sig_mods = apply(gwas_enrich_pval, 1, min) < 0.01
# mat = mat[idx, ]

# sum(idx)

pdf("co-expression/annotate/plots/GWAS_enrichment.pdf", height=5)
# mat = pval_mats$GWAS[mod_order, ]
mat = pval_mats$GWAS[sig_mods, ][mod_order, ]
heatmap.2(
	-log10(t(mat)),
	Colv=FALSE,
	trace="none",
	col=colorRampPalette(brewer.pal(9, "Blues"))(100),
	breaks=seq(0, 5, length.out=101),
	margins=c(12, 8),
	cexCol=0.3,
	cexRow=0.5,
	key.xlab="-log10 p",
	xlab="Module", ylab="GWAS"
)
dev.off()


pdf("co-expression/annotate/plots/GO_enrichment.pdf", height=5)
# mat = pval_mats$GWAS[mod_order, ]
mat = pval_mats$GO[sig_mods, ][mod_order, ]
heatmap.2(
	-log10(t(mat)),
	Colv=FALSE,
	trace="none",
	col=colorRampPalette(brewer.pal(9, "Greens"))(100),
	breaks=seq(0, 5, length.out=101),
	margins=c(12, 8),
	cexCol=0.3,
	cexRow=0.5,
	key.xlab="-log10 p",
	xlab="Module",
	ylab="GO"
)
dev.off()

pdf("co-expression/annotate/plots/CIBERSORT_enrichment.pdf", height=5)
# mat = pval_mats$GWAS[mod_order, ]
mat = pval_mats$CIBERSORT[sig_mods, ][mod_order, ]
heatmap.2(
	-log10(t(mat)),
	Colv=FALSE,
	trace="none",
	col=colorRampPalette(brewer.pal(9, "Purples"))(100),
	breaks=seq(0, 5, length.out=101),
	margins=c(12, 8),
	cexCol=0.3,
	cexRow=0.5,
	key.xlab="-log10 p",
	xlab="Module",
	ylab="CIBERSORT"
)
dev.off()


tissue_col = brewer.pal(9, "Set1")[-6]

pdf("co-expression/annotate/plots/tissue_distribution.pdf", height=2.5)
mat = between_stats$tissue_counts
colnames(mat) = 1:ncol(mat)
mat = mat[, sig_mods][, mod_order]
barplot(mat,
	ylim=c(0, 2000),
	border=NA,
	names.arg=colnames(mat),
	cex.names=0.3,
	las=2,
	col=tissue_col)
legend("topright", legend=rownames(between_stats$tissue_counts),
	pch=15,
	col=tissue_col,
	cex=0.4
)
dev.off()




# Secreted protein enrichment


cross_tissue = mod_tab$purity < 0.95  # Cross tissue definition
# cross_tissue = mod_tab$purity < 0.99  # Cross tissue definition

d = list(
	"Cross-tissue"=-log10(sec_tab$secreted_protein_qvalue[cross_tissue]),
	"Tissue-specific"=-log10(sec_tab$secreted_protein_qvalue[!cross_tissue])
)

d = lapply(d, function(x) {
	x[is.infinite(x)] = 16
	return(x)
})

boxplot(d)

# wilcox.test(d[[1]], d[[2]])



pdf("co-expression/plots/module_purity_tests.pdf", width=12)
# Selectio criteria for testing cross-tissue fraction
sel = list(
	"Secreted"=mod_tab$secreted_protein_qvalue < 0.1,
	"CAD"=mod_tab$CAD_qvalue < 0.1,

	"SYNTAX"=mod_tab$pval_syntax_score < 0.05,
	"DUKE"=mod_tab$pval_DUKE < 0.05,
	"ndv"=mod_tab$pval_ndv < 0.05,
	"lesions"=mod_tab$pval_lesions< 0.05,

	"BMI"=mod_tab$pval_BMI < 0.05,
	"HbA1c"=mod_tab$pval_HbA1c < 0.05,
	"LDL"=mod_tab$pval_LDL < 0.05,
	"HDL"=mod_tab$pval_HDL < 0.05,
	"CRP"=mod_tab$pval_CRP < 0.05,
	"p_chol"=mod_tab$pval_p_chol < 0.05
)

par(mfrow=c(2, floor(length(sel)/2)))

for (i in 1:length(sel)) {

	d = list(
		"-"=(1-mod_tab$purity[!sel[[i]]]) * 100,
		"+"=(1-mod_tab$purity[sel[[i]]]) * 100
	)

	boxplot(d,
		outline=FALSE, ylab="Cross-tissue %",
		main=paste0(
			names(sel)[i],
			" p=",
			format(wilcox.test(d[[1]], d[[2]])$p.value, digits=4)
		),
		border="grey"
	)
	points(
		jitter(rep(1, length(d[[1]])), 4),
		d[[1]])
	points(
		jitter(rep(2, length(d[[2]])), 4),
		d[[2]])

}
dev.off()



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