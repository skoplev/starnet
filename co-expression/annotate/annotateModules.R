rm(list=ls())

library(RColorBrewer)
library(gplots)
library(qvalue)
library(data.table)

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

# Parses transcript ids from meta_genes table.
# Changes ensembl and gene_symbol column if input data frame and returns new data frame.
# Supports gene symbols with "_".
parseTranscriptId = function(meta_genes) {
	# Get ensembl IDs as last "_" separated string
	meta_genes$ensembl = sapply(
		strsplit(as.character(meta_genes$transcript_id), "_"),
		function(x) {
			x[length(x)]  # last element
	})

	meta_genes$gene_symbol = sapply(
		strsplit(as.character(meta_genes$transcript_id), "_"),
		function(x) {
			if (length(x) == 2) {
				# Gene symbols does not contain "_"
				return(x[1])
			} else {
				# Gene symbol contains "_"
				gene_symbol = paste(x[-length(x)], collapse="_")  # not last element
				# message(gene_symbol)
				return(gene_symbol)
			}
	})
	return(meta_genes)
}

# Parses module data
parseModuleData = function(mod_env)  {
	# Clusters as integers
	mod_env$clust = as.integer(factor(mod_env$bwnet$colors))
	mod_env$meta_genes = parseTranscriptId(mod_env$meta_genes)

	# Rename eigengene matrix
	eigen_gene_names = substring(colnames(mod_env$bwnet$eigengenes), 3)
	eigen_gene_n = match(eigen_gene_names, levels(factor(mod_env$bwnet$colors)))

	colnames(mod_env$bwnet$eigengenes) = eigen_gene_n

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
phenotypes = c("syntax_score", "ndv", "lesions", "DUKE", "BMI", "HbA1c", "LDL", "HDL", "CRP", "p_chol", "bl.glucose")
between_pheno_cor = phenoCorTest(
	mat=between$bwnet$eigengenes,
	# mat=complete$bwnet$eigengenes,
	pheno_matched=pheno_matched,
	phenotypes=phenotypes
)

# pmat = sapply(between_pheno_cor, function(x) x$pval)
# rownames(pmat) = 1:nrow(pmat)

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

i = 150
plot(between$bwnet$eigengenes[, i], pheno_matched$DUKE)
# cor.test(between$bwnet$eigengenes[, i], pheno_matched$syntax_score)



# GO enrichment of modules
# ---------------------------------------------
between_go_enrich = enrichmentGO(between)
complete_go_enrich = enrichmentGO(complete)

between_go_enrich_filter = filterCommonGOTerms(between_go_enrich, 50)
complete_go_enrich_filter = filterCommonGOTerms(complete_go_enrich, 50)

# range(between_go_enrich$bestPTerms$BP$enrichment$enrichmentP)

head(between_go_enrich$bestPTerms$BP$enrichment)
# head(between_go_enrich$bestPTerms$BP$enrichment)
head(between_go_enrich$bestPTerms$BP$enrichment, 30)

i = 98  # Ariella's liver module
# i = 2
# i = 155
# i = 82
# i = 20
# i = 35
# i = 37
# i = 121
# i = 127
# i = 136
# i = 112
# i = 189
i = 78
i = 218
i = 64
idx = between_go_enrich_filter$bestPTerms$BP$enrichment$module == i
between_go_enrich_filter$bestPTerms$BP$enrichment[idx, ]

term_name = "immune response"
term_name = "lipid metabolic process"
idx = between_go_enrich_filter$bestPTerms$BP$enrichment$termName == term_name
between_go_enrich_filter$bestPTerms$BP$enrichment[idx, ]


# head(between_go_enrich$bestPTerms$BP$enrichment$termName, 100)


sort(table(between_go_enrich$bestPTerms$BP$enrichment$termID))
barplot(sort(table(between_go_enrich$bestPTerms$BP$enrichment$termName)))

# head(go_enrich$bestPTerms$BP$enrichment)
# head(go_enrich$bestPTerms$MF$enrichment)
# head(go_enrich$bestPTerms$CC$enrichment)
