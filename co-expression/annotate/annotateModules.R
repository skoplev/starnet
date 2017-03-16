rm(list=ls())

library(RColorBrewer)
library(gplots)

library(devtools)  # for installing heatmap.3
# Load heatmap.3
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")

library(compiler)
enableJIT(3)

data_dir = "/Users/sk/DataProjects/cross-tissue"  # root of data directory
setwd("/Users/sk/Google Drive/projects/cross-tissue")

source("src/models/regr.R")  # regression models
source("src/models/cor.R")

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

enrichmentGO = function(modules) {
	# Map ensembl IDs
	ensembl = sapply(strsplit(modules$meta_genes$ensembl, "[.]"), function(x) x[1])

	# Load map of Ensembl -> ENTREX IDs
	entrez_map = select(org.Hs.eg.db, ensembl, "ENTREZID", "ENSEMBL")

	entrez = entrez_map$ENTREZID[match(ensembl, entrez_map$ENSEMBL)]

	go_enrich = GOenrichmentAnalysis(between$clust,
		entrez,
		# between$meta_genes$ensembl_nosuf,
		organism="human",
		removeDuplicates=FALSE,
		nBestP=50,
		# pCut=0.05  # doesn't work
	)

	return(go_enrich)
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


# between_fits_pheno = fitLinearEigenPheno(pheno_matched, between$bwnet$eigengenes)
# between_fits_pheno[["syntax_score"]]
# between_fits_pheno[["LDL"]]
# summary(between_fits_pheno[["syntax_score"]])




phenotypes = c("syntax_score", "ndv", "lesions", "DUKE", "BMI", "HbA1c", "LDL", "HDL", "CRP", "p_chol", "bl.glucose")
between_pheno_cor = phenoCorTest(
	mat=between$bwnet$eigengenes,
	# mat=complete$bwnet$eigengenes,
	pheno_matched=pheno_matched,
	phenotypes=phenotypes
)

# pheno_cor[[phenotype]] = tab[order(tab$pval), ]

# qmat = sapply(between_pheno_cor, function(x) x$qval)


# points(rep(1, sum(sig_mods)), between_stats$purity[sig_mods], pch=16)
# points(rep(2, sum(!sig_mods)), between_stats$purity[!sig_mods], pch=16)
# 	# between_stats$purity[!sig_mods])

# t.test(between_stats$purity[sig_mods], between_stats$purity[!sig_mods])

# ks.test(between_stats$purity[sig_mods], between_stats$purity[!sig_mods])
# wilcox.test(between_stats$purity[sig_mods], between_stats$purity[!sig_mods])
# wilcox.test(between_stats$purity[sig_mods], between_stats$purity[!sig_mods])
# # kruskal.test(between_stats$purity[sig_mods], between_stats$purity[!sig_mods])




# boxplot(-log10(between_stats$purity[sig_mods]), -log10(between_stats$purity[!sig_mods]))

# wilcox.test(between_stats$purity[sig_mods], between_stats$purity[!sig_mods])
# t.test(-log10(between_stats$purity[sig_mods]), -log10(between_stats$purity[!sig_mods]))
# # ks.test(between_stats$purity[sig_mods], between_stats$purity[!sig_mods])

# boxplot(between_stats$purity[sig_mods], between_stats$purity[!sig_mods])

# boxplot(between_stats$n_tissues[sig_mods], between_stats$n_tissues[!sig_mods])

# boxplot(between_stats$purity[sig_mods & between_stats$size < 2000], between_stats$purity[!sig_mods & between_stats$size < 2000])
# t.test(between_stats$purity[sig_mods & between_stats$size < 2000], between_stats$purity[!sig_mods & between_stats$size < 2000])

# barplot(
# 	c(mean(between_stats$purity[sig_mods]),
# 		mean(between_stats$purity[!sig_mods]))
# )



# heatmap.2(log10(t(tissue_counts) + 1),
# 	col=colorRampPalette(brewer.pal(9, "Blues"))(100),
# 	trace="none")

# col=rep("black", length(sig_mods))
# col[sig_mods] = "red"
# plot(log10(between_stats$size), between_stats$purity, col=col, pch=16)


between_pheno_cor = lapply(between_pheno_cor, function(tab) {
	tab[order(tab$pval), ]
})


between_pheno_cor[["syntax_score"]]
between_pheno_cor[["lesions"]]
between_pheno_cor[["DUKE"]]
between_pheno_cor[["ndv"]]

between_pheno_cor[["LDL"]]
between_pheno_cor[["Age"]]
between_pheno_cor[["bl.glucose"]]
between_pheno_cor[["HbA1c"]]
between_pheno_cor[["CRP"]]
between_pheno_cor[["Age"]]
between_pheno_cor[["p_chol"]]
between_pheno_cor[["SBP"]]
# between_pheno_cor[["TG"]]


# i = 45
# i = 98
# i = 120
# i = 115
# i = 186
# i = 104

# i = 20

i = 150
plot(between$bwnet$eigengenes[, i], pheno_matched$DUKE)
# cor.test(between$bwnet$eigengenes[, i], pheno_matched$syntax_score)

# i = 78
i = 220
plot(between$bwnet$eigengenes[, i], pheno_matched$lesions)


# i = 194
# plot(between$bwnet$eigengenes[, i], pheno_matched$ndv)



i = 98
plot(between$bwnet$eigengenes[, i], pheno_matched$syntax_score)

cor.test(between$bwnet$eigengenes[, i], pheno_matched$syntax_score)



# Remove GO terms that are prevalent across modules
removePrevalentTerms = function(enrich_tab, max_preval) {
	mod_counts = table(enrich_tab$termID)
	include_idx = mod_counts < max_preval
	message("Excluding GO terms: ", sum(!include_idx))
	include_terms = names(which(include_idx))
	return(enrich_tab[enrich_tab$termID %in% include_terms, ])
}

# Filter out GO terms from WGCNA GO enrichment analysis
filterCommonGOTerms = function(go, max_preval=50) {
	go$bestPTerms$BP$enrichment = removePrevalentTerms(go$bestPTerms$BP$enrichment, max_preval)
	go$bestPTerms$CC$enrichment = removePrevalentTerms(go$bestPTerms$CC$enrichment, max_preval)
	go$bestPTerms$MF$enrichment = removePrevalentTerms(go$bestPTerms$MF$enrichment, max_preval)
	return(go)
}


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
