# Locate and write co-expression module GO analysis from WGCNA
rm(list=ls())

data_dir = "~/DataProjects/cross-tissue"  # root of data directory
setwd("~/GoogleDrive/projects/STARNET/cross-tissue")

load(file.path(data_dir, "R_workspaces/annotateModules.RData"), verbose=TRUE)




names(between_go_enrich)

write.table(between_go_enrich$bestPTerms$BP$enrichment,
	"co-expression/tables/mod_GO_BP.tsv",
	sep="\t")

write.table(between_go_enrich$bestPTerms$CC$enrichment,
	"co-expression/tables/mod_GO_CC.tsv",
	sep="\t")

write.table(between_go_enrich$bestPTerms$MF$enrichment,
	"co-expression/tables/mod_GO_MF.tsv",
	sep="\t")

# Write complete p-value matrix
write.table(between_go_enrich$enrichmentP,
	"co-expression/tables/mod_GO_pmat.tsv",
	sep="\t"
)



# Analysis of top-100 terms
# ----------------------------------

# Module selection criteria
mod_tab = fread("co-expression/tables/module_tab.csv")

cross_tissue = mod_tab$purity < 0.95
tissue_specific = mod_tab$purity >= 0.95

sum(cross_tissue)
sum(tissue_specific)


pheno_pval = read.table("pheno/tables/pheno_pval.csv",
	sep=",",
	header=TRUE,
	check.names=FALSE
)

features = c("syntax_score", "DUKE", "case_control_DEG")
pheno_padj = matrix(
	p.adjust(
		data.matrix(pheno_pval[, features]),
		method="BH"),
	ncol=3)
colnames(pheno_padj) = c("SYNTAX", "DUKE", "Case-Ctrl DEG")

fdr = 0.01

cad_modules_idx = apply(pheno_padj < fdr,
	1,
	function(x) sum(x) >= 2)

cad_modules = which(cad_modules_idx)
length(cad_modules)


between_go_enrich$bestPTerms$BP$enrichment
table(between_go_enrich$bestPTerms$BP$enrichment$module)

# sum(between_go_enrich$bestPTerms$BP$enrichment$BonferoniP > 0.05)

go_enrich_top100 = between_go_enrich$bestPTerms$BP$enrichment[between_go_enrich$bestPTerms$BP$enrichment$BonferoniP < 0.01, ]

term = "response to stimulus"
# term = "secretion"

nrow(go_enrich_top100[go_enrich_top100$module %in% which(cross_tissue) & go_enrich_top100$termName == term, ])  # CT
nrow(go_enrich_top100[go_enrich_top100$module %in% cad_modules & go_enrich_top100$termName == term, ])  # CAD modules

nrow(go_enrich_top100[go_enrich_top100$module %in% which(tissue_specific) & go_enrich_top100$termName == term, ])  # TS
