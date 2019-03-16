# Locate and write co-expression module GO analysis from WGCNA
rm(list=ls())

data_dir = "~/DataProjects/cross-tissue"  # root of data directory
setwd("~/Google Drive/projects/STARNET/cross-tissue")

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


