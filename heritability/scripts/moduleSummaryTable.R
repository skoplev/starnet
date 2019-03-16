# Some module statistics

rm(list=ls())
library(data.table)

setwd("~/GoogleDrive/projects/STARNET/cross-tissue")


mod_tab = fread("co-expression/tables/module_tab.csv")


modules = fread("co-expression/tables/modules.csv")

eqtl = fread("heritability/eQTL/STARNET_eQTL_FDR05_tissues_modules.csv")

tab = mod_tab[, 1:2]
colnames(tab)[1] = "mod_id"

tab$genes = sapply(tab$mod_id, function(k) {
	message(k)
	n_genes = length(unique(na.omit(modules$ensembl[modules$clust == k])))
	return(n_genes)
})

# Get number of IDs
tab$eqtl_genes = sapply(tab$mod_id, function(k) {
	message(k)
	n_module_eqtl_genes = length(unique(na.omit(eqtl$gene[eqtl$module == k])))
	return(n_module_eqtl_genes)
})

tab$eqtl_perc = (tab$eqtl_genes / tab$genes) * 100

tab$type[mod_tab$purity < 0.95] = "cross-tissue"
tab$type[mod_tab$purity >= 0.95] = "tissue-specific"

write.csv(tab, "heritability/eQTL/module_eqtl.csv",
	row.names=FALSE)


# total number of unique genes
length(unique(na.omit(modules$ensembl)))

# Number of unique module genes that have eQTL
length(unique(na.omit(eqtl$gene[(!is.na(eqtl$module))])))


sum(!is.na(eqtl$module))
