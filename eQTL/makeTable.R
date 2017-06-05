rm(list=ls())

library(data.table)
library(biomaRt)

setwd("~/Google Drive/projects/cross-tissue")

sel = "CAD"

eqtl_dir = file.path("eQTL/adjusted", sel)
eqtl_files = list.files(eqtl_dir)

tissues = sapply(strsplit(eqtl_files, "[.]"), function(x) x[4])

eqtls = lapply(file.path(eqtl_dir, eqtl_files), fread)

# Add tissue information to table
for (i in 1:length(eqtls)) {
	eqtls[[i]]$tissue = tissues[i]
}

# Combine all tissue
eqtls = rbindlist(eqtls)

eqtls$ensembl_gene_id = sapply(strsplit(eqtls$gene, "[.]"), function(x) x[1])

# Gene annotations
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
eqtl_annot = getBM(
	filters="ensembl_gene_id",
	attributes=c("ensembl_gene_id", "hgnc_symbol"),
	values=eqtls$ensembl_gene_id,
	mart=ensembl
)

eqtls$hgnc_symbol = eqtl_annot$hgnc_symbol[
	match(eqtls$ensembl_gene_id, eqtl_annot$ensembl_gene_id)
]


# Include only significant eQTLs
eqtls = eqtls[eqtls[["adj.p-value"]] < 0.1, ]

# Sort by p-values
eqtls = eqtls[order(eqtls[["p-value"]]), ]

write.csv(eqtls, file.path("eQTL/tables", paste0(sel, ".csv")),
	row.names=FALSE
)

dim(eqtls)
