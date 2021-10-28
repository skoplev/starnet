# eQTL table for Nelson SNPs and new STARNET eQTL (estimated by Vamsi).

rm(list=ls())

library(data.table)
library(dplyr)

setwd("~/GoogleDrive/projects/STARNET/cross-tissue")


# Load module data
mod_tab = fread("co-expression/tables/modules.csv")
mod_tab$ensembl_base = sapply(strsplit(mod_tab$ensembl, "[.]"), function(x) x[1])
mod_tab$tissue_ensembl = paste(mod_tab$tissue, mod_tab$ensembl_base, sep="_")


# Load eQTL data from Vamsi run of MatrixEQTL
eqtl_dir = "~/DataProjects/STARNET/vamsi_eQTL/adjusted.final"

# List eQTL table files
eqtl_files = list.files(eqtl_dir, "*.tbl")

tissues = sapply(strsplit(eqtl_files, "_"), function(x) x[1])

# Rename tissue codes
tissues[tissues == "SKM"] = "SKLM"
tissues[tissues == "SUF"] = "SF"
tissues[tissues == "BLO"] = "BLOOD"

# Load tables
eqtls = lapply(file.path(eqtl_dir, eqtl_files), fread)
names(eqtls) = tissues

# Add tissue information to table
for (i in 1:length(eqtls)) {
	eqtls[[i]]$tissue = tissues[i]
}

# Exclude macrophage eQTL
eqtls = eqtls[-which(names(eqtls) == "MAC")]

# Combine tables
eqtls = rbindlist(eqtls)

eqtls$tissue_ensembl = paste(eqtls$tissue, eqtls$gene, sep="_")  # tissue ensembl IDs for matching with module assignments

# Replace V1 colname duplicate, generating valid data frame that can be used with merge()
colnames(eqtls)[1] = "rown"


# Add module assigments
eqtls = merge(
	eqtls,
	mod_tab[, c("tissue_ensembl", "clust", "gene_symbol")],
	all.x=TRUE,
	by="tissue_ensembl"
)

# Global sort by p-value
eqtls = eqtls[order(eqtls[["p-value"]]), ]



# Load nelson SNPs 
nelson_snp = read.table("~/DataProjects/cross-tissue/GWAS/Nelson/variants_q05.txt")[, 1]
nelson_snp = as.character(nelson_snp)

nelson_snp = gsub("\\*", "", nelson_snp)  # remove wildcard indications extracted from supplementary table

length(grep("rs", nelson_snp))  # number of regular SNP IDs


eqtls_nelson = eqtls[eqtls$snpid %in% nelson_snp, ]

eqtls_nelson_format = eqtls_nelson[, c("snpid", "gene", "gene_symbol", "tissue", "clust", "beta", "t-stat", "p-value", "padj_fdr")]


write.csv(eqtls_nelson_format, "eQTL/tables/nelson_eQTL.csv",
	row.names=FALSE)


barplot(table(eqtls_nelson_format$tissue))
length(unique(eqtls_nelson_format$snpid))
