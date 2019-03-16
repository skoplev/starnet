# Script that collects a table with all eSNPs in modules
# Annotate the tissue of eQTL association and whether modules are cross-tissue or tissue-specific
# SNP and chromosome location

rm(list=ls())

library(data.table)
library(biomaRt)
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
ensembl_map = getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"), mart=ensembl)


setwd("~/Google Drive/projects/STARNET/cross-tissue")

# Load table with modules for each transcript
mod = fread("co-expression/tables/modules.csv")

# Decapitalize blood tissue encoding
mod$tissue[mod$tissue == "BLOOD"] = "Blood"

mod$tissue_ensembl = paste(mod$tissue,
	sapply(strsplit(mod$ensembl, "[.]"), function(x) x[1]),
	sep="_")  # for identifying tissue-transcript pairs


# Load module annotation table
mod_annot = fread("co-expression/tables/module_tab.csv")



# Load genpotype/SNP metadata
snp_info = fread("~/DataProjects/STARNET/genotype/STARNET.v3.genotype.dose.info")

colnames(snp_info)[colnames(snp_info) == "marker_id"] = "SNP"

# Load eQTL table
tissues = c("AOR", "MAM", "LIV", "SKLM", "VAF", "SF", "Blood")
eqtl = lapply(tissues, function(t) {
	d = fread(
		paste0(
			"~/DataProjects/STARNET/eQTL/",
			"STARNET.eQTLs.MatrixEQTL.",
			t,
			".cis.tbl"
		)
	)
	d$tissue = t
	return(d)
})

eqtl = rbindlist(eqtl)

fdr = 0.05

# Filter by genome-wide significance, per tissue
eqtl = eqtl[eqtl$adj.p < fdr, ]


# Gene annotation
ensembl_id = sapply(strsplit(eqtl$gene, "[.]"), function(x) x[1])  # ensembl without version
eqtl$tissue_ensembl = paste(eqtl$tissue,
	ensembl_id,
	sep="_")

# Map to HGNC gene symbols
eqtl$hgnc_symbol = ensembl_map$hgnc_symbol[match(ensembl_id, ensembl_map$ensembl_gene_id)]


# Annotate eSNP table with modules.
# Find module number for each transcript.
idx = match(eqtl$tissue_ensembl, mod$tissue_ensembl)
eqtl$module = mod$clust[idx]


# Add cross-tissue, tissue-specific information for each tissue
module_group = rep(NA, nrow(mod_annot))
module_group[mod_annot$purity < 0.95] = "cross-tissue"
module_group[mod_annot$purity >= 0.95] = "tissue-specific"
eqtl$module.group = module_group[eqtl$module]


# Add genotype/SNP information
eqtl_SNP = merge(eqtl, snp_info, by="SNP", all.x=TRUE)

# Order by significance
eqtl_SNP = eqtl_SNP[order(eqtl_SNP[["p-value"]]), ]


write.csv(eqtl_SNP, "heritability/eQTL/STARNET_eQTL_FDR05_tissues_modules.csv",
	row.names=FALSE,
	quote=FALSE
)


head(eqtl[eqtl$module == 98, ], 100)

head(eqtl_SNP[eqtl_SNP$module == 98, ], 100)
head(eqtl_SNP[eqtl_SNP$module == 78, ], 100)
head(eqtl_SNP[eqtl_SNP$module == 199, ], 100)

