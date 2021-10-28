# Some module statistics

rm(list=ls())
library(data.table)

setwd("~/GoogleDrive/projects/STARNET/cross-tissue")


# Load STARNET cis-eQTL data, estimated by Vamsi
# -----------------------------------------------------
getEqtlNew = function() {

	cis_eqtl_dir = "~/DataProjects/STARNET/vamsi_eQTL/adjusted.final"
	cis_eqtl_files = list.files(cis_eqtl_dir, pattern="*.tbl")

	tissues = sapply(strsplit(cis_eqtl_files, "_"), function(x) x[1])

	# Rename tissue codes
	tissues[tissues == "SKM"] = "SKLM"
	tissues[tissues == "SUF"] = "SF"
	tissues[tissues == "BLO"] = "BLOOD"


	cis_eqtl = lapply(
		cis_eqtl_files,
		function(file_name) {
			d = fread(file.path(cis_eqtl_dir, file_name))
			# d = d[d$padj_fdr < 0.05, ]  # FDR < 5%
			# d = d[d$padj_fdr < 0.01, ]  # FDR < 1%
			# d = d[d$padj_fdr < 0.001, ]  # FDR < .1%
			d = d[d$padj_fdr < 0.0001, ]  # FDR < .01%
			d = d[order(d[["p-value"]]), ]
			return(d)
	})
	names(cis_eqtl) = tissues

	# Add tissue information to table
	for (i in 1:length(cis_eqtl)) {
		cis_eqtl[[i]]$tissue = tissues[i]
	}

	# Exclude macrophage eQTL
	cis_eqtl = cis_eqtl[-which(names(cis_eqtl) == "MAC")]

	# Combine tables
	cis_eqtl = rbindlist(cis_eqtl)

	cis_eqtl$tissue_ensembl_id = paste(cis_eqtl$tissue, cis_eqtl$gene, sep="_")  # tissue ensembl IDs for matching with module assignments

	return(unique(cis_eqtl$tissue_ensembl_id))
}

cis_eqtl_all = getEqtlNew()  # Vamsi's eQTL


# Load TF definition from Lambert et al
# -------------------------------------------------------
tf_symbols = as.character(read.table("transcription-factors/lambert/TF_names_v_1.01.txt")$V1)



# Load key driver analysis results
# --------------------------------
# kda = fread("co-expression/annotate/grn_vamsi_eqtl/kda/modules.results.txt")
kda = fread("co-expression/annotate/grn_vamsi_eqtl/kda/modules.directed.results.txt")


kda = kda[kda$FDR < 0.05, ]
# kda = kda[kda$FDR < 0.0001, ]


# Load module table
mod_tab = fread("co-expression/tables/module_tab.csv")

# Load meta gene table
modules = fread("co-expression/tables/modules.csv")

modules$tissue_ensembl_id = paste0(modules$tissue, "_",
	sapply(strsplit(modules$ensembl, "[.]"), function(x) x[1])
)


# Regulator status in GENIE3 analysis. From geneRegulatoryNetworkInference.R script
modules$regulator = FALSE
modules$regulator[modules$gene_symbol %in% tf_symbols] = TRUE
sum(modules$regulator)

modules$regulator[modules$tissue_ensembl_id %in% cis_eqtl_all] = TRUE
sum(modules$regulator)

mean(modules$regulator)


tab = mod_tab[, 1:2]
colnames(tab)[1] = "mod_id"


tab$n_regulators_TF_eSNP = table(modules$clust, modules$regulator)[, 2]



kda[kda$MODULE == 1, ]

kda_numbers = melt(table(kda$MODULE))
colnames(kda_numbers) = c("mod_id", "n_key_drivers")

tab = merge(tab, kda_numbers, all.x=TRUE)


tab$type[mod_tab$purity < 0.95] = "cross-tissue"
tab$type[mod_tab$purity >= 0.95] = "tissue-specific"

# write.csv(tab, "heritability/eQTL/module_eqtl_TF_KD.csv",
# 	row.names=FALSE)

write.csv(tab, "heritability/eQTL/module_eqtl_TF_KD_directed.csv",
	row.names=FALSE)



# tab$n_key_drivers / tab$n_regulators_TF_eSNP * 100

sum(tab$n_key_drivers, na.rm=TRUE) / sum(tab$n_regulators_TF_eSNP)
sum(tab$n_key_drivers[tab$type == "cross-tissue"], na.rm=TRUE) / sum(tab$n_regulators_TF_eSNP[tab$type == "cross-tissue"])
sum(tab$n_key_drivers[tab$type == "tissue-specific"], na.rm=TRUE) / sum(tab$n_regulators_TF_eSNP[tab$type == "tissue-specific"])



tab$n_key_drivers / tab$n_regulators_TF_eSNP * 100

hist(tab$n_key_drivers / tab$mod_size * 100, breaks=20)
