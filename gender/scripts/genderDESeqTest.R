rm(list=ls())

library(dplyr)
library(DESeq2)


data_dir = "/Users/sk/DataProjects/cross-tissue"

emat_dir = file.path(data_dir, "STARNET/gene_exp/matrices")

setwd("/Users/sk/GoogleDrive/projects/STARNET/cross-tissue")

# Returns names most common occarance in vector
mostCommonValue = function(values) {
	if (all(is.na(values))) {
		return("unknown")  # dummy variable
	} else {
		return(names(sort(table(values), decreasing=TRUE))[1])
	}
}


# Load phenotype data
# -------------------------------
pheno = fread(
	"~/GoogleDrive/projects/STARNET/phenotype/data/current/STARNET_main_phenotype_table.2017_12_03.tsv"
)

colnames(pheno) = make.names(colnames(pheno))  # syntactically valid names

# Clean gender data
table(pheno$Sex)
pheno$Sex[pheno$Sex == ""] = NA
pheno$Sex[pheno$Sex == "patient nr male7female"] = NA
pheno$Sex[pheno$Sex == "male (sama mis 677)"] = "male"




# Load gene expression count matrices
# ----------------------------
# gene symbol, exp
expr_files = list.files(emat_dir, "*.mat")
expr_files = expr_files[c(-3, -4, -6)]

# Load all data matrices
expr_mats = lapply(expr_files, function(file_name) {
	message(file_name)
	d = fread(
		file.path(emat_dir, file_name))

	mat = d[,-1, with=FALSE]
	mat = data.matrix(mat)
	rownames(mat) = d$id

	return(mat)
})
names(expr_mats) = sapply(strsplit(expr_files, "[.]"), function(x) x[4])

# Fix tissue names
names(expr_mats) = recode(
	names(expr_mats),
	BLO="BLOOD",
	SUF="SF",
	SKM="SKLM")


# Gender association
results = list()
for (i in 1:length(expr_mats)) {
	mat = expr_mats[[i]]

	# Match patient IDs
	patient_ids = sapply(strsplit(colnames(mat), "_"), function(x) x[2])

	idx = match(patient_ids, pheno$starnet.ID)


	pheno_matched = pheno[idx, ]

	pheno_matched$Age[is.na(pheno_matched$Age)] = median(pheno_matched$Age, na.rm=TRUE)
	pheno_matched$Sex[is.na(pheno_matched$Sex)] = mostCommonValue(pheno_matched$Sex)


	dds = DESeqDataSetFromMatrix(
		countData=mat,
		colData=pheno_matched,
		design=as.formula("~Age + Sex")
	)

	dds = DESeq(dds)


	res = results(dds)
	res = res[order(res$pvalue), ]

	# save results in list
	results[[names(expr_mats)[i]]] = res

}

saveRDS(results, file="gender/tables/gender_deseq_tables_by_tissue.rds")
