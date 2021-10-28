library(data.table)
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(plyr)

rm(list=ls())

data_dir = "/Users/sk/DataProjects/cross-tissue"

emat_dir = file.path(data_dir, "STARNET/gene_exp/matrices")

setwd("/Users/sk/GoogleDrive/projects/STARNET/cross-tissue")


# Load phenotype data
# -------------------------------
pheno = fread(
	"~/GoogleDrive/projects/STARNET/phenotype/data/current/STARNET_main_phenotype_table.2017_12_03.tsv"
)

colnames(pheno) = make.names(colnames(pheno))  # syntactically valid names



# Load covariate table, copied from normalize.R
# ----------------------------------
# First covariate table
covar = fread(file.path(
	"~/GoogleDrive/projects/STARNET/phenotype/data/Oscar",
	"covariates.cases_and_controls.April_12_2016.txt"
))
covar = rename(covar, c("sample"="id"))

# Second covariate table
covar2 = fread(file.path(
	"~/GoogleDrive/projects/STARNET/phenotype/data/Oscar",
	"covariates.tbl"
))

# Merge covariate tables
covar_merged = merge(covar, covar2,
	by=c("id", "sex", "age"),
	all=TRUE)
# Fix discrepant read_lenghts. Default read_length to first covar table.
covar_merged = rename(covar_merged, c("read_length.x"="read_length"))
# Fill missing read_length fields from covar2
covar_merged$read_length[is.na(covar_merged$read_length)] = covar_merged$read_length.y[is.na(covar_merged$read_length)]
# Drop read_length.y
covar_merged = covar_merged[, names(covar_merged) != "read_length.y", with=FALSE]
# Drop incomplete fields
covar_merged = covar_merged[, !names(covar_merged) %in% c("subject", "batch"), with=FALSE]


covar_merged$starnet.ID = sapply(strsplit(covar_merged$id, "_"), function(x) x[2])
covar_merged$tissue = sapply(strsplit(covar_merged$id, "_"), function(x) x[1])
covar_merged$tissue = toupper(covar_merged$tissue)


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


# 

# expr_mats[[1]][1:10, 1:10]

# names(expr_mats)

phenotypes = c(
	"syntax_score",
	"DUKE",
	"fP-TG(mmol/l)",
	"HbA1c(%)",
	"Waist/Hip",
	"BMI(kg/m2)",
	"fP-LDL-Chol(mmol/l)",
	"fP-HDL-Chol(mmol/l)",
	"CRP(mg/l)",
	"P-Chol(mmol/l)")
phenotypes = make.names(phenotypes)


# Returns names most common occarance in vector
mostCommonValue = function(values) {
	if (all(is.na(values))) {
		return("unknown")  # dummy variable
	} else {
		return(names(sort(table(values), decreasing=TRUE))[1])
	}
}


deseq_results = lapply(phenotypes, function(feature) {
	message(feature)

	# Loop over all expression matrices
	results = list()
	for (i in 1:length(expr_mats)) {

		message(names(expr_mats)[i])


		mat = expr_mats[[i]]
		# mat = mat[1:100, ]  # test case

		# Merge phenotype and covariate table for selected tissue
		pheno_covar = merge(pheno,
			filter(covar_merged, tissue == !!names(expr_mats)[i]),
			by="starnet.ID")

		# Match patient IDs
		patient_ids = sapply(strsplit(colnames(mat), "_"), function(x) x[2])

		idx = match(patient_ids, pheno_covar$starnet.ID)

		if (sum(is.na(idx)) > 0) {
			message("WARNING: missing phenotype data for n=", sum(is.na(idx)))
		}

		pheno_matched = pheno_covar[idx, ]

		# Impute missing phenotype and covariate data
		pheno_matched[[feature]][is.na(pheno_matched[[feature]])] = median(pheno_matched[[feature]], na.rm=TRUE)
		pheno_matched$age[is.na(pheno_matched$age)] = median(pheno_matched$age, na.rm=TRUE)

		pheno_matched$lab[is.na(pheno_matched$lab)] = mostCommonValue(pheno_matched$lab)
		pheno_matched$sex[is.na(pheno_matched$sex)] = mostCommonValue(pheno_matched$sex)


		if (length(unique(pheno_matched$lab) == 1)) {
			warning("Skipping adjusting for lab.")
			dds = DESeqDataSetFromMatrix(
				countData=mat,
				colData=pheno_matched,
				design=as.formula(paste0("~age + sex + ", feature))
			)
		} else {
			dds = DESeqDataSetFromMatrix(
				countData=mat,
				colData=pheno_matched,
				design=as.formula(paste0("~lab + age + sex + ", feature))
			)
		}

		dds = DESeq(dds)


		res = results(dds)
		res = res[order(res$pvalue), ]

		# save results in list
		results[[names(expr_mats)[i]]] = res
	}

	# names(results) = names(expr_mats)
	return(results)
})
names(deseq_results) = phenotypes

save(deseq_results, file=file.path(data_dir, "STARNET/pheno_cor/deseq.RData"))

