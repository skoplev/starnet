rm(list=ls())

library(qvalue)
library(devtools)
library(compiler)
enableJIT(3)
library(data.table)

library(RcppEigen)  # fastLm
library(parallel)

# heatmap.3
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")


data_dir = "/Users/sk/DataProjects/cross-tissue"

setwd("/Users/sk/Google Drive/projects/STARNET/cross-tissue")
source("src/base.R")

# Load data
# ---------------------------------------------------------------

# STARNET phenotype data
# pheno = fread(file.path(
# 	"/Volumes/SANDY/phenotype_data",
# 	"STARNET_main_phenotype_table.cases.Feb_29_2016.tbl"
# ))

pheno = fread(
	"~/Google Drive/projects/STARNET/phenotype/data/current/STARNET_main_phenotype_table.2017_12_03.tsv"
)

# covar = fread(file.path(
# 	"/Volumes/SANDY/phenotype_data",
# 	"covariates.cases_and_controls.April_12_2016.txt"
# ))



# Load batch corrected expression data
load(file.path(data_dir, "STARNET/gene_exp_norm_batch/all.RData"))
expr_mats_batch = expr_mats_batch[!is.na(expr_mats_batch)]  # remove missing entries

# Rename loaded normalized gene expression matrices
names(expr_mats_batch) = sapply(
	strsplit(names(expr_mats_batch), "[.]"),
	function(x) x[4]
)

# Match phenotype data to selected gene expression matrix
pheno_match = lapply(expr_mats_batch, function(mat) {
	patient_ids = sapply(
		strsplit(colnames(mat), "_"),
		function(x) x[2]
	)

	pmatched = pheno[match(patient_ids, pheno$starnet.ID), ]

	return(pmatched)
})

# Fit linear regression models where SYNTAX score is the first 
syntax_regr = lapply(1:length(expr_mats_batch), function(i) {
	message("Fitting ", i, " out of ", length(expr_mats_batch))
	mat = expr_mats_batch[[i]]
	mat = as.matrix(mat)

	fits = apply(mat, 1, function(row) {
		fastLm(
			# row~syntax_score + Sex + BMI + Age,
			row~syntax_score + Sex + BMI + Age + LDL,
			data=pheno_match[[i]],
			na.action=na.omit)
	})

	return(fits)
})

syntax_pvals = lapply(syntax_regr, function(regr_fits) {
	pvals = sapply(regr_fits, function(fit) {
		summary(fit)$coefficients[2, 4]  # p-value of SYNTAX
	})
})

syntax_fdr = lapply(syntax_pvals, function(pvals) qvalue(pvals)$lfdr)

sapply(syntax_fdr, function(fdr) {
	sum(fdr < 0.5)
})

# i = 8
png("pheno/plots/account_gender_bmi_age_ldl.png", width=8, height=8, units="in", res=600)
par(mfrow=c(3, 3))
for (i in 1:9) {
	plot(-log10(syntax_cor[[i]]$pval), -log10(syntax_pvals[[i]]),
		cex=0.5,
		xlab=expression(log[10] * "p" * " (cor.)"),
		ylab=expression(log[10] * "p" * " (regr.)"),
		# ylab="Regr",
		col=rgb(0, 0, 0, 0.5),
		main=names(expr_mats_batch)[i])
	abline(a=0, b=1, col="red")
	abline(h=-log10(0.05), lty=3, col="grey")
	abline(v=-log10(0.05), lty=3, col="grey")
}
dev.off()

# syntax_pvals = sapply(syntax_regr[[1]], function(fit) {
# 	summary(fit)$coefficients[2, 4]
# })

sort(syntax_pvals[[4]])[1:20]


syntax_fdr = qvalue(syntax_pvals)$lfdr

sum(syntax_pvals < 0.001)
sum(syntax_fdr < 0.2)
sum(syntax_fdr < 0.5)

rownames(mat)[syntax_fdr < 0.4]

