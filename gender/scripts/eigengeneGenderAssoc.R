rm(list=ls())

library(data.table)

setwd("~/GoogleDrive/projects/STARNET/cross-tissue")

# Load phenotype data
# -----------------------------------------
# STARNET phenotype data
pheno = fread("~/GoogleDrive/projects/STARNET/phenotype/data/current/STARNET_main_phenotype_table.2017_12_03.tsv")

eig_mat = read.csv("co-expression/tables/eigengene_mat.csv")

rownames(eig_mat) = eig_mat[, 1]  # patient IDs in rows
eig_mat = eig_mat[, -1]

# Match pheno to eig_mat
pheno_matched = pheno[match(rownames(eig_mat), pheno$starnet.ID), ]

table(pheno_matched$Sex)

idx_male = pheno_matched$Sex == "male"
idx_female = pheno_matched$Sex == "female"


# Test eigen genes for gender association
tests_gender = lapply(1:ncol(eig_mat), function(i) {
	test = t.test(eig_mat[idx_male, i], eig_mat[idx_female, i])
	return(test)
})

stats_gender = data.frame(
	module=1:length(tests_gender),
	pvals=sapply(tests_gender, function(x) x$p.value)
)

stats_gender = stats_gender[order(stats_gender$pvals), ]



hist(pheno_matched$BMI)
idx_BMI_high = pheno_matched$BMI >= median(pheno_matched$BMI, na.rm=TRUE)
idx_BMI_low = pheno_matched$BMI < median(pheno_matched$BMI, na.rm=TRUE)


# Test eigen genes for gender association
tests_BMI = lapply(1:ncol(eig_mat), function(i) {
	test = t.test(eig_mat[idx_BMI_high, i], eig_mat[idx_BMI_low, i])
	return(test)
})

stats_BMI = data.frame(
	module=1:length(tests_BMI),
	pvals=sapply(tests_BMI, function(x) x$p.value)
)

stats_BMI = stats_BMI[order(stats_BMI$pvals), ]

pdf("gender/plots/eigenassoc_gender_BMI_ttest.pdf", width=6, height=8)
par(mfrow=c(2, 1))
plt_range = c(0, 55)
plot(-log10(stats_gender$pvals),
	xlab="Ranked co-expression modules",
	ylab="Gender (-log10 p)",
	ylim=plt_range,
	cex=0.5,
	pch=16)
abline(h=-log10(0.05), col="grey")
plot(-log10(stats_BMI$pvals),
	xlab="Ranked co-expression modules",
	ylab="High/low BMI (-log10 p)",
	ylim=plt_range,
	cex=0.5,
	pch=16)
abline(h=-log10(0.05), col="grey")
dev.off()
