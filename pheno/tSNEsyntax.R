rm(list=ls())

library(tsne)
library(data.table)
library(RColorBrewer)
library(sva)  # ComBat
library(parallel)
library(fields)

library(compiler)
enableJIT(3)

data_dir = "/Users/sk/DataProjects/cross-tissue"

setwd("/Users/sk/Google Drive/projects/cross-tissue")
source("src/base.R")

# Load data
# ---------------------------------------------------------------

# STARNET phenotype data
pheno = fread(file.path(
	"/Volumes/SANDY/phenotype_data",
	"STARNET_main_phenotype_table.cases.Feb_29_2016.tbl"
))

covar = fread(file.path(
	"/Volumes/SANDY/phenotype_data",
	"covariates.cases_and_controls.April_12_2016.txt"
))


# expr_mats = lapply(list.files(file.path(data_dir, "STARNET/")))

load(file.path(data_dir, "STARNET/gene_exp_norm/all.RData"))

# Rename loaded normalized gene expression matrices
names(expr_mats_norm) = sapply(
	strsplit(names(expr_mats_norm), "[.]"),
	function(x) x[4]
)

# Rename columns to STARNET patient IDs
# expr_mats_norm = lapply(expr_mats_norm, function(mat) {
# 	colnames(mat) = sapply(
# 		strsplit(colnames(mat), "_"),
# 		function(x) x[2])
# 	return(mat)
# })


# Variance filter and batch correction
min_sd = 0.5
expr_mats_batch = lapply(expr_mats_norm, function(mat) {
	# Filter genes based on standard deviation
	mat = mat[apply(mat, 1, sd, na.rm=T) > min_sd, ]

	# Match phenotype data to selected gene expression matrix
	patient_ids = sapply(
		strsplit(colnames(mat), "_"),
		function(x) x[2]
	)

	pheno_matched = pheno[match(patient_ids, pheno$starnet.ID), ]
	covar_matched = covar[match(colnames(mat), covar$sample), ]

	# Correct batch effects
	# Model matrix taking into acount primariy covariates
	options(na.action="na.pass")  # keep NA rows in model matrix
	modcombat = model.matrix(~syntax_score + BMI + LDL + Age, data=pheno_matched)

	# Impute model input to median
	modcombat[is.na(modcombat[, 2]), 2] = median(modcombat[, 2], na.rm=T)  # syntax score
	modcombat[is.na(modcombat[, 3]), 3] = median(modcombat[, 3], na.rm=T)  # BMI
	modcombat[is.na(modcombat[, 4]), 4] = median(modcombat[, 4], na.rm=T)  # BMI
	modcombat[is.na(modcombat[, 5]), 5] = median(modcombat[, 5], na.rm=T)  # BMI

	# Identify batches
	batch = covar_matched$read_length
	if (all(is.na(batch))) {
		batch[is.na(batch)] = 100  # default read length if all are missing
	}
	batch[is.na(batch)] = median(batch, na.rm=T)  # 

	# Batch corrected expression matrix
	tryCatch({
		batch_mat = ComBat(mat,
			batch=batch,
			mod=modcombat  # model maintaining covariates
		)
		return(batch_mat)
	}, error=function(e){
		warning(e)
		# Use expression matrix as batch corrected data
		return(mat)
	})
})


# # read length and <90 patient ids
# batch=factor(
# 	paste(covar_matched$read_length, covar_matched$subject <= 89)
# ),

# # Detect additional batches
# clust = kmeans(t(batch_mat), 3)
# # clust = kmeans(t(batch_mat), 10)

# batch_mat = ComBat(batch_mat,
# 	batch=clust$clust,
# 	mod=modcombat  # model to maintain
# )


# Run tSNE for each expression matrix
embed = mclapply(expr_mats_batch, function(mat) {
	dmat = dist(t(mat))
	res = tsne(dmat, max_iter=2000, perplexity=10)
	return(res)
}, mc.cores=3)


# Selected phenotype for tSNE plots
phenotype = "syntax_score"
# phenotype = "BMI"
# phenotype = "ID"
# phenotype = "LDL"

pdf(paste0("pheno/plots/tSNE_", phenotype, ".pdf"), width=16)
par(mfrow=c(2, 5))
for (i in 1:length(expr_mats_batch)) {
	print(names(expr_mats_batch)[i])
	# i = 2
	mat = expr_mats_batch[[i]]
	res = embed[[i]]

	# Match phenotype data to selected gene expression matrix
	patient_ids = sapply(
		strsplit(colnames(mat), "_"),
		function(x) x[2]
	)

	pheno_matched = pheno[match(patient_ids, pheno$starnet.ID), ]
	covar_matched = covar[match(colnames(mat), covar$sample), ]

	# Pick color scheme
	if (phenotype == "syntax_score") {
		# # SYNTAX score
		# col = colorGradient(pheno_matched$syntax_score,
		# 	gradlim=c(0, 100),
		# 	# range(pheno$syntax_score, na.rm=T),
		# 	colors=rev(brewer.pal(9, "Spectral")))
		col = colorGradient(pheno_matched$syntax_score,
			gradlim=c(0, 100),
			# range(pheno$syntax_score, na.rm=T),
			colors=brewer.pal(9, "YlGnBu"))
	} else if (phenotype == "LDL") {
		# LDL 
		col = colorGradient(as.numeric(pheno_matched$LDL),
			gradlim=range(pheno$LDL),
			colors=rev(brewer.pal(9, "Spectral")))

	} else if (phenotype == "BMI") {
		# BMI
		col = colorGradient(as.numeric(pheno_matched$BMI),
			gradlim=range(pheno$BMI, na.rm=T),
			colors=rev(brewer.pal(9, "Spectral")))
	} else if (phenotype == "ID") {
		col = colorGradient(
			as.numeric(pheno_matched$starnet.ID),
			colors=brewer.pal(9, "Spectral"))
	} else {
		stop("Invalid phenotype")
	}
	# col=colorGradient(pheno_matched$DUKE, colors=rev(brewer.pal(9, "Spectral")))
	# col=colorGradient(pheno_matched$lesions, colors=rev(brewer.pal(9, "Spectral")))
	# col=colorGradient(as.numeric(pheno_matched$Smoking.Years), colors=rev(brewer.pal(9, "Spectral")))
	# col=colorGradient(as.numeric(covar_matched$read_length), colors=brewer.pal(9, "Spectral"))
	# col=colorGradient(as.numeric(pheno_matched$Age), colors=brewer.pal(9, "Spectral"))
	plot(res[,1], res[,2],
		col=col,
		pch=16,
		cex=1.5,
		main=names(expr_mats_batch)[i],
		xlab="", ylab=""
	)
	# text(res[,1], res[,2], labels=pheno_matched$starnet.ID, cex=0.5)
	# text(res[,1], res[,2], labels=pheno_matched$syntax_score, cex=0.5)
	# text(res[,1], res[,2], labels=pheno_matched$DUKE, cex=0.5)
}
dev.off()


# Plot color bars for selected tSNE plots
pdf("pheno/plots/legends.pdf", width=4.0)
par(mfrow=c(1, 2))
plotColorBar(
	colorGradient(seq(0, 100, length.out=100),
		colors=brewer.pal(9, "YlGnBu")),
	min=0, max=100, title="SYNTAX")

plotColorBar(
	colorGradient(
		seq(
			min(as.numeric(pheno$starnet.ID), na.rm=T),
			max(as.numeric(pheno$starnet.ID), na.rm=T),
			length.out=100),
	colors=rev(brewer.pal(9, "Spectral"))),
	min=0,
	max=max(as.numeric(pheno$starnet.ID), na.rm=T),
	# min=0,
	# max=100,
	title="STARNET ID")
dev.off()
