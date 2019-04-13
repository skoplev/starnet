# Analysis of olink CVD panel
# Pulled from annotateModules.R
#

rm(list=ls())

library(data.table)
library(RColorBrewer)

data_dir = "/Users/sk/DataProjects/cross-tissue"  # root of data directory
setwd("~/GoogleDrive/projects/STARNET/cross-tissue")

source("src/models/cor.R")


# Load co-expression data
between = new.env()
load(file.path(data_dir, "modules/between_within-cross-tissue.RData"),
	between,
	verbose=TRUE)


# Load olink data
olink = fread("~/GoogleDrive/projects/STARNET/olink/blood/STARNET_protein_uppsala/untitled.csv",
	# sep=",",
	header=TRUE)
olink[olink == "NAN"] = NA
# olink = data.matrix(olink)





# Load phenotype data
# -----------------------------------------
# STARNET phenotype data
pheno = fread("~/GoogleDrive/projects/STARNET/phenotype/data/current/STARNET_main_phenotype_table.2017_12_03.tsv")


# Match data
# ----------------------------------------
patient_ids = between$patient_ids

pheno_matched = pheno[match(patient_ids, pheno$starnet.ID), ]

olink_matched = olink[match(patient_ids, olink$NPX)]


# Global olink correlation with eigengenes
# -----------------------------------------

olink_cor = phenoCorTest(
	mat=between$bwnet$eigengenes,
	pheno_matched=olink_matched)


olink_pmat = sapply(olink_cor[-1], function(x) x$pval)
olink_qmat = sapply(olink_cor[-1], function(x) x$qval)

colnames(olink_pmat) = sapply(strsplit(colnames(olink_pmat), "_"), function(x) x[2])


# Screen for co-expresion module 78
# ------------------------------------------
i = 78
pdf(paste0("co-expression/annotate/plots/olink/mod", i, ".pdf"), width=10, height=5)
par(mfrow=c(1, 2))

olink_pmat[, colnames(olink_pmat) == "LEP"]

barplot(-log10(sort(olink_pmat[i, olink_qmat[i, ] < 0.1])),
	main=paste("Module", i),
	ylab=expression("-log" [10] * " p"),
	las=2)

plot(
	olink_matched[["189_LEP"]], 
	between$bwnet$eigengenes[, i],
	xlab="LEP",
	ylab=paste("Module", i),
	pch=16,
	cex=0.5,
)
cor.test(
	as.numeric(olink_matched[["189_LEP"]]), 
	between$bwnet$eigengenes[, i]
)
dev.off()

scatterPlot = function(x, y, col="black", ...) {
	plot(x, y,
		pch=16,
		cex=0.5,
		col=col,
		...
	)

	cor_test = cor.test(x, y)
	print(cor_test)

	fit = lm(y~x)

	abline(fit,
		col=col,
		lwd=2.0)

	legend("topleft",
		legend=c(
			paste0("r=", format(cor_test$estimate, digits=3)),
			paste0("P=", format(cor_test$p.value, digits=3))
		),
		cex=1.0,
		bty="n")
}


# Separate scatter plots
pdf("co-expression/annotate/plots/mod_98_78_assocation_scatter.pdf", height=3.5, width=6)
par(mfcol=c(2, 4), mar=c(4, 4, 2, 1))

i = 78
j = 98

colors = brewer.pal(9, "Set1")

scatterPlot(
	as.numeric(olink_matched[["189_LEP"]]), 
	between$bwnet$eigengenes[, i],
	col=colors[1],
	xlab="Plasma Leptin (log2 NPX)",
	ylab=paste("Eigengene", i)
)

scatterPlot(
	as.numeric(olink_matched[["189_LEP"]]), 
	between$bwnet$eigengenes[, j],
	col=colors[2],
	xlab="Plasma Leptin (log2 NPX)",
	ylab=paste("Eigengene", j)
)

# BMI
scatterPlot(
	as.numeric(pheno_matched[["BMI(kg/m2)"]]), 
	between$bwnet$eigengenes[, i],
	col=colors[1],
	xlab="BMI (kg/m2)",
	ylab=paste("Eigengene", i)
)

scatterPlot(
	as.numeric(pheno_matched[["BMI(kg/m2)"]]), 
	between$bwnet$eigengenes[, j],
	col=colors[2],
	xlab="BMI (kg/m2)",
	ylab=paste("Eigengene", j)
)

# HbAC1
scatterPlot(
	as.numeric(pheno_matched[["HbA1c(%)"]]), 
	between$bwnet$eigengenes[, i],
	col=colors[1],
	xlab="HbA1c(%)",
	xlim=c(4.5, 10),  # don't show <4 outlier
	ylab=paste("Eigengene", i)
)

scatterPlot(
	as.numeric(pheno_matched[["HbA1c(%)"]]), 
	between$bwnet$eigengenes[, j],
	col=colors[2],
	xlab="HbA1c (%)",
	xlim=c(4.5, 10),
	ylab=paste("Eigengene", j)
)

# Total cholesterol
scatterPlot(
	as.numeric(pheno_matched[["P-Chol(mmol/l)"]]), 
	between$bwnet$eigengenes[, i],
	col=colors[1],
	xlab="P-Chol (mmol/l)",
	ylab=paste("Eigengene", i)
)

scatterPlot(
	as.numeric(pheno_matched[["P-Chol(mmol/l)"]]), 
	between$bwnet$eigengenes[, j],
	col=colors[2],
	xlab="P-Chol (mmol/l)",
	ylab=paste("Eigengene", j)
)
dev.off()


# # LDL
# scatterPlot(
# 	as.numeric(pheno_matched[["fP-LDL-Chol(mmol/l)"]]), 
# 	between$bwnet$eigengenes[, i],
# 	col=colors[1],
# 	xlab="fP-LDL-Chol(mmol/l)",
# 	ylab=paste("Eigengene", i)
# )

# scatterPlot(
# 	as.numeric(pheno_matched[["fP-LDL-Chol(mmol/l)"]]), 
# 	between$bwnet$eigengenes[, j],
# 	col=colors[2],
# 	xlab="fP-LDL-Chol(mmol/l)",
# 	ylab=paste("Eigengene", j)
# )


# # HDL
# scatterPlot(
# 	as.numeric(pheno_matched[["fP-HDL-Chol(mmol/l)"]]), 
# 	between$bwnet$eigengenes[, i],
# 	col=colors[1],
# 	xlab="fP-HDL-Chol(mmol/l)",
# 	ylab=paste("Eigengene", i)
# )

# scatterPlot(
# 	as.numeric(pheno_matched[["fP-HDL-Chol(mmol/l)"]]), 
# 	between$bwnet$eigengenes[, j],
# 	col=colors[2],
# 	xlab="fP-HDL-Chol(mmol/l)",
# 	ylab=paste("Eigengene", j)
# )
