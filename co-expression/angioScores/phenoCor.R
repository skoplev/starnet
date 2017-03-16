rm(list=ls())

library(RColorBrewer)
library(gplots)
library(qvalue)
library(data.table)


library(devtools)  # for installing heatmap.3
# Load heatmap.3
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")

library(compiler)
enableJIT(3)

data_dir = "/Users/sk/DataProjects/cross-tissue"  # root of data directory
setwd("/Users/sk/Google Drive/projects/cross-tissue")

source("src/models/regr.R")  # regression models
source("src/models/cor.R")



# Load cross-tissue modules
# ------------------------------
between = new.env()
load(file.path(data_dir, "modules/between_within-cross-tissue.RData"),
	between,
	verbose=TRUE)

complete = new.env()
load(file.path(data_dir, "modules/complete-cross-tissue.RData"),
	complete,
	verbose=TRUE)

# Parse module data
between = parseModuleData(between)
complete = parseModuleData(complete)

between_stats = countModuleTissueStat(between)
complete_stats = countModuleTissueStat(complete)


# Load phenotype data
# -----------------------------------

# STARNET phenotype data
pheno = fread(file.path(
	"/Volumes/SANDY/phenotype_data",
	"STARNET_main_phenotype_table.cases.Feb_29_2016.tbl"
))

# Brainshake phenotype data
brainshake = fread(file.path(
	"/Volumes/SANDY/phenotype_data",
	"tz.mat"
))


# Match phenotype to sample order
if (!all(between$patient_ids == complete$patient_ids)) {
	stop("Patient ID mismatch between module detection methods.")
} # else the same patient match can be used for both methods
patient_ids = between$patient_ids

# Match phenotype data tables to
pheno_matched = pheno[match(patient_ids, pheno$starnet.ID), ]
brainshake_matched = brainshake[match(patient_ids, brainshake$id), ]


# Test eigengene correlations with angio indicators
phenotypes = c("syntax_score", "ndv", "lesions", "DUKE")
between_pheno_cor = phenoCorTest(
	mat=between$bwnet$eigengenes,
	pheno_matched=pheno_matched,
	phenotypes=phenotypes
)

# Get pmat
pmat = sapply(between_pheno_cor, function(x) x$pval)
rownames(pmat) = 1:nrow(pmat)

# Nominally significant modules
sig_mods = apply(pmat, 1, min) < 0.05

pmat = pmat[sig_mods, ]

pdf("co-expression/angioScores/plots/angio_association_heatmap.pdf", height=16)
hmap = heatmap.2(-log10(pmat), trace="none",
	cexRow=1.2,
	# cexCol=0.8,
	key.title="",
	key.xlab=expression("-log"[10] * " p"),
	margins=c(10, 8),
	col=colorRampPalette(brewer.pal(9, "YlGnBu"))(100)
)
dev.off()

pdf("co-expression/angioScores/plots/angio_association_tissue_counts.pdf", height=16)

tissue_counts = between_stats$tissue_counts[, sig_mods]
colnames(tissue_counts) = which(sig_mods)

tissue_counts = tissue_counts[, hmap$rowInd]

colors = brewer.pal(9, "Set1")[-6]
barplot(tissue_counts, horiz=TRUE,
	las=2,
	xlim=c(0, 2000),
	col=colors)
legend("topright", legend=rownames(tissue_counts), pch=15, col=colors)
dev.off()

# Violin plot comparing cross
pdf("co-expression/angioScores/plots/angio_association_cross-tissue_fraction.pdf", height=4, width=3)
x = 1 - between_stats$purity[!sig_mods]
y = 1 - between_stats$purity[sig_mods]
par(bty="n")
vioplot(x, y,
	# main=t.test(x, y)$p.value,
	# col=brewer.pal(9, "Set1")
	col=brewer.pal(9, "Set1")[2],
	names=c("Non-angio", "Angio mods")
	# las=2
)
abline(h=0.05, col="grey", lty=3)
title(ylab="Cross-tissue fraction", main=t.test(x, y)$p.value)
dev.off()
