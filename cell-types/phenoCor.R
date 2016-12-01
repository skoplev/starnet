rm(list=ls())

library(RColorBrewer)
library(gplots)


data_dir = "/Users/sk/DataProjects/cross-tissue"  # root of data directory
setwd("/Users/sk/Google Drive/projects/cross-tissue")

source("src/base.R")
source("src/parse/io.R")
source("src/models/regr.R")  # regression models


# Load data
# ------------------------------------
# Load STARNET phenotype data
pheno = fread(file.path(
	"/Volumes/SANDY/phenotype_data",
	"STARNET_main_phenotype_table.cases.Feb_29_2016.tbl"
))

# Load CIBERSORT frequency data
cibersort_freq = loadCibersortFreq(file.path(data_dir, "CIBERSORT/out_freq"))

cibersort_freq_matched = matchTrimCibersortFreq(cibersort_freq, pheno$starnet.ID)


exclude_rows = c("starnet.ID")
cmat = cor(pheno, cibersort_freq_matched[[1]],
	use="pairwise.complete.obs")

cmat = cmat[!rownames(cmat) %in% exclude_rows,]

apply(cmat, 1, function(row) all(is.na(row)))


cmat[is.na(cmat)] = 0.0

heatmap.2(cmat,
	trace="none",
	col=colorRampPalette(brewer.pal(9, "RdBu"))(100))


pdf("cell-types/plots/cibersort_pheno.pdf",
	# width=4.5
	width=10.0
)

par(mfrow=c(2, 2))
plotPhenoCibersortCor(pheno, cibersort_freq_matched,
	phenotype="syntax_score", ciber_i=1, k=5, bar_col=brewer.pal(9, "Set1")[1]
)
plotPhenoCibersortCor(pheno, cibersort_freq_matched,
	positive=TRUE,
	phenotype="syntax_score", ciber_i=7, k=5, bar_col=brewer.pal(9, "Set1")[5]
)

plotPhenoCibersortCor(pheno, cibersort_freq_matched,
	phenotype="BMI", ciber_i=9, k=5, bar_col=brewer.pal(9, "Set1")[2]
)

plotPhenoCibersortCor(pheno, cibersort_freq_matched,
	positive=TRUE,
	phenotype="Age", 
	ciber_i=2, k=5, bar_col=brewer.pal(9, "Set1")[3]
)
dev.off()


plotPhenoCibersortCor(pheno, cibersort_freq_matched,
	positive=TRUE,
	phenotype="BMI", 
	ciber_i=2, k=10, bar_col=brewer.pal(10, "Set1")[2]
)

plotPhenoCibersortCor(pheno, cibersort_freq_matched,
	positive=FALSE,
	phenotype="BMI", ciber_i=9, k=5, bar_col=brewer.pal(9, "Set1")[2]
)

plotPhenoCibersortCor(pheno, cibersort_freq_matched,
	positive=FALSE,
	phenotype="syntax_score", ciber_i=1, k=5, bar_col=brewer.pal(9, "Set1")[1]
)


plotCor = function(x, y, ... ) {
	fit = lm(y~x)
	coef(fit)
	summary(fit)

 	cor_test = cor.test(x, y)

	plot(x, y, cex=0.7, pch=16,
		main=round(cor_test$p.value, 5),
		...)
	abline(coef=coef(fit), col="red")
}

# x = cibersort_freq_matched[[1]][["bone marrow:granulocyte macrophage progenitor"]]
# y = pheno$syntax_score

cell_type = "bone marrow:granulocyte macrophage progenitor"
phenotype = "syntax_score"
plotCor(
	cibersort_freq_matched[[1]][[cell_type]],
	pheno[[phenotype]],
	xlab=cell_type,
	ylab=phenotype
)

cell_type = "blood:dendritic cell, myeloid, immature"
cell_type = "blood:macrophage"
phenotype = "BMI"
plotCor(
	cibersort_freq_matched[[9]][[cell_type]],
	pheno[[phenotype]],
	xlab=cell_type,
	ylab=phenotype
)




plot(
	cibersort_freq_matched[[9]][["blood:dendritic cell, myeloid, immature"]],
	pheno$BMI)

plot(cibersort_freq_matched[[9]][["blood:macrophage"]], pheno$BMI)


cor.test(cibersort_freq_matched[[1]][["bone marrow:granulocyte macrophage progenitor"]], pheno_matched$syntax_score, use="pairwise.complete.obs")
