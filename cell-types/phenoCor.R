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

# Rename
names(cibersort_freq) = sapply(
	strsplit(names(cibersort_freq_matched), "[.]"),
	function(x) x[4]
)

tissue_col = brewer.pal(9, "Set1")  # tissue colors

cibersort_freq_matched = matchTrimCibersortFreq(cibersort_freq, pheno$starnet.ID)

# cor_thresh = 0.15
cor_thresh = 0.2
for (i in 1:length(cibersort_freq_matched)) {
	pdf(paste0("cell-types/plots/pheno-cor/", names(cibersort_freq_matched)[i], ".pdf"))
	cmat = cor(pheno, cibersort_freq_matched[[i]],
		use="pairwise.complete.obs")

	cmat = cmat[!rownames(cmat) %in% exclude_rows,]

	# Filter out all NA rows and cols
	cmat = cmat[
		apply(cmat, 1, function(row) any(!is.na(row))),
		apply(cmat, 2, function(col) any(!is.na(col)))
	]

	# Filter out minimum
	cmat = cmat[
		apply(cmat, 1, function(row) max(abs(row), na.rm=TRUE)) > cor_thresh,
		apply(cmat, 2, function(col) max(abs(col), na.rm=TRUE)) > cor_thresh
	]

	cmat[is.na(cmat)] = 0.0
	cmat = t(cmat)

	heatmap.2(cmat,
		trace="none",
		col=colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
		breaks=seq(-0.5, 0.5, length.out=101),  # cap of coloring 
		mar=c(10, 16),
		cexRow=0.7,
		main=names(cibersort_freq_matched)[i]
	)
	dev.off()
}

# fits_pheno = fitLinearEigenPheno(pheno, cibersort_freq_matched[[1]])

# 
pdf("cell-types/plots/cibersort_pheno.pdf", width=10.0)
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


plotCor = function(x, y, main="", ... ) {
	fit = lm(y~x)
	coef(fit)
	summary(fit)

 	cor_test = cor.test(x, y)

	plot(x, y, cex=0.7, pch=16,
		main=paste0(main, " p=", format(cor_test$p.value, digits=3, scientific=TRUE)),
		...)
	abline(coef=coef(fit), col="black", lwd=1.5)
}

# x = cibersort_freq_matched[[1]][["bone marrow:granulocyte macrophage progenitor"]]
# y = pheno$syntax_score


svg("cell-types/plots/cor_sel.svg", width=4)
par(mfrow=c(3, 2))
cell_type = "bone marrow:granulocyte macrophage progenitor"
phenotype = "syntax_score"
i = 1  # AOR
plotCor(
	cibersort_freq_matched[[i]][[cell_type]],
	pheno[[phenotype]],
	main=names(cibersort_freq_matched)[i],
	xlab=cell_type,
	col=tissue_col[i],
	ylab=phenotype
)

cell_type = "coronary artery:smooth muscle cell"
phenotype = "syntax_score"
i = 4
plotCor(
	cibersort_freq_matched[[i]][[cell_type]],
	pheno[[phenotype]],
	main=names(cibersort_freq_matched)[i],
	xlab=cell_type,
	col=tissue_col[i],
	ylab=phenotype
)

cell_type = "blood:dendritic cell, myeloid, immature"
cell_type = "blood:macrophage"
phenotype = "BMI"
i = 9
plotCor(
	cibersort_freq_matched[[i]][[cell_type]],
	pheno[[phenotype]],
	main=names(cibersort_freq_matched)[i],
	col=tissue_col[i],
	xlab=cell_type,
	ylab=phenotype
)

cell_type = "blood:monocyte"
phenotype = "Age"
i = 2
plotCor(
	cibersort_freq_matched[[i]][[cell_type]],
	pheno[[phenotype]],
	main=names(cibersort_freq_matched)[i],
	col=tissue_col[i],
	xlab=cell_type,
	ylab=phenotype
)

cell_type = "blood:natural killer cell"
phenotype = "Smoking.Years"
i = 2
plotCor(
	cibersort_freq_matched[[i]][[cell_type]],
	as.numeric(pheno[[phenotype]]),
	main=names(cibersort_freq_matched)[i],
	col=tissue_col[i],
	xlab=cell_type,
	ylab=phenotype
)

cell_type = "blood:neutrophil"
phenotype = "LDL"
i = 2
plotCor(
	cibersort_freq_matched[[i]][[cell_type]],
	as.numeric(pheno[[phenotype]]),
	main=names(cibersort_freq_matched)[i],
	col=tissue_col[i],
	xlab=cell_type,
	ylab=phenotype
)
dev.off()

