rm(list=ls())

data_dir = "/Users/sk/DataProjects/cross-tissue"  # root of data directory
setwd("/Users/sk/Google Drive/projects/cross-tissue")

source("src/base.R")
source("src/parse/io.R")

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


plotPhenoCibersortCor(pheno_matched, cibersort_freq_matched,
	positive=TRUE,
	phenotype="BMI", 
	ciber_i=2, k=10, bar_col=brewer.pal(10, "Set1")[2]
)

plotPhenoCibersortCor(pheno_matched, cibersort_freq_matched,
	positive=FALSE,
	phenotype="BMI", ciber_i=9, k=5, bar_col=brewer.pal(9, "Set1")[2]
)

plotPhenoCibersortCor(pheno_matched, cibersort_freq_matched,
	positive=FALSE,
	phenotype="syntax_score", ciber_i=1, k=5, bar_col=brewer.pal(9, "Set1")[1]
)



plot(
	cibersort_freq_matched[[1]][["bone marrow:granulocyte macrophage progenitor"]],
	pheno_matched$syntax_score)


plot(
	cibersort_freq_matched[[9]][["blood:dendritic cell, myeloid, immature"]],
	pheno_matched$BMI)

plot(cibersort_freq_matched[[9]][["blood:macrophage"]], pheno_matched$BMI)


cor.test(cibersort_freq_matched[[1]][["bone marrow:granulocyte macrophage progenitor"]], pheno_matched$syntax_score, use="pairwise.complete.obs")
