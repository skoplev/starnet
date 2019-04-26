# Associations between CIBERSORT cell type fractions and phenotype
rm(list=ls())

setwd("~/GoogleDrive/projects/STARNET/cross-tissue")


scatterPlot = function(x, y,
	alternative="two.sided", corner="topleft", col="black", method="pearson", pch=16, ...) {
	cor_test = cor.test(x, y, alternative=alternative, method=method)

	plot(x, y,
		pch=pch,
		col=col,
		...)

	legend(corner, legend=c(
			paste0("r=", format(cor_test$estimate, digits=3)),
			paste0("P=", format(cor_test$p.value, digits=3))
		),
		bty="n"
	)

	fit = lm(y~x)
	abline(fit, col=col)
}

# Load phenotype data
pheno = fread("~/GoogleDrive/projects/STARNET/phenotype/data/current/STARNET_main_phenotype_table.2017_12_03.tsv")


# Load CIBERSORT
file_names = list.files("deconvolute/out_freq")
cibersort = lapply(file_names, function(file_name) {
	fread(file.path("deconvolute/out_freq", file_name))
})

names(cibersort) = sapply(strsplit(file_names, "[.]"), function(x) x[4])

# Match to starnet IDs
cibersort_matched = lapply(cibersort, function(tab) {
	id = sapply(strsplit(tab$V1, "_"), function(x) x[2])
	tab[match(pheno$starnet.ID, id), ]
})



cibersort$AOR

colnames(cibersort$AOR)[grep("macro", colnames(cibersort$AOR))]
colnames(cibersort$AOR)[grep("mono", colnames(cibersort$AOR))]
colnames(cibersort$AOR)[grep("artery", colnames(cibersort$AOR))]



features = colnames(cibersort$AOR)[grep("artery", colnames(cibersort$AOR))]

features = c(
	"blood:macrophage",
	"bone marrow:granulocyte macrophage progenitor",
	"blood:monocyte",
	"NA:monocyte"
)


par(mfrow=c(4, 4))
for (feat in features) {
	scatterPlot(cibersort_matched$AOR[[feat]] * 100, pheno$syntax_score,
		xlab=feat,
		ylab="SYNTAX",
		col=brewer.pal(9, "Set1")[2]
	)
}

for (feat in features) {
	scatterPlot(cibersort_matched$AOR[[feat]] * 100, pheno$DUKE,
		xlab=feat,
		ylab="DUKE",
		col=brewer.pal(9, "Set1")[1]
	)
}


features = colnames(cibersort_matched$AOR)[2:147]
cor_tests = lapply(features, function(feat) {
	# cor.test(cibersort_matched$AOR[[feat]], pheno$syntax_score)
	cor.test(cibersort_matched$AOR[[feat]], pheno[["BMI(kg/m2)"]])
})
names(cor_tests) = features

sort(sapply(cor_tests, function(x) x$p.value))
