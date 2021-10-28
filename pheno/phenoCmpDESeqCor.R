library(dplyr)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)

rm(list=ls())

setwd("/Users/sk/GoogleDrive/projects/STARNET/cross-tissue")

data_dir = "/Users/sk/DataProjects/cross-tissue"

load(file.path(data_dir, "STARNET/pheno_cor/deseq.RData"), verbose=TRUE)

load(file.path(data_dir, "STARNET/pheno_cor/pheno_cor2.RData"), verbose=TRUE)
names(pheno_cor) = make.names(names(pheno_cor))

for (i in 1:length(pheno_cor)) {
	names(pheno_cor[[i]]) = recode(
		names(pheno_cor[[i]]),
		BLO="BLOOD",
		SKM="SKLM",
		SUF="SF"
	)
}


tissues = names(deseq_results[[1]])

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

head(pheno_cor[[1]][[1]])
head(deseq_results[[1]][[1]])

head(data.frame(deseq_results[["syntax_score"]][["VAF"]]), 20)
head(data.frame(deseq_results[["syntax_score"]][["LIV"]]), 20)

n_deseq = matrix(NA, nrow=length(tissues), ncol=length(phenotypes))
n_cor = matrix(NA, nrow=length(tissues), ncol=length(phenotypes))
n_overlap = matrix(NA, nrow=length(tissues), ncol=length(phenotypes))

rownames(n_deseq) = tissues
colnames(n_deseq) = phenotypes

rownames(n_cor) = tissues
colnames(n_cor) = phenotypes

rownames(n_overlap) = tissues
colnames(n_overlap) = phenotypes



fdr_cutoff = 0.01
r_cutoff = 0.2


for (i in 1:length(tissues)) {
	for (j in 1:length(phenotypes)) {
		tis = tissues[i]
		phe = phenotypes[j]

		# message(tis)
		# message(phe)

		cor_tab = pheno_cor[[phe]][[tis]]
		deseq_tab = deseq_results[[phe]][[tis]]
		deseq_tab = data.frame(deseq_tab)

		cor_genes = cor_tab$transcript_id[abs(cor_tab$cor) > r_cutoff]
		cor_genes = as.character(cor_genes)

		deseq_genes = rownames(deseq_tab)[deseq_tab$padj < fdr_cutoff]
		deseq_genes = na.omit(deseq_genes)
		deseq_genes = as.character(deseq_genes)


		n_deseq[i, j] = length(deseq_genes)
		n_cor[i, j] = length(cor_genes)
		n_overlap[i, j] = length(intersect(deseq_genes, cor_genes))
	}
}


# n_deseq
# n_cor

n_deseq_flat = melt(n_deseq, value.name="n_DESeq")
n_cor_flat = melt(n_cor, value.name="n_cor")
n_overlap_flat = melt(n_overlap, value.name="n_overlap")

d = merge(merge(n_deseq_flat, n_cor_flat), n_overlap_flat)
colnames(d)[1] = "tissue"
colnames(d)[2] = "phenotype"

write.csv(d, "pheno/tables/cmpCorDESeq.csv", row.names=FALSE)

# d$tissue

ggplot(d,
	aes(
		x=log10(n_DESeq + 1),
		y=log10(n_cor + 1),
		color=tissue, size=n_overlap,
		shape=phenotype)) +
	# stat_cor(method = "pearson") +
	# geom_smooth() +
	geom_point() +
	scale_shape_manual(values=1:nlevels(d$phenotype)) +
	# geom_smooth(method = "lm") +
	theme_classic()

plot(log10(n_deseq + 1), log10(n_cor + 1))
cor.test(log10(n_deseq + 1), log10(n_cor + 1))



pdf("pheno/plots/DESeq/pheno_cmp_DESeq_cor_all.pdf", width=2.5, height=8)

par(
	mfrow=c(nrow(n_deseq), 1),
	mar=c(1, 4, 1, 1)
)

for (i in 1:nrow(n_deseq)) {
	bar_counts = rbind(
		n_deseq[i, ] - n_overlap[i, ],
		n_overlap[i, ],
		n_cor[i, ] - n_overlap[i, ]
	)

	if (i == nrow(n_deseq)) {
		# last
		rownames(bar_counts) = c("DESeq2", "overlap", "Pearson cor")
	}

	# bar_counts = bar_counts[, 1:2]

	barplot(bar_counts,
		las=2,
		legend.text=rownames(bar_counts),
		col=c(
			brewer.pal(9, "Set1")[1],
			"black",
			brewer.pal(9, "Set1")[2]
		),
		ylab=paste0(rownames(n_deseq)[i], " genes")
	)
}
dev.off()