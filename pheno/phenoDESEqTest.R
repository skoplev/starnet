library(data.table)
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(plyr)

rm(list=ls())

data_dir = "/Users/sk/DataProjects/cross-tissue"

emat_dir = file.path(data_dir, "STARNET/gene_exp/matrices")

setwd("/Users/sk/GoogleDrive/projects/STARNET/cross-tissue")


# Load phenotype data
# -------------------------------
pheno = fread(
	"~/GoogleDrive/projects/STARNET/phenotype/data/current/STARNET_main_phenotype_table.2017_12_03.tsv"
)


# Load covariate table, copied from normalize.R
# ----------------------------------
# First covariate table
covar = fread(file.path(
	"~/GoogleDrive/projects/STARNET/phenotype/data/Oscar",
	"covariates.cases_and_controls.April_12_2016.txt"
))
covar = rename(covar, c("sample"="id"))

# Second covariate table
covar2 = fread(file.path(
	"~/GoogleDrive/projects/STARNET/phenotype/data/Oscar",
	"covariates.tbl"
))

# Merge covariate tables
covar_merged = merge(covar, covar2,
	by=c("id", "sex", "age"),
	all=TRUE)
# Fix discrepant read_lenghts. Default read_length to first covar table.
covar_merged = rename(covar_merged, c("read_length.x"="read_length"))
# Fill missing read_length fields from covar2
covar_merged$read_length[is.na(covar_merged$read_length)] = covar_merged$read_length.y[is.na(covar_merged$read_length)]
# Drop read_length.y
covar_merged = covar_merged[, names(covar_merged) != "read_length.y", with=FALSE]
# Drop incomplete fields
covar_merged = covar_merged[, !names(covar_merged) %in% c("subject", "batch"), with=FALSE]


covar_merged$starnet.ID = sapply(strsplit(covar_merged$id, "_"), function(x) x[2])
covar_merged$tissue = sapply(strsplit(covar_merged$id, "_"), function(x) x[1])


# Load gene expression count matrices
# ----------------------------
# gene symbol, exp
expr_files = list.files(emat_dir, "*.mat")
expr_files = expr_files[c(-3, -4, -6)]

# Load all data matrices
expr_mats = lapply(expr_files, function(file_name) {
	message(file_name)
	d = fread(
		file.path(emat_dir, file_name))

	mat = d[,-1, with=FALSE]
	mat = data.matrix(mat)
	rownames(mat) = d$id

	return(mat)
})
names(expr_mats) = sapply(strsplit(expr_files, "[.]"), function(x) x[4])

lapply(expr_mats, dim)

expr_mats[[1]][1:10, 1:10]



names(expr_mats)
# Select tissue
# i = 1  # AOR
i = 3  # LIV
# i = 7  # VAF
tissue = names(expr_mats)[i]


# Select feature
feature = "SYNTAX"
# feature = "DUKE"
# feature = "BMI(kg/m2)"
# feature = "fP-LDL-Chol(mmol/l)"

# Merge phenotype and covariate table for selected tissue
pheno_covar = merge(pheno,
	filter(covar_merged, tissue == !!names(expr_mats)[i]),
	by="starnet.ID")

# Match patient IDs
patient_ids = sapply(strsplit(colnames(expr_mats[[i]]), "_"), function(x) x[2])

idx = match(patient_ids, pheno_covar$starnet.ID)

if (sum(is.na(idx)) > 0) {
	message("WARNING: missing phenotype data for n=", sum(is.na(idx)))
}

pheno_matched = pheno_covar[idx, ]

colnames(pheno_matched) = make.names(colnames(pheno_matched))  # syntactically valid names

# Returns names most common occarance in vector
mostCommonValue = function(values) {
	names(sort(table(values), decreasing=TRUE))[1]
}

# Impute missing phenotype and covariate data
pheno_matched$syntax_score[is.na(pheno_matched$syntax_score)] = median(pheno_matched$syntax_score, na.rm=TRUE)  # impute to median
pheno_matched$DUKE[is.na(pheno_matched$DUKE)] = median(pheno_matched$DUKE, na.rm=TRUE)  # impute to median
pheno_matched$age[is.na(pheno_matched$age)] = median(pheno_matched$age, na.rm=TRUE)  # impute to median
pheno_matched$BMI.kg.m2.[is.na(pheno_matched$BMI.kg.m2.)] = median(pheno_matched$BMI.kg.m2., na.rm=TRUE)  # impute to median
pheno_matched$fP.LDL.Chol.mmol.l.[is.na(pheno_matched$fP.LDL.Chol.mmol.l.)] = median(pheno_matched$fP.LDL.Chol.mmol.l., na.rm=TRUE)  # impute to median

pheno_matched$lab[is.na(pheno_matched$lab)] = mostCommonValue(pheno_matched$lab)
pheno_matched$protocol[is.na(pheno_matched$protocol)] = mostCommonValue(pheno_matched$protocol)
pheno_matched$sex[is.na(pheno_matched$sex)] = mostCommonValue(pheno_matched$sex)


# Fit DESeq model
if (feature == "SYNTAX") {
	dds = DESeqDataSetFromMatrix(
		countData=expr_mats[[i]],
		colData=pheno_matched,
		design=~lab + age + sex + syntax_score
	)
} else if (feature == "DUKE") {
	dds = DESeqDataSetFromMatrix(
		countData=expr_mats[[i]],
		colData=pheno_matched,
		design=~lab + age + sex + DUKE
	)
} else if (feature == "BMI(kg/m2)") {
	dds = DESeqDataSetFromMatrix(
		countData=expr_mats[[i]],
		colData=pheno_matched,
		design=~lab + age + sex + BMI.kg.m2.
	)
} else if (feature == "fP-LDL-Chol(mmol/l)") {
	dds = DESeqDataSetFromMatrix(
		countData=expr_mats[[i]],
		colData=pheno_matched,
		design=~lab + age + sex + fP.LDL.Chol.mmol.l.
	)
}


dds = DESeq(dds)


res = results(dds)
res = res[order(res$pvalue), ]

head(data.frame(res), 50)

res$gene_symbol = sapply(strsplit(rownames(res), "_"), function(x) x[1])

res$gene_symbol_label = res$gene_symbol

# Show only labels for top 20 at most
res$gene_symbol_label[21:nrow(res)] = ""

# res$gene_symbol_label[res$padj > 0.1] = ""
res$gene_symbol_label[res$padj > 0.01] = ""


sum(res$padj < 0.01, na.rm=TRUE)
sum(res$padj < 0.05, na.rm=TRUE)
# summary(res$log2FoldChange)


head(data.frame(res), 20)

# pdf("pheno/plots/DESeq/AOR_SYNTAX_DESeq2.pdf")
# pdf("pheno/plots/DESeq/AOR_SYNTAX_DESeq2.pdf", height=12)

if (feature == "SYNTAX" & i == 1) {
	# AOR SYNTAX
	plot_height = 12
} else if (feature == "DUKE") {
	plot_height = 6
} else {
	# Default
	plot_height = 6
}


pdf(paste0("pheno/plots/DESeq/", tissue, "_", make.names(feature), "_DESeq2.pdf"), height=plot_height)
ggplot(data.frame(res), aes(x=log2FoldChange, y=-log10(res$padj), label=gene_symbol_label)) +
	ggtitle(tissue) +
	ylab("-log10 p (BH)") +
	xlab("SYNTAX log2FoldChange") +
	xlab(paste0(feature, " log2FoldChange")) +
	geom_point() +
	geom_text_repel() +
	theme_classic()
dev.off()


# plot(res$log2FoldChange, -log10(res$padj))

# j = which(rownames(expr_mats[[i]]) == "RAB3B_ENSG00000169213.6") 
# plot(expr_mats[[i]][j, ], pheno_matched$syntax_score)


sum(res$padj < 0.05, na.rm=TRUE)
sum(res$padj < 0.2, na.rm=TRUE)

library(limma)

if (feature == "SYNTAX" & tissue == "AOR") {
	# fdr_thresh = 0.2
	# idx = res$padj < fdr_thresh

	# deseq_genes = res$gene_symbol[which(idx)]

	deseq_genes_top50 = res$gene_symbol[1:50]
	cor_genes_top50 = fread("pheno/tables/aor_syntax50.tsv")$transcript_id

	all_genes = unique(c(deseq_genes_top50, cor_genes_top50))
	intersect(deseq_genes_top50, cor_genes_top50)

	bool_mat = data.frame(
		DESeq2=all_genes %in% deseq_genes_top50,
		Pearson_cor_norm=all_genes %in% cor_genes_top50)

	pdf(paste0("pheno/plots/DESeq/DESeq2_cor_comparison_venn_", tissue, "_", feature, ".pdf"))
	vennDiagram(
		vennCounts(bool_mat),
		circle.col=brewer.pal(9, "Set1"),
		main="SYNTAX AOR top-50 genes"
	)
	dev.off()

}


# # Compare to correlation statistics
# # -----------------------------------------------
# data_dir = "/Users/sk/DataProjects/cross-tissue"
# load(file.path(data_dir, "STARNET/pheno_cor/pheno_cor2.RData"))


# if (feature == "SYNTAX") {


	# pheno_cor$syntax_score$AOR$gene_symbol = sapply(strsplit(as.character(pheno_cor$syntax_score$AOR$transcript_id), "_"), function(x) x[1])
	# idx = pheno_cor$syntax_score$AOR$qval < fdr_thresh

# 	cor_genes = pheno_cor$syntax_score$AOR$gene_symbol[idx]

# }




