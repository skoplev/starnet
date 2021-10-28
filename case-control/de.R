# Differential expression analysis of case-control STARNET samples

rm(list=ls())

library(data.table)

library(DESeq2)
library(gplots)
library(RColorBrewer)
library(biomaRt)
library(edgeR)

# Load biomaRt database for annotating the genes
ensembl = useEnsembl(biomart="ensembl",
	dataset="hsapiens_gene_ensembl", version=89)

gene_map = getBM(
	attributes=c("ensembl_gene_id", "hgnc_symbol", "gene_biotype"),
	mart=ensembl)


setwd("/Users/sk/GoogleDrive/projects/STARNET/cross-tissue")


# Load and parse phenotype data
# ------------------------------------------------
pheno = fread(
	"~/GoogleDrive/projects/STARNET/phenotype/data/current/STARNET_main_phenotype_table.2017_12_03.tsv"
)
pheno$Sex[pheno$Sex == ""] = NA


# Metabolic features for covariate analysis
metabolic_feat = c(
	"fP-TG(mmol/l)",
	"HbA1c(%)",
	"Waist/Hip",
	"BMI(kg/m2)",
	"fP-LDL-Chol(mmol/l)",
	"fP-HDL-Chol(mmol/l)",
	"CRP(mg/l)",
	"P-Chol(mmol/l)")

# Drug treatment features for covariate analysis
treatment_feat = c(
	"β-blocker",
	"ThrombInhib",
	"LongActNitr",
	"ShortActNitr",
	"ACE-Inhib/ARB",
	"LoopDiuretic",
	"ThiazidDiuretic",
	"LipidLowerer",
	"Ca-ChannelBlocker",
	"Carvedilol",
	"AntiCoag",
	"OralAntiDiabetics",
	"Insulin")

# Convert treatment strings to lower case
for (feat in treatment_feat) {
	pheno[[feat]] = tolower(pheno[[feat]])
}

# Prints
data.frame(sort(table(unlist(pheno[, treatment_feat, with=FALSE])), decreasing=TRUE))

# Recode 
no_values = c("np", "n", "bo", "on", "when needed", "if needed")
yes_values = c("yez", "tes", "ys", "ues")
na_values = c("", "x", "<NA>")


for (feat in treatment_feat) {
	pheno[[feat]][pheno[[feat]] %in% no_values] = "no"
	pheno[[feat]][pheno[[feat]] %in% yes_values] = "yes"
	pheno[[feat]][pheno[[feat]] %in% na_values] = NA
	pheno[[feat]] = factor(pheno[[feat]])
}
	

# Impute drug response data to most commonly observed value
for (feat in treatment_feat) {
	message(feat, " impute #:", sum(is.na(pheno[[feat]])))
	pheno[[feat]][is.na(pheno[[feat]])] = names(which.max(table(pheno[[feat]])))  # most common
}



# Impute metabolic features.
for (feat in metabolic_feat) {
	message(feat, " impute: ", sum(is.na(pheno[[feat]])))
	pheno[[feat]][is.na(pheno[[feat]])] = median(pheno[[feat]], na.rm=TRUE)
}


# Load MDS genetic components from Ke Hao
genotype_mds = fread("case-control/data/genotype_MDS_hao/PCs.plink.mds.csv")
genotype_mds =genotype_mds[, -1]
colnames(genotype_mds)[1] = "starnet.ID"
colnames(genotype_mds)[2:5] = paste0("MDS_", colnames(genotype_mds)[2:5])

# Add to pheno table
pheno = merge(pheno, genotype_mds, all.x=TRUE)

# Impute missing MDS to median values
for (feat in paste0("MDS_C", 1:4)) {
	pheno[[feat]][is.na(pheno[[feat]])] = median(pheno[[feat]], na.rm=TRUE)
}

# sum(is.na(pheno$MDS_C1))

# Load RNA-seq gene counts
# -----------------------------------------------
file_paths = file.path("case-control/data/feature_counts",
	c(
		"gene_counts_AOR.txt",
		"gene_counts_LIV.txt",
		"gene_counts_VAF.txt",
		"gene_counts_SF.txt",
		"gene_counts_SKLM.txt"
	)
)

# load output table from featureCount
counts = lapply(file_paths, function(path) {
	fread(path)
})

count_col = 7
min_total_counts = 1e6

# Excluded samples
exclude_samples = lapply(counts, function(d) {
	mat = data.matrix(d[, count_col:ncol(d)])
	rownames(mat) = d$Geneid

	# Parse bam file names
	bam_file_names = colnames(mat)
	starnet_id = sapply(strsplit(bam_file_names, "[/_]+"), function(elem) {
		paste(elem[2:8], collapse="_")
	})
	colnames(mat) = starnet_id

	exclude = apply(mat, 2, sum) <= min_total_counts
	return(colnames(mat)[exclude])
})

exclude_samples = unlist(exclude_samples)
# write.table(exclude_samples, "case-control/data/exclude_samples.txt",
# 	quote=FALSE, col.names=FALSE, row.names=FALSE)

# apply(counts[[1]], 2, sum)

length(exclude_samples) / sum(sapply(counts, ncol))

# Parse into data matrices
count_mats = lapply(counts, function(d) {
	mat = data.matrix(d[, count_col:ncol(d)])
	rownames(mat) = d$Geneid

	# Parse bam file names
	bam_file_names = colnames(mat)
	starnet_id = sapply(strsplit(bam_file_names, "[/_]+"), function(elem) {
		paste(elem[2:8], collapse="_")
	})
	colnames(mat) = starnet_id


	include = apply(mat, 2, sum) > min_total_counts
	mat = mat[, include]

	return(mat)
})
names(count_mats) = sapply(strsplit(file_paths, "[/_.]+"), function(x) x[7])


# calculate RPKM for cutoffs based on absolute expression
# Equivalent to 10% samples cutoff
# Recommend RPKM > 1.0
rpkm_mats = lapply(1:length(count_mats), function(i) {
	rpkm(count_mats[[i]], gene.length=counts[[i]]$Length)
})

rpkm_90perc = lapply(rpkm_mats, function(mat) {
	apply(mat, 1, quantile, 0.9)
})


# Look at sample sizes
# sapply(count_mats, dim)

total_aligned_reads = lapply(count_mats, function(x) {
	apply(x, 2, sum)
})

sapply(total_aligned_reads, median)

# hist(total_aligned_reads[[1]])

# Get total case-control sample sizes
lapply(count_mats, function(mat) {
	message("First column: ", colnames(mat)[1])
	out = list()

	# Match phenotype table to included samples
	# starnet_id = sapply(strsplit(colnames(count_mats[[i]]), "_"), function(x) x[3])
	starnet_id = sapply(strsplit(colnames(mat), "_"), function(x) x[3])
	pheno_matched = pheno[match(starnet_id, pheno$starnet.ID), ]

	# Specify main phenotype to ensure correct sign, + == up in CAD
	pheno_matched$CAD.status = factor(pheno_matched$CAD.status,
		levels=c("control", "case"))

	return(table(pheno_matched$CAD.status))
})


renamePheno = function(names) gsub("/|-|%|\\(|\\)", ".", names)

# Annot is optional annotation data for output table
fitDiffModel = function(mat, pheno, formula, annot=NA, parallel=TRUE) {
	message("First column: ", colnames(mat)[1])
	message("Fitting: ", formula)
	out = list()

	# Match phenotype table to included samples
	# starnet_id = sapply(strsplit(colnames(count_mats[[i]]), "_"), function(x) x[3])
	starnet_id = sapply(strsplit(colnames(mat), "_"), function(x) x[3])
	pheno_matched = pheno[match(starnet_id, pheno$starnet.ID), ]
	pheno_matched = data.frame(pheno_matched)

	# Specify main phenotype to ensure correct sign, + == up in CAD
	pheno_matched$CAD.status = factor(pheno_matched$CAD.status,
		levels=c("control", "case"))

	out$pheno_matched = pheno_matched

	message("Imputing missing age for #samples: ", sum(is.na(pheno_matched$Age)))
	out$impute_age_samples = sum(is.na(pheno_matched$Age))
	pheno_matched$Age[is.na(pheno_matched$Age)] = median(pheno_matched$Age, na.rm=TRUE)

	out$impute_gender_samples = sum(is.na(pheno_matched$Sex))
	pheno_matched$Sex[is.na(pheno_matched$Sex)] = names(which.max(table(pheno_matched$Sex)))  # most common gender

	colnames(pheno_matched) = renamePheno(colnames(pheno_matched))

	# treatment_feat_sel = renamePheno(setdiff(treatment_feat, "ShortActNitr"))  # Removes all same value drug

	# DESeq2 differential analysis
	dds = DESeqDataSetFromMatrix(
		countData=mat,
		colData=pheno_matched,
		design=formula
	)

	# Fit model
	dds = DESeq(dds,
		parallel=parallel)

	# Calculate DE statistics
	res = results(dds,
		independentFiltering=FALSE,
		cooksCutoff=Inf,
		parallel=parallel)


	# Add gene annotations to table
	res = cbind(res,
		gene_map[match(rownames(res), gene_map$ensembl_gene_id), c("hgnc_symbol", "gene_biotype")]
	)

	if (!is.na(annot)) {
		res = cbind(res, annot)
	}

	res = res[order(res$pvalue),]

	out$size_factors = sizeFactors(dds)
	out$deg = res
	# Statistics
	out$covariate_stats = table(pheno_matched$Sex, pheno_matched$CAD.status)
	out$dds = dds

	gc()
	return(out)
}


# Fit DESeq2 models
# ----------------------------------------------

# # Original code for first DE analysis
# deseq = lapply(count_mats, function(mat) {
# 	message("First column: ", colnames(mat)[1])
# 	out = list()

# 	# Match phenotype table to included samples
# 	# starnet_id = sapply(strsplit(colnames(count_mats[[i]]), "_"), function(x) x[3])
# 	starnet_id = sapply(strsplit(colnames(mat), "_"), function(x) x[3])
# 	pheno_matched = pheno[match(starnet_id, pheno$starnet.ID), ]

# 	# Specify main phenotype to ensure correct sign, + == up in CAD
# 	pheno_matched$CAD.status = factor(pheno_matched$CAD.status,
# 		levels=c("control", "case"))

# 	out$pheno_matched = pheno_matched

# 	message("Imputing missing age for #samples: ", sum(is.na(pheno_matched$Age)))
# 	out$impute_age_samples = sum(is.na(pheno_matched$Age))
# 	pheno_matched$Age[is.na(pheno_matched$Age)] = median(pheno_matched$Age, na.rm=TRUE)

# 	out$impute_gender_samples = sum(is.na(pheno_matched$Sex))
# 	pheno_matched$Sex[is.na(pheno_matched$Sex)] = names(which.max(table(pheno_matched$Sex)))  # most common gender

# 	# DESeq2 differential analysis
# 	dds = DESeqDataSetFromMatrix(
# 		countData=mat,
# 		colData=pheno_matched,
# 		design=~Sex + Age + CAD.status
# 	)

# 	# Fit model
# 	dds = DESeq(dds,
# 		parallel=TRUE)

# 	# Calculate DE statistics
# 	res = results(dds,
# 		independentFiltering=FALSE,
# 		cooksCutoff=Inf,
# 		parallel=TRUE)


# 	# Add gene annotations to table
# 	res = cbind(res,
# 		gene_map[match(rownames(res), gene_map$ensembl_gene_id), c("hgnc_symbol", "gene_biotype")]
# 	)

# 	res = res[order(res$pvalue),]

# 	out$size_factors = sizeFactors(dds)
# 	out$deg = res
# 	# Statistics
# 	out$covariate_stats = table(pheno_matched$Sex, pheno_matched$CAD.status)
# 	out$dds = dds

# 	gc()
# 	return(out)
# })

# 

# deseq = lapply(count_mats, function(mat) {
# 	fitDiffModel(mat, pheno, formula(~Sex + Age + CAD.status))
# })

deseq = lapply(1:length(count_mats), function(i) {
	stopifnot(rownames(count_mats[[i]]) == names(rpkm_90perc[[i]]))
	fitDiffModel(count_mats[[i]],
		pheno,
		formula(~Sex + Age + CAD.status),
		parallel=FALSE,
		annot=data.frame(RPKM_q90=rpkm_90perc[[i]])
	)
})
names(deseq) = names(count_mats)



# Fit DESeq2 models adjusted for drug treatment
# ----------------------------------------------

treatment_feat_sel = renamePheno(setdiff(treatment_feat, "ShortActNitr"))  # Removes all same value drug
deseq_drug = lapply(count_mats, function(mat) {
	fitDiffModel(mat, pheno,
		formula=as.formula(paste(
			"~",
			paste(treatment_feat_sel, collapse="+"),
			"+ Sex + Age + CAD.status"))
	)
})



# Fit DESeq2 adjusting for metabolic parameters
# --------------------------------------------------

deseq_pheno = lapply(count_mats, function(mat) {
	fitDiffModel(mat, pheno, 
		formula=as.formula(paste("~",
			paste(renamePheno(metabolic_feat), collapse="+"),
			"+ Sex + Age + CAD.status"))
	)
})



# drug + pheno
treatment_feat_sel = renamePheno(setdiff(treatment_feat, "ShortActNitr"))  # Removes all same value drug
deseq_drug_pheno = lapply(count_mats, function(mat) {
	fitDiffModel(mat, pheno,
		formula=as.formula(paste("~",
			paste(c(renamePheno(metabolic_feat), treatment_feat_sel), collapse="+"),
			"+ Sex + Age + CAD.status"))
	)
})

# drug + pheno + genotypeMDS
mds = paste0("MDS_C", 1:4)
deseq_drug_pheno_mds = lapply(count_mats, function(mat) {
	fitDiffModel(mat, pheno,
		formula=as.formula(paste("~",
			paste(c(renamePheno(metabolic_feat), treatment_feat_sel, mds), collapse="+"),
			"+ Sex + Age + CAD.status"))
	)
})

# Save DESeq results
# save(deseq, deseq_drug, deseq_pheno, treatment_feat_sel, deseq_drug_pheno_mds, file="case-control/deseq_results.RData")


# Fit individual feature adjustments for single variable comparison signatures.
# not included in main run (DE2.RData file).
# -----------------------------------------------------------------------------
treatment_feat_sel = renamePheno(setdiff(treatment_feat, "ShortActNitr"))  # Removes all same value drug

deseq_treat_single = lapply(renamePheno(treatment_feat_sel), function(feat) {
	# Fit model per tissue
	lapply(count_mats, function(mat) {
		fit = fitDiffModel(mat, pheno,
			formula=as.formula(paste("~", feat, "+ Sex + Age + CAD.status"))
		)
		return(fit$deg)  # only DEG table for memory efficiency
	})
})
names(deseq_treat_single) = treatment_feat_sel
# save(deseq_treat_single, file="~/DataProjects/STARNET/DEG/deseq_treat_single.RData")

deseq_metabolic_single = lapply(renamePheno(metabolic_feat), function(feat) {
	# Fit model per tissue
	lapply(count_mats, function(mat) {
		fit = fitDiffModel(mat, pheno,
			formula=as.formula(paste("~", feat, "+ Sex + Age + CAD.status"))
		)
		return(fit$deg)  # only DEG table for memory efficiency
	})
})
names(deseq_metabolic_single) = metabolic_feat


# DEGs accounted for by individual covariates.
# -----------------------------------------
# Compare DEG per tissue, make signatures
tissue = "AOR"
feat = "ACE.Inhib.ARB"

deg_main = deseq[[tissue]]$deg
deg_adjust = deseq_treat_single[[feat]][[tissue]]

getAdjustGenes = function(deg_main, deg_adjust, fdr=0.01, fc_lim=log2(1.3)) {
	# match tables
	deg_adjust = deg_adjust[match(rownames(deg_adjust), rownames(deg_main)), ]

	# sig_main_idx = deg_main$padj < fdr & abs(deg_main$log2FoldChange) > fc_lim
	sig_main_idx = deg_main$padj < fdr & abs(deg_main$log2FoldChange) > fc_lim & deg_main$gene_biotype == "protein_coding"
	sum(sig_main_idx, na.rm=TRUE)

	sig_adjust_idx = sig_main_idx & deg_adjust$padj >= fdr  # not significant when accounted for feat
	sum(sig_adjust_idx, na.rm=TRUE)
	genes = unique(deg_main$hgnc_symbol[which(sig_adjust_idx)])

	genes = genes[genes != ""]
}


tissues = names(deseq)
for (tissue in tissues) {
	for (feat in names(deseq_treat_single)) {
		genes = getAdjustGenes(
			deg_main=deseq[[tissue]]$deg, 
			deg_adjust=deseq_treat_single[[feat]][[tissue]]
		)
		write(genes, paste0("case-control/signatures/single_adjust/", feat, "_", tissue, ".txt"))
	}
}

# write(genes, paste0("case-control/signatures/single_adjust/", feat, "_", tissue, ".txt"))




plotCmpDE = function(deg1, deg2, fdr=0.01, fc_lim=log2(1.3), prim_col="black", ...) {

	# Filter: significant genes in method 1 
	deg1 = deg1[which(abs(deg1$log2FoldChange) > fc_lim & deg1$padj < fdr), ]
	print(nrow(deg1))

	deg2 = deg2[match(rownames(deg1), rownames(deg2)), ]

	stopifnot(all(rownames(deg1) == rownames(deg2)))

	pts_col = rep("grey", nrow(deg1))
	pts_col[deg2$padj < fdr] = prim_col

	# plot(deg1$log2FoldChange, deg2$log2FoldChange)
	plot(-log10(deg1$padj), -log10(deg2$padj),
		cex=0.3,
		col=pts_col,
		pch=16,
		...)
	# abline(1, 1, col=brewer.pal(9, "Set1")[1])  # red identity line
	abline(1, 1, col="black")  # black identity line

	abline(h=-log10(fdr), lty=2, col="grey")  # line indicating reproducible FDR threshold after adjusting with method2

	# Number of transcripts that are not significant after adjusting for method2
	message(sum(deg2$padj > fdr, na.rm=TRUE))

	legend("bottomright",
		text.col=c(prim_col, "grey"),
		legend=c(
			sum(deg2$padj < fdr, na.rm=TRUE),  # reproduced
			sum(deg2$padj >= fdr, na.rm=TRUE)  # not reproduced
		),
		bty="n"
	)

	cor_test = cor.test(-log10(deg1$padj), -log10(deg2$padj))
	legend("topleft",
		paste0("R2=", format(cor_test$estimate^2, digits=4)),
		bty="n"
	)
}



tissues = names(deseq)

# pdf("case-control/plots/treatment_metabolic_adjustment.pdf", width=12, height=8)
pdf("case-control/plots/treatment_metabolic_adjustment_v2.pdf", width=12, height=11)
tis_cols = brewer.pal(9, "Set1")[c(1, 7, 4, 8, 3)]
names(tis_cols) = tissues
tis_cols = as.list(tis_cols)

par(mfrow=c(4, length(tissues)))

for(t in tissues) {
	plotCmpDE(deseq[[t]]$deg, deseq_pheno[[t]]$deg,
		xlab="-log10 p (BH, age+gender adj.)",
		ylab="-log10 p (+metabolic adj.)",
		prim_col=tis_cols[[t]],
		main=t)
}

for(t in tissues) {
	plotCmpDE(deseq[[t]]$deg, deseq_drug[[t]]$deg,
		xlab="-log10 p (BH, age+gender adj.)",
		ylab="-log10 p (+treatment adj.)",
		prim_col=tis_cols[[t]],
		main=t)
}

for(t in tissues) {
	plotCmpDE(deseq[[t]]$deg, deseq_drug_pheno[[t]]$deg,
		xlab="-log10 p (BH, age+gender adj.)",
		ylab="-log10 p (+metabolic+treatment adj.)",
		prim_col=tis_cols[[t]],
		main=t)
}

for(t in tissues) {
	plotCmpDE(deseq[[t]]$deg, deseq_drug_pheno_mds[[t]]$deg,
		xlab="-log10 p (BH, age+gender adj.)",
		ylab="-log10 p (+metabolic+treatment+genotype adj.)",
		prim_col=tis_cols[[t]],
		main=t)
}

dev.off()


# plotCmpDE(deseq$AOR$deg, deseq_drug$AOR$deg)
# plotCmpDE(deseq$AOR$deg, deseq_pheno$AOR$deg)
# plotCmpDE(deseq$AOR$deg, deseq_drug$SKLM$deg)
# plotCmpDE(deseq$AOR$deg, deseq_drug_pheno$AOR$deg)
# plotCmpDE(deseq$SKLM$deg, deseq_drug$SKLM$deg)



# Secreted protein annotations
sec = fread("~/GoogleDrive/projects/STARNET/STARNET-endocrine/from_marcus/secreted_proteins/uniprot_human_secreted_proteins.tab")
sec_symbols = unique(sec[["Gene names  (primary )"]])
sec_symbols = sec_symbols[sec_symbols != ""]

# Add secreted protein column
deseq = lapply(deseq, function(d) {
	d$deg$secreted = d$deg$hgnc_symbol %in% sec_symbols
	return(d)
})




# Write DESeq2 output tables
fdr = 0.01
for (i in 1:length(deseq)) {
	tissue = names(deseq)[i]
	deg = deseq[[i]]$deg

	write.csv(deg,
		file=paste0("case-control/data/deseq/deseq_full_", tissue, ".csv")
	)

	deg = deg[which(deg$padj < fdr), ]
	write.csv(deg,
		file=paste0("case-control/data/deseq/deseq_FDR001_", tissue, ".csv")
	)

	# Only secreted
	deg = deg[deg$secreted, ]
	write.csv(deg,
		file=paste0("case-control/data/deseq/deseq_secreted_", tissue, ".csv")
	)
}



# sum(d$deg$secreted)
# head(data.frame(deseq$AOR$deg), 50)
# head(data.frame(deseq$AOR$deg)[deseq$AOR$deg$secreted, ], 50)
# head(data.frame(deseq$LIV$deg)[deseq$LIV$deg$secreted, ], 100)


# Plot DEG statistics for selected genes
# -------------------------------------
# deloukas = read.csv("~/DataProjects/cross-tissue/GWAS/Deloukas/ng.csv")
deloukas = fread("~/DataProjects/cross-tissue/GWAS/Deloukas/ng.csv", skip=1)
CAD_genes = paste(deloukas$Loci_Nearest_Transcript, collapse="/")
CAD_genes = strsplit(CAD_genes, "/")[[1]]
CAD_genes = unique(CAD_genes)

# CAD_genes = "PCSK9"

CAD_DEG = list()
CAD_DEG$log2_fc = sapply(deseq, function(d) {
	d$deg$log2FoldChange[match(CAD_genes, d$deg$hgnc_symbol)]
})
rownames(CAD_DEG$log2_fc) = CAD_genes

CAD_DEG$padj = sapply(deseq, function(d) {
	d$deg$padj[match(CAD_genes, d$deg$hgnc_symbol)]
})
rownames(CAD_DEG$padj) = CAD_genes



library(devtools)  # for installing heatmap.3
# Load heatmap.3
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")

idx = apply(CAD_DEG$padj, 1, function(x) any(!is.na(x)))
# idx = idx & apply(CAD_DEG$padj, 1, function(p) any(p < 0.001))

mat = -log10(CAD_DEG$padj[idx, ])
# mat = CAD_DEG$log2_fc[idx, ]

# sec_idx = rownames(mat) %in% sec_symbols
# rownames(mat)[sec_idx] = paste0(rownames(mat)[sec_idx], "*")

rlab = data.frame(Secreted=rep("white", nrow(mat)))
rlab$Secreted = as.character(rlab$Secreted)
rlab$Secreted[which(rownames(mat) %in% sec_symbols)] = brewer.pal(9, "Set1")[2]

pdf("case-control/plots/deloukas_secreted_heatmap.pdf")
heatmap.3(mat,
	trace="none",
	col=colorRampPalette(brewer.pal(9, "YlOrRd"))(100),
	breaks=seq(0, 10, length.out=100 + 1),
	# col=colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
	# breaks=seq(-1, 1, length.out=100 + 1),
	cexRow=0.4,
	cexCo=0.8,
	mar=c(5, 30),
	main="Case-control",
	ylab="Deloukas CAD GWAS loci",
	# key.title="",
	KeyValueName="FDR",
	RowSideColors=t(rlab),
	density="density"
)
dev.off()


# padj =

# deseq$LIV$deg[deseq$LIV$deg$hgnc_symbol == "PCSK9", ]
# deseq$VAF$deg[deseq$VAF$deg$hgnc_symbol == "PCSK9", ]
# deseq$AOR$deg[deseq$AOR$deg$hgnc_symbol == "PCSK9", ]
# deseq$VAF$deg[deseq$VAF$deg$hgnc_symbol == "PCSK9", ]
# deseq$VAF$deg[deseq$VAF$deg$hgnc_symbol == "PCSK9", ]






# names(deseq[[1]])

# 
res = deseq$AOR$deg
write.csv(res$hgnc_symbol[which(res$padj < 0.01 & abs(res$log2FoldChange) > 1 & res$gene_biotype == "protein_coding")],
	"case-control/signatures/AOR01.txt",
	row.names=FALSE,
	col.names=FALSE,
	quote=FALSE
)


res = deseq$LIV$deg
write.csv(res$hgnc_symbol[which(res$padj < 0.01 & abs(res$log2FoldChange) > 1 & res$gene_biotype == "protein_coding")],
	"case-control/signatures/LIV01.txt",
	row.names=FALSE,
	col.names=FALSE,
	quote=FALSE
)

	
# Volcano plots.
# Note idx selection below.
# --------------------------------------------

#
# png("case-control/plots/case_control_volcano.png")
# pdf("case-control/plots/case_control_volcano.pdf", width=4.5)
# pdf("case-control/plots/case_control_volcano.pdf", width=8, height=12)
# pdf("case-control/plots/case_control_volcano2.pdf", width=8, height=6)
# pdf("case-control/plots/case_control_volcano3.pdf", width=7.5, height=5.5)
pdf("case-control/plots/case_control_volcano4.pdf", width=7.5, height=5.5)
# pdf("case-control/plots/case_control_volcano3_RPKMq90-05.pdf", width=7.5, height=5.5)
par(mfrow=c(2, 3))
fdr = 0.01
fc_lim = log2(1.3)  # +- 30%, same as signature definition


# colors = c(brewer.pal(8, "Set1"), brewer.pal(8, "Dark2"))
# colors = c("black", brewer.pal(8, "Set1"), brewer.pal(8, "Dark2"))
colors = c(brewer.pal(8, "Set1")[2:1], brewer.pal(8, "Dark2"))
# colors = c("grey", "black")

# biotypes = names(sort(table(deseq[[1]]$deg$gene_biotype), decreasing=TRUE)[1:14])

biotypes = c("non-coding", "protein_coding")
# biotypes = c("protein_coding", "lincRNA", "Non-coding")

# Most common biotypes
# biotypes = names(sort(table(deg$gene_biotype), decreasing=TRUE)[1:8])

# For plotting range
min_adj_p = sapply(deseq, function(d) {
	min(d$deg$padj, na.rm=TRUE)
})
min_adj_p = min(min_adj_p)

max_fc = sapply(deseq, function(d) {
	max(abs(d$deg$log2FoldChange), na.rm=TRUE)
})
max_fc = max(max_fc)

for (i in 1:length(deseq)) {
	deg = deseq[[i]]$deg
	tissue = names(deseq)[i]

	# col_pts = colors[match(deg$gene_biotype, biotypes)]

	col_pts = rep(colors[1], length(deg$gene_biotype))  # non-coding

	# Define coding RNAs
	col_pts[deg$gene_biotype == biotypes[2]] = colors[2]

	col_pts[is.na(col_pts)] = "black"
	col_pts[deg$padj > fdr] = "grey"

	# Significance
	idx = deg$padj < fdr & abs(deg$log2FoldChange) > fc_lim
	# idx = deg$padj < fdr & abs(deg$log2FoldChange) > fc_lim & deg$RPKM_q90 > 0.5  # RPKM q96


	# Coding and non-coding counts
	n_up_coding = sum(idx[deg$log2FoldChange > 0 & deg$gene_biotype == "protein_coding"], na.rm=TRUE)
	n_up_noncoding = sum(idx[deg$log2FoldChange > 0 & deg$gene_biotype != "protein_coding"], na.rm=TRUE)

	n_down_coding = sum(idx[deg$log2FoldChange < 0 & deg$gene_biotype == "protein_coding"], na.rm=TRUE)
	n_down_noncoding = sum(idx[deg$log2FoldChange < 0 & deg$gene_biotype != "protein_coding"], na.rm=TRUE)


	plot(deg$log2FoldChange[idx], -log10(deg$padj)[idx],
		pch=16,
		# cex=0.4,
		cex=0.2,
		main=tissue,
		col=col_pts,
		ylim=c(0, -log10(min_adj_p)),
		xlim=c(-max_fc, max_fc),
		ylab=expression("-log"[10] * " p (BH)"),
		xlab=expression("log"[2] * " FC")
	)

	# up-regulated
	legend("topright",
		legend=c(n_up_noncoding, n_up_coding),
		text.col=colors,
		bty="n"
	)

	# down-regulated
	legend("topleft",
		legend=c(n_down_noncoding, n_down_coding),
		text.col=colors,
		bty="n"
	)


	abline(h=-log10(fdr), lty=3, col="black")
}

plot(0, 0, type="n",
	xaxt="n",
	yaxt="n",
	bty="n",
	xlab="",
	ylab=""
	)
legend("topright", pch=16,
	xpd=TRUE,
	legend=biotypes,
	col=colors,
	cex=1.0
)
dev.off()


cleanGeneSymbols = function(symbols) {
	symbols[symbols == ""] = NA
	return(unique(na.omit(symbols)))
}

# Write DEG signatures to files
fdr = 0.01
fc_lim = log2(1.3)  # +- 30%, same as signature definition

for (i in 1:length(deseq)) {
	# Select table and tissue name
	deg = deseq[[i]]$deg
	tissue = names(deseq)[i]


	# significantly differentially expressed genes (boolean vectors)
	idx_up_coding = deg$padj < fdr &
		abs(deg$log2FoldChange) > fc_lim &
		deg$log2FoldChange > 0 &
		deg$gene_biotype == "protein_coding"
	sum(idx_up_coding, na.rm=TRUE)

	idx_up_noncoding = deg$padj < fdr &
		abs(deg$log2FoldChange) > fc_lim &
		deg$log2FoldChange > 0 &
		deg$gene_biotype != "protein_coding"
	sum(idx_up_noncoding, na.rm=TRUE)

	idx_down_coding = deg$padj < fdr &
		abs(deg$log2FoldChange) > fc_lim &
		deg$log2FoldChange < 0 &
		deg$gene_biotype == "protein_coding"
	sum(idx_down_coding, na.rm=TRUE)

	idx_down_noncoding = deg$padj < fdr &
		abs(deg$log2FoldChange) > fc_lim &
		deg$log2FoldChange < 0 &
		deg$gene_biotype != "protein_coding"
	sum(idx_down_noncoding, na.rm=TRUE)


	# Write signatures to files
	write.table(
		cleanGeneSymbols(deg$hgnc_symbol[idx_up_coding]),
		paste0("case-control/signatures/lists/DEG_", tissue, "_up_coding.txt"),
		quote=FALSE, row.names=FALSE, col.names=FALSE)

	write.table(
		cleanGeneSymbols(deg$hgnc_symbol[idx_up_noncoding]),
		paste0("case-control/signatures/lists/DEG_", tissue, "_up_noncoding.txt"),
		quote=FALSE, row.names=FALSE, col.names=FALSE)

	write.table(
		cleanGeneSymbols(deg$hgnc_symbol[idx_down_coding]),
		paste0("case-control/signatures/lists/DEG_", tissue, "_down_coding.txt"),
		quote=FALSE, row.names=FALSE, col.names=FALSE)

	write.table(
		cleanGeneSymbols(deg$hgnc_symbol[idx_down_noncoding]),
		paste0("case-control/signatures/lists/DEG_", tissue, "_down_noncoding.txt"),
		quote=FALSE, row.names=FALSE, col.names=FALSE)
}


# Counts of differentially expressed genes, accounting for expression level
# -----------------------------------------------------

fdr = 0.01
fc_lim = log2(1.3)  # +- 30%, same as signature definition
RPKM_q90_min = 0.5
base_mean_min = 50

deg_counts = sapply(deseq, function(de) {
	deg = de$deg

	idx = deg$padj < fdr & abs(deg$log2FoldChange) > fc_lim

	# Coding and non-coding counts
	n_up_coding = sum(idx[
			deg$log2FoldChange > 0 &
			deg$gene_biotype == "protein_coding"],
		na.rm=TRUE)
	n_up_coding_rpkm = sum(idx[
			deg$log2FoldChange > 0 &
			deg$gene_biotype == "protein_coding" &
			deg$RPKM_q90 > RPKM_q90_min],
		na.rm=TRUE)
	n_up_coding_base = sum(idx[
			deg$log2FoldChange > 0 &
			deg$gene_biotype == "protein_coding" &
			deg$baseMean > base_mean_min],
		na.rm=TRUE)

	n_up_noncoding = sum(idx[
			deg$log2FoldChange > 0 &
			deg$gene_biotype != "protein_coding"],
		na.rm=TRUE)
	n_up_noncoding_rpkm = sum(idx[
			deg$log2FoldChange > 0 &
			deg$gene_biotype != "protein_coding" &
			deg$RPKM_q90 > RPKM_q90_min],
		na.rm=TRUE)
	n_up_noncoding_base = sum(idx[
			deg$log2FoldChange > 0 &
			deg$gene_biotype != "protein_coding" &
			deg$baseMean > base_mean_min],
		na.rm=TRUE)

	n_down_coding = sum(idx[
			deg$log2FoldChange < 0 &
			deg$gene_biotype == "protein_coding"],
		na.rm=TRUE)
	n_down_coding_rpkm = sum(idx[
			deg$log2FoldChange < 0 &
			deg$gene_biotype == "protein_coding" &
			deg$RPKM_q90 > RPKM_q90_min],
		na.rm=TRUE)
	n_down_coding_base = sum(idx[
			deg$log2FoldChange < 0 &
			deg$gene_biotype == "protein_coding" &
			deg$baseMean > base_mean_min],
		na.rm=TRUE)

	n_down_noncoding = sum(idx[
			deg$log2FoldChange < 0 &
			deg$gene_biotype != "protein_coding"],
		na.rm=TRUE)
	n_down_noncoding_rpkm = sum(idx[
			deg$log2FoldChange < 0 &
			deg$gene_biotype != "protein_coding" &
			deg$RPKM_q90 > RPKM_q90_min],
		na.rm=TRUE)
	n_down_noncoding_base = sum(idx[
			deg$log2FoldChange < 0 &
			deg$gene_biotype != "protein_coding" &
			deg$baseMean > base_mean_min],
		na.rm=TRUE)

	return(data.frame(
		# Up-regulated
		n_up_coding,
		n_up_coding_rpkm,
		n_up_coding_base,

		n_up_noncoding,
		n_up_noncoding_rpkm,
		n_up_noncoding_base,

		# Down-regulatd
		n_down_coding,
		n_down_coding_rpkm,
		n_down_coding_base,

		n_down_noncoding,
		n_down_noncoding_rpkm,
		n_down_noncoding_base
	))
})

deg_counts = apply(deg_counts, 1, unlist)
deg_counts = as.data.frame(deg_counts)

deg_counts

# barplot(t(data.matrix(deg_counts)),
# 	las=2,
# 	beside=TRUE)

deg_frac_rpkm = data.frame(
	up_coding_rpkm=deg_counts$n_up_coding_rpkm / deg_counts$n_up_coding,
	up_noncoding_rpkm=deg_counts$n_up_noncoding_rpkm / deg_counts$n_up_noncoding,

	down_coding_rpkm=deg_counts$n_down_coding_rpkm / deg_counts$n_down_coding,
	down_noncoding_rpkm=deg_counts$n_down_noncoding_rpkm / deg_counts$n_down_noncoding,

	row.names=rownames(deg_counts)
)

deg_frac_base = data.frame(
	up_coding_base=deg_counts$n_up_coding_base / deg_counts$n_up_coding,
	up_noncoding_base=deg_counts$n_up_noncoding_base / deg_counts$n_up_noncoding,

	down_coding_base=deg_counts$n_down_coding_base / deg_counts$n_down_coding,
	down_noncoding_base=deg_counts$n_down_noncoding_base / deg_counts$n_down_noncoding,

	row.names=rownames(deg_counts)
)


pdf("case-control/plots/case_control_deg_filtering_barplot.pdf", width=7.5, height=5.5)

tis_cols = brewer.pal(9, "Set1")[c(1, 7, 4, 8, 3)]

par(mfrow=c(2, 1))
barplot(data.matrix(deg_frac_rpkm),
	main=paste0("RPKM q90 > ", RPKM_q90_min),
	ylab="DEG fraction",
	las=2,
	col=tis_cols,
	beside=TRUE)

barplot(data.matrix(deg_frac_base),
	main=paste0("base mean > ", base_mean_min),
	ylab="DEG fraction",
	las=2,
	col=tis_cols,
	beside=TRUE)
dev.off()

# sum(res$padj[res$gene_biotype == "protein_coding"] < 0.1, na.rm=TRUE)


# norm_counts = list()
# norm_counts = sweep(count_mats, 2, sizeFactors(dds), "/")

# Normalize count matrices by DESEq size factors
norm_counts = list()
for (i in 1:length(count_mats)) {
	norm_counts[[i]] = sweep(count_mats[[i]], 2, deseq[[i]]$size_factors, "/")
}
names(norm_counts) = names(count_mats)


# PCA
# -----------------------------------------------------------------
min_median_counts = 5
pca = lapply(norm_counts, function(mat) {
	include = apply(mat, 1, median) >= min_median_counts
	mat_filter = mat[include, ]

	pca_ = prcomp(t(mat_filter), center=TRUE, scale=TRUE)
	# pca_ = prcomp(t(log10(mat_filter+1)), center=TRUE, scale=TRUE)

	return(pca_)
})


stopifnot(all(rownames(norm_counts[[1]]) == rownames(norm_counts[[2]])))
stopifnot(all(rownames(norm_counts[[1]]) == rownames(norm_counts[[3]])))
stopifnot(all(rownames(norm_counts[[1]]) == rownames(norm_counts[[4]])))
stopifnot(all(rownames(norm_counts[[1]]) == rownames(norm_counts[[5]])))

norm_counts_all = Reduce(cbind, norm_counts)

include = apply(norm_counts_all, 1, median) >= min_median_counts

pca_all = prcomp(t(norm_counts_all[include, ]),
	center=TRUE,
	scale=TRUE
)

# pca_all = prcomp(t(log10(norm_counts_all[include, ] + 1)),
# 	center=TRUE,
# 	scale=TRUE
# )



pdf("case-control/plots/case_control_pca.pdf", width=8, height=12)
par(mfrow=c(3, 2))
for (i in 1:length(pca)) {

	col_pts = rep("black", ncol(norm_counts[[i]]))
	col_pts[deseq[[i]]$pheno_matched$CAD.status == "case"] = brewer.pal(9, "Set1")[1]

	pch_pts = rep(16, ncol(norm_counts[[i]]))
	pch_pts[deseq[[i]]$pheno_matched$CAD.status == "case"] = 17

	plot(pca[[i]]$x[, 1], pca[[i]]$x[, 2],
		pch=pch_pts,
		col=col_pts,
		xlab="PC1", ylab="PC2",
		main=names(count_mats)[i])
}

legend("topright",
	legend=c("ctrl", "case"),
	pch=c(16, 17),
	cex=1.0,
	col=c("black", brewer.pal(9, "Set1")[1])
)


tissue = sapply(strsplit(colnames(norm_counts_all), "_"), function(x) x[2])
starnet_id = sapply(strsplit(colnames(norm_counts_all), "_"), function(x) x[3])

pheno_matched = pheno[match(starnet_id, pheno$starnet.ID), ]

tissues = names(deseq)


col_pts = colors[match(tissue, tissues)]

pch_pts = rep(16, ncol(norm_counts_all))
pch_pts[pheno_matched$CAD.status == "case"] = 17


plot(pca_all$x[, 1], pca_all$x[, 2],
	cex=0.6,
	# pch=16,
	pch=pch_pts,
	col=col_pts,
	xlab="PC1", ylab="PC2"
)

legend("topleft", legend=tissues, pch=16, col=colors, cex=1.0)
dev.off()



# Check expression of genes
#------------------------------------


deseq$LIV$deg

# gene = "ENSG00000241985"
gene = "ENSG00000270641"
tissue = "AOR"



# "ENSG00000160613"

gene = "ENSG00000160613"
# gene = "ENSG00000280063"
# gene = "ENSG00000250644"
# gene = "ENSG00000264940"
# gene = "ENSG00000250644"
# gene = "ENSG00000070831"
# gene = "ENSG00000125870"
tissue = "LIV"


gene = "ENSG00000273108"
tissue = "VAF"

gene = "ENSG00000281490"
tissue = "SF"

gene = "ENSG00000279476"

gene = "ENSG00000160613"
tissue = "SKLM"



x = norm_counts[[tissue]][rownames(norm_counts[[tissue]]) == gene]
pheno_matched = deseq[[tissue]]$pheno_matched

expr = list()
expr$ctrl = x[pheno_matched$CAD.status == "control"]
expr$case = x[pheno_matched$CAD.status == "case"]

boxplot(expr, main=gene, ylab="Norm. expression")
points(jitter(rep(1, length(expr$ctrl))), expr$ctrl,
	pch=16)
points(jitter(rep(2, length(expr$case))), expr$case,
	pch=16)

t.test(expr$ctrl, expr$case)


sapply(expr, mean, na.rm=TRUE)


summary(deseq$AOR$dds)


lapply(deseq, function(d) {
	d$deg[which(abs(d$deg$log2FoldChange) > 1), ]
})

deseq$AOR$deg
deseq$LIV$deg
deseq$SKLM$deg

deseq$AOR$pheno_matched$CAD.status
deseq$LIV$pheno_matched$CAD.status


plotDispEsts(deseq$AOR$dds, main="Dispersion plot")



