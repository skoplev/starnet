rm(list=ls())

library(qvalue)
library(devtools)
library(RColorBrewer)
library(data.table)

# heatmap.3
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")

data_dir = "/Users/sk/DataProjects/cross-tissue"

setwd("~/GoogleDrive/projects/STARNET/cross-tissue")
source("src/base.R")

# Load data
# ---------------------------------------------------------------

# STARNET phenotype data
# pheno = fread(file.path(
# 	"/Volumes/SANDY/phenotype_data",
# 	"STARNET_main_phenotype_table.cases.Feb_29_2016.tbl"
# ))

pheno = fread("~/Google Drive/projects/STARNET/phenotype/data/current/STARNET_main_phenotype_table.2017_12_03.tsv")


# covar = fread(file.path(
# 	"/Volumes/SANDY/phenotype_data",
# 	"covariates.cases_and_controls.April_12_2016.txt"
# ))


# Load batch corrected expression data
load(file.path(data_dir, "STARNET/gene_exp_norm_batch/all.RData"))
expr_mats_batch = expr_mats_batch[!is.na(expr_mats_batch)]  # remove missing entries

# # Rename loaded normalized gene expression matrices
# names(expr_mats_batch) = sapply(
# 	strsplit(names(expr_mats_batch), "[.]"),
# 	function(x) x[4]
# )

# Match phenotype data to selected gene expression matrix
pheno_match = lapply(expr_mats_batch, function(mat) {
	patient_ids = sapply(
		strsplit(colnames(mat), "_"),
		function(x) x[2]
	)

	pmatched = pheno[match(patient_ids, pheno$starnet.ID), ]

	return(pmatched)
})


# Test of phenotype correlations

# phenotypes = c("syntax_score", "LDL", "HDL")

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


pheno_cor = list()
for (phenotype in phenotypes) {
	message(phenotype)
	pheno_cor[[phenotype]] = lapply(1:length(expr_mats_batch), function(i) {
		mat = expr_mats_batch[[i]]
		mat = as.matrix(mat)

		# mat[mat < 0] = 0

		# Match phenotype data to selected gene expression matrix
		patient_ids = sapply(
			strsplit(colnames(mat), "_"),
			function(x) x[2]
		)

		pheno_matched = pheno[match(patient_ids, pheno$starnet.ID), ]

		# Correlation tests
		cor_tests = lapply(1:nrow(mat), function(k) {
			tryCatch({
				cor_test = cor.test(mat[k,], pheno_matched[[phenotype]], method="pearson")
				return(cor_test)
			}, error=function(e) {
				# Null
				cor_test = cor.test(c(0, 0, 0), c(0, 0, 0))
				return(cor_test)
			})
		})

		cor_pvals = sapply(cor_tests, function(x) x$p.value)
		cor_coef = sapply(cor_tests, function(x) x$estimate)

		# qvalue, local correction for multiple hypothesis testing
		cor_qvals = qvalue(p=cor_pvals)

		tab = data.frame(
			transcript_id=rownames(mat),
			cor=cor_coef,
			# cor_tests=cor_tests,
			pval=cor_pvals,
			qval=cor_qvals$lfdr
		)
		return(tab)
	})

	names(pheno_cor[[phenotype]]) = names(expr_mats_batch)

	# Order correlation tables based on p-values
	pheno_cor[[phenotype]] = lapply(pheno_cor[[phenotype]], function(x) {
		x[order(x$pval), ]
	})
}

save(pheno_cor, file=file.path(data_dir, "STARNET/pheno_cor/pheno_cor2.RData"))



# Density plots of distributions
# ------------------------------------------------

pdf("pheno/plots/pheno_cor_densities.pdf", width=8.0, height=3.5)
tissue_order = c("AOR", "MAM", "LIV", "VAF", "SUF", "SKM", "BLO")
tissue_col = brewer.pal(9, "Set1")[c(1, 5, 7, 4, 8, 3, 2)]

par(mfcol=c(2, 5), mar=c(4, 1.5, 2, 1.5))
for (i in 1:length(pheno_cor)) {

	# tissue_cor = pheno_cor[[i]][!names(pheno_cor[[i]]) %in% exclude_tissues]

	tissue_cor = pheno_cor[[i]][match(tissue_order, names(pheno_cor[[i]]))]

	densities = lapply(tissue_cor, function(x) {
		density(x$cor)
	})

	min_max_val = sapply(densities, function(d) c(min(d$x), max(d$x)))
	density_height = sapply(densities, function(d) c(max(d$y)))

	cor_range = range(min_max_val)

	delta = 2.0

	plot(0, 0, type="n",
		xlim=range(min_max_val),
		ylim=c(0, max(density_height) + delta*length(densities)),
		bty="n", yaxt="n",
		xlab=paste(names(pheno_cor)[i], "cor."),
		ylab=""
	)

	for (j in length(densities):1) {
		x = densities[[j]]$x
		y = densities[[j]]$y

		polygon(x, y + delta * (j - 1),
			col=addAlpha(tissue_col[j], 0.4),
			border=tissue_col[j])

		# Endpoints
		points(x[1],
			y[1] + delta * (j - 1),
			pch=16, col=tissue_col[j]
		)

		points(x[length(x)],
			y[length(y)] + delta * (j - 1),
			pch=16, col=tissue_col[j]
		)
	}
}
dev.off()



# Correlations based on permutations
# -------------------------------------------------------
K = 1000  # number of permutions

cor_coef = list()
cor_coef_perm = list()
for (i in 1:length(expr_mats_batch)) {
	message("Permutations test for: ", names(expr_mats_batch)[i])

	mat = expr_mats_batch[[i]]
	mat = as.matrix(mat)

	# Match phenotype data to selected gene expression matrix
	patient_ids = sapply(
		strsplit(colnames(mat), "_"),
		function(x) x[2]
	)

	pheno_matched = pheno[match(patient_ids, pheno$starnet.ID), ]

	cor_coef[[i]] = cor(t(mat), as.matrix(pheno_matched$syntax_score), use="pairwise.complete")

	# Permute phenotypes
	pheno_perm = lapply(1:K, function(k) {
		pheno_matched[sample(1:nrow(pheno_matched)), ]
	})

	# Permute test
	cor_coef_perm[[i]] = matrix(NA, nrow=nrow(mat), ncol=K)
	for (k in 1:K) {
		if (k %% 50 == 0) {
			message("	Perm iter: ", k)
		}

		cor_coef_perm[[i]][, k] = cor(t(mat), as.matrix(pheno_perm[[k]]$syntax_score), use="pairwise.complete")
	}
}
names(cor_coef) = names(expr_mats_batch)
names(cor_coef_perm) = names(expr_mats_batch)


# Specify tissue order for plotting
tissue_order = c("AOR", "MAM", "BLO", "SUF", "VAF", "LIV", "SKM", "MAC", "FOC")
tissue_idx = match(tissue_order, names(expr_mats_batch))

tissue_cols = c(
	brewer.pal(9, "Set1")[1:9 != 6],
	brewer.pal(8, "Dark2")
)


# Count significant correlations at
sig_cor = sapply(syntax_cor, function(x) {
	# sum(x$qval < 0.2)
	sum(x$pval < 0.001)
})
names(sig_cor) = names(expr_mats_batch)

# Write as tables

# AOR
aor_syntax_tab = syntax_cor[[which(names(syntax_cor) == "AOR")]]
aor_syntax_tab$transcript_id = sapply(
	strsplit(as.character(aor_syntax_tab$transcript_id), "_"),
	function(x) x[1])

write.table(head(aor_syntax_tab, 50),
	file="pheno/tables/aor_syntax50.tsv",
	row.names=FALSE, quote=FALSE, sep="\t")

# SUF
suf_syntax_tab = syntax_cor[[which(names(syntax_cor) == "SUF")]]
suf_syntax_tab$transcript_id = sapply(
	strsplit(as.character(suf_syntax_tab$transcript_id), "_"),
	function(x) x[1])

write.table(head(suf_syntax_tab, 50),
	file="pheno/tables/suf_syntax50.tsv",
	row.names=FALSE, quote=FALSE, sep="\t")


# LIV
liv_syntax_tab = syntax_cor[[which(names(syntax_cor) == "LIV")]]
liv_syntax_tab$transcript_id = sapply(
	strsplit(as.character(liv_syntax_tab$transcript_id), "_"),
	function(x) x[1])

write.table(head(liv_syntax_tab, 50),
	file="pheno/tables/liv_syntax50.tsv",
	row.names=FALSE, quote=FALSE, sep="\t")



pdf("pheno/plots/SYNTAX_cor_barplot.pdf", width=3.5, height=4)
barplot(sig_cor[tissue_idx],
	# sig_cor[!names(sig_cor) %in% c("COR", "FOC", "MAC")],
	main="p<0.001",
	# col=brewer.pal(9, "Set1"),
	col=tissue_cols,
	las=2,
	ylab="SYNTAX correlated RNAs")
dev.off()


# Count number of samples with SYNTAX score
nsamples = sapply(expr_mats_batch, function(mat) {
	mat = as.matrix(mat)
	# Match phenotype data to selected gene expression matrix
	patient_ids = sapply(
		strsplit(colnames(mat), "_"),
		function(x) x[2]
	)

	pheno_matched = pheno[match(patient_ids, pheno$starnet.ID), ]

	return(sum(!is.na(pheno_matched$syntax_score)))
})

# Wilcox test, comparing the correlation distribution with empirical null distribution
# wilcox = list()
# for (i in 1:length(cor_coef)) {
# 	message(i, " of ", length(cor_coef))
# 	wilcox[[i]] = wilcox.test(cor_coef[[i]], as.vector(cor_coef_perm[[i]]))
# }

# Seperate statistic for each random sample
wilcox_ensemble = list()
for (i in 1:length(cor_coef)) {
	message(i, " of ", length(cor_coef))
	wilcox_ensemble[[i]] = apply(cor_coef_perm[[i]], 2, function(rand_sample) {
		wilcox.test(cor_coef[[i]], rand_sample)
	})
}

# mean(sapply(wilcox_ensemble[[i]], function(x) x$p.value))

# plot(-log10(sort(sapply(wilcox_ensemble[[1]], function(x) x$p.value))))
# lines(-log10(sort(sapply(wilcox_ensemble[[2]], function(x) x$p.value))))

# sort(sapply(wilcox_ensemble[[i]], function(x) x$p.value))


pdf("pheno/plots/syntax_cor_distr.pdf", height=8, width=16)
# svg("pheno/plots/syntax_cor_distr.svg", height=5)
par(mfcol=c(3, 3), mar=c(8, 9, 4, 4))
# for (i in 1:length(cor_coef)) {
for (i in 1:length(tissue_idx)) {
	idx = tissue_idx[i]
	mean_wilcox_p = mean(sapply(wilcox_ensemble[[idx]], function(x) x$p.value))

	# par(lwd=0.5)
	histo = hist(cor_coef[[idx]], breaks=100,
		main=paste0(
			names(cor_coef)[idx],
			" n=", nsamples[idx],
			# " p=", format(wilcox[[idx]]$p.value, digits=3, scientific=TRUE)
			# " p=", format(mean_wilcox_p, digits=3, scientific=TRUE)
			" p=", format(mean_wilcox_p, digits=3)
			),
		cex.main=2.5,
		xlab="", ylab="",
		prob=TRUE,
		border=tissue_cols[i],
		xpd=TRUE,
		col=tissue_cols[i])

	# lines(c(histo$breaks, max(histo$breaks)), c(0, histo$density, 0), type="S")
	# barplot(histo$density, histo$breaks, space=0, col=tissue_cols[i], border=NA)

	lines(density(cor_coef_perm[[idx]]), 
		lwd=2.5,
		# xpd=TRUE,
		# col=rgb(0, 0, 0, 0.5)
	)
}
dev.off()


# Calculate empirical p-values for each gene

# # Gene-specific counts
# indicator_mat = sweep(abs(cor_coef_perm), 1, abs(cor_coef), "<")
# empi_counts = apply(indicator_mat, 1, sum)
# empi_pval = empi_counts / K

# cor_coef_perm_abs = abs(cor_coef_perm)
# # COunts for all genes
# empi_counts = sapply(cor_coef, function(stat) {
# 	sum(abs(stat) > cor_coef_perm_abs)
# })

# empi_pval = empi_counts / (K * nrow(mat))


# load(file=file.path(data_dir, "STARNET/pheno_cor/syntax_cor.RData"))



cors_to_plot = data.frame(
	tissue=c("AOR", "BLO", "LIV", "VAF", "FOC", "FOC"),
	transcript=c("FAM13C_ENSG00000148541.12", "GMEB1_ENSG00000162419.12", "INA_ENSG00000148798.9", "WFDC6_ENSG00000243543.8", "TP53BP2_ENSG00000143514.16", "PINX1_ENSG00000254093.8"),
	color=brewer.pal(9, "Set1")[c(1, 2, 3, 7, 9, 9)]
)

par(mfrow=c(3, 2))
for (i in 1:nrow(cors_to_plot)) {
	# kth expression data entry
	k = which(names(expr_mats_batch) == cors_to_plot$tissue[i])

	mat = expr_mats_batch[[k]]
	mat = as.matrix(mat)

	x = mat[rownames(mat) == cors_to_plot$transcript[i], ]
	y = pheno_match[[k]]$syntax_score

	complete = complete.cases(x, y)
	x = x[complete]
	y = y[complete]

	plot(x, y,
		col=as.character(cors_to_plot$color[i]),
		# col="#E41A1C",
		xlab=cors_to_plot$transcript[i],
		ylab="SYNTAX",
		pch=16, cex=0.6
	)

	# Fit line
	fit = lm(y~x)
	abline(coef(fit))
}

i = 1  # tissue
mat = expr_mats_batch[[i]]
mat = as.matrix(mat)


# Standardize expression matrix
mat = t(scale(t(mat)))


# Match phenotype data to selected gene expression matrix
patient_ids = sapply(
	strsplit(colnames(mat), "_"),
	function(x) x[2]
)

pheno_matched = pheno[match(patient_ids, pheno$starnet.ID), ]



# sig_transcripts = syntax_cor[[i]]$transcript_id[syntax_cor[[i]]$qval < 0.5]
sig_transcripts = syntax_cor[[i]]$transcript_id[syntax_cor[[i]]$qval < 0.2]
# sig_transcripts = syntax_cor[[i]]$transcript_id[syntax_cor[[i]]$pval < 0.001]

pdf("pheno/plots/cor_AOR_FDR02.pdf", height=8, width=12)
# rownames(mat)[cor_qvals$lfdr < 0.2]

submat = mat[rownames(mat) %in% sig_transcripts, ]
rownames(submat) = sapply(strsplit(rownames(submat), "_"), function(x) x[1])

rlab = colorGradient(pheno_matched$syntax_score,
	gradlim=c(0, 100),
	# colors=brewer.pal(9, "YlGnBu"))
	colors=brewer.pal(9, "Blues"))
rlab = as.matrix(rlab)
colnames(rlab) = "SYNTAX"

heatmap.3(
	submat,
	trace="none",
	col=colorRampPalette(rev(brewer.pal(9, "Spectral")))(100),
	breaks=seq(-3, 3, length.out=101),  # cap of coloring 
	ColSideColors=rlab,
	ColSideColorsSize=1,
	margins=c(6, 12),
	keysize=0.9,
	# cexRow=0.4,
	# KeyValueName=expression("-log"[10] * " p"),
	KeyValueName=expression("z-score"),
	xlab="Patients"
)
dev.off()
