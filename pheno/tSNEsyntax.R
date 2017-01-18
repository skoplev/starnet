rm(list=ls())

library(tsne)
library(amap)
library(data.table)
library(RColorBrewer)
library(sva)  # ComBat
library(parallel)
library(fields)
library(readxl)
library(dendextend)
library(plyr)

library(compiler)
enableJIT(3)

data_dir = "/Users/sk/DataProjects/cross-tissue"

setwd("/Users/sk/Google Drive/projects/cross-tissue")
source("src/base.R")

# tSNE plot config
tsnePlot = function(embed, ...) {
	plot(embed[,1], embed[,2],
		xlab="", ylab="",
		xaxt="n", yaxt="n",
		...
	)
}

# emat is genes x samples
# Filters transcripts based on SYNTAX correlation.
tsneSyntaxCmp = function(emat, pheno_matched, row_meta) {

	# Calculate correlations with SYNTAX score
	syntax_cor = cor(t(emat), pheno_matched$syntax_score, use="pairwise.complete")

	# Select highly correlated transcripts
	sel_transcripts = which(abs(syntax_cor) > 0.1)

	message("Selecting: ", length(sel_transcripts), " transcripts.")
	print(table(row_meta$tissue[sel_transcripts]))

	# Parallelized calculation of distance matrix
	message("Estimating distance matrix")
	dmat = Dist(t(emat[sel_transcripts, ]), method="euclidean", nbproc=6)

	# Run tSNE
	message("Running tSNE")
	embed = tsne(dmat, max_iter=20000, perplexity=20)

	return(embed)
}

# Plot tSNE panels colored by phenotype
tsnePlotSyntax = function(embed, pheno_matched) {
	# Plot results
	col = colorGradient(pheno_matched$syntax_score,
		gradlim=c(0, 100),
		colors=rev(brewer.pal(9, "Spectral")))
	tsnePlot(embed,
		col=col,
		pch=16,
		main="SYNTAX"
	)
	legendCol(colorRampPalette(rev(brewer.pal(9, "Spectral")))(20), c(0, 100))
	# text(embed[,1], embed[,2], labels=pheno_matched$syntax_score, cex=0.3)

	col = colorGradient(pheno_matched$DUKE,
		gradlim=c(0, 100),
		colors=rev(brewer.pal(9, "Spectral")))
	tsnePlot(embed,
		col=col,
		pch=16,
		main="DUKE"
	)
	# legendCol(rev(brewer.pal(9, "Spectral")), c(0, 100))
	legendCol(colorRampPalette(rev(brewer.pal(9, "Spectral")))(20), c(0, 100))
	# text(embed[,1], embed[,2], labels=pheno_matched$DUKE, cex=0.3)

	col = colorGradient(pheno_matched$ndv,
		# gradlim=range(pheno$ndv, na.rm=TRUE),
		gradlim=c(0, 3),
		colors=rev(brewer.pal(9, "Spectral")))
	tsnePlot(embed,
		col=col,
		pch=16,
		main="ndv"
	)
	legendCol(colorRampPalette(rev(brewer.pal(9, "Spectral")))(20), c(0, 3))
	# text(embed[,1], embed[,2], labels=pheno_matched$ndv, cex=0.3)

	col = colorGradient(pheno_matched$lesions,
		gradlim=c(0, 9),
		colors=rev(brewer.pal(9, "Spectral")))
	tsnePlot(embed,
		col=col,
		pch=16,
		main="lesions"
	)
	legendCol(colorRampPalette(rev(brewer.pal(9, "Spectral")))(20), c(0, 9))
}

barplotErr = function(numbers, barcol=rgb(0.93, 0.93, 0.93), pt_col=rgb(31, 120, 180, 200, maxColorValue=255), ...) {
	# Calculate stats
	df = data.frame(
		mean=sapply(numbers, mean, na.rm=TRUE),
		sd=sapply(numbers, sd, na.rm=TRUE),
		n=sapply(numbers, function(x) sum(!is.na(x)))
	)

	# bars = barplot(df$mean, ylim=c(0, max(df$mean + df$sd, na.rm=TRUE)), ...)
	bars = barplot(df$mean, col=barcol, ...)


	for (i in 1:length(numbers)) {
		points(
			jitter(
				rep(bars[i], length(numbers[[i]])),
				amount=0.2
			),
			numbers[[i]],
			xpd=TRUE,
			cex=0.4,
			pch=16,
			col=pt_col
		)
	}

	l_ply(seq_along(bars), function(i) {
		arrows(x0 = bars[i],
			# col="red",
			y0 = df$mean[i],
			x1 = bars[i], 
		    y1 = df$mean[i] + 1.96 * df$sd[i]/sqrt(df$n[i]),
		    code = 2,
		    length = 0.04, 
		    angle = 90,
		    lwd = 1.5)
	})
}


# Load data
# ---------------------------------------------------------------

# STARNET phenotype data
pheno = fread(file.path(
	"/Volumes/SANDY/phenotype_data",
	"STARNET_main_phenotype_table.cases.Feb_29_2016.tbl"
))

covar = fread(file.path(
	"/Volumes/SANDY/phenotype_data",
	"covariates.cases_and_controls.April_12_2016.txt"
))


# Senescence markers identified across 6 cell lines
senescent = read_excel("~/Google Drive/projects/senescence-markers/targets/from_zhidong/1368 genes significant across 6 cell lines.xlsx")

# Load imputed recast gene expression matrix
# load(file.path(data_dir, "STARNET/gene_exp_norm_batch_imp/all.RData"), verbose=TRUE)

# Load batch corrected expression data for each tissue
# load(file.path(data_dir, "STARNET/gene_exp_norm_batch/all.RData"), verbose=TRUE)

load(file.path(data_dir, "STARNET/gene_exp_norm_reshape/expr_recast.RData"), verbose=TRUE)
mat = expr_recast[, 3:ncol(expr_recast)]
mat = data.matrix(mat)
row_meta = expr_recast[, 1:2]


# Match phenotype data to selected gene expression matrix
pheno_matched = pheno[match(colnames(mat), pheno$starnet.ID), ]


# # Rename loaded normalized gene expression matrices
# names(expr_mats_batch) = sapply(
# 	strsplit(names(expr_mats_batch), "[.]"),
# 	function(x) x[4]
# )


# Rename columns to STARNET patient IDs
# expr_mats_norm = lapply(expr_mats_norm, function(mat) {
# 	colnames(mat) = sapply(
# 		strsplit(colnames(mat), "_"),
# 		function(x) x[2])
# 	return(mat)
# })

scaled_mat = t(scale(t(mat)))

scaled_mat[is.na(scaled_mat)] = 0

# # Principal components
# pca = prcomp(t(scaled_mat))

# embed = pca$x[, 1:2]

# par(mfrow=c(2, 2))
# tsnePlotSyntax(pca$x[, 1:2], pheno_matched)


# col = colorGradient(
# 	gradlim=range(as.numeric(pheno$starnet.ID), na.rm=T),
# 	# as.numeric(pheno_matched$starnet.ID),
# 	as.numeric(pheno_matched$starnet.ID),
# 	colors=rev(brewer.pal(9, "Spectral")))
# tsnePlot(pca$x[, 1:2],
# 	col=col,
# 	pch=16,
# 	main="STARNET ID"
# )
# legendCol(colorRampPalette(rev(brewer.pal(9, "Spectral")))(20), range(as.numeric(pheno$starnet.ID), na.rm=T))


# All transcripts combined across tissues
# -------------------------------------
dmat = Dist(t(scaled_mat), method="euclidean", nbproc=6)

embed = tsne(dmat, max_iter=2000, perplexity=30)

par(mfrow=c(3, 2))
tsnePlotSyntax(embed, pheno_matched)

col = colorGradient(
	gradlim=range(as.numeric(pheno$starnet.ID), na.rm=T),
	# as.numeric(pheno_matched$starnet.ID),
	as.numeric(pheno_matched$starnet.ID),
	colors=rev(brewer.pal(9, "Spectral")))
tsnePlot(embed,
	col=col,
	pch=16,
	main="STARNET ID"
)
legendCol(colorRampPalette(rev(brewer.pal(9, "Spectral")))(20), range(as.numeric(pheno$starnet.ID), na.rm=T))


# Tissue-specific tSNE plots
dmat_tissue = list()
embed_tissue = list()
# i = 1
for (i in 1:length(unique(row_meta$tissue))) {
	idx = row_meta$tissue == unique(row_meta$tissue)[i]

	dmat_tissue[[i]] = Dist(t(scaled_mat[idx, ]), method="euclidean", nbproc=6)
	embed_tissue[[i]] = tsne(dmat_tissue[[i]], max_iter=1000, perplexity=30)
}


names(dmat_tissue) = unique(row_meta$tissue)
names(embed_tissue) = unique(row_meta$tissue)


i = 6

names(embed_tissue)[i]
par(mfrow=c(2, 2))
tsnePlotSyntax(embed_tissue[[i]], pheno_matched)







# Senescent selection of transcripts
row_meta$gene_symbol = sapply(strsplit(as.character(row_meta$transcript_id), "_"), function(x) x[1])
sel_transcripts = row_meta$gene_symbol %in% senescent$Gene.Symbol



# ---------------------------------------

submat = t(scaled_mat[sel_transcripts & row_meta$tissue %in% c("AOR", "MAM"), ])  # this one!

# Randomize
# submat = submat[sample(nrow(submat)),]

clust = kmeans(submat, centers=20)$cluster

# Drop clusters with less than 3 members
clust[clust %in% which(table(clust) < 3)] = NA



dmat = Dist(submat, method="euclidean", nbproc=6)

# embed = tsne(dmat, max_iter=5000, perplexity=5)
embed = tsne(submat, max_iter=5000, perplexity=30)


pdf("pheno/plots/patient_clusters/aor_mam_senescence_tSNE.pdf", width=6)
par(mfrow=c(3, 2), mar=c(2, 2, 2, 3))
tsnePlotSyntax(embed, pheno_matched)

clust_cols = c(brewer.pal(9, "Set1"), brewer.pal(8, "Dark2"), brewer.pal(8, "Accent"))
col=clust_cols[as.numeric(factor(clust))]
tsnePlot(embed,
	col=col,
	pch=16,
	main="k-means clusters"
)

col = colorGradient(
	gradlim=range(as.numeric(pheno$starnet.ID), na.rm=T),
	# as.numeric(pheno_matched$starnet.ID),
	as.numeric(pheno_matched$starnet.ID),
	colors=rev(brewer.pal(9, "Spectral")))
tsnePlot(embed,
	col=col,
	pch=16,
	main="STARNET ID"
)
legendCol(colorRampPalette(rev(brewer.pal(9, "Spectral")))(20), range(as.numeric(pheno$starnet.ID), na.rm=T))
dev.off()


# Get phenotype  data for each cluster
syntax_scores = lapply(unique(clust), function(k) {
	pheno_matched$syntax_score[clust == k]
})

duke_scores = lapply(unique(clust), function(k) {
	pheno_matched$DUKE[clust == k]
})
lesions = lapply(unique(clust), function(k) {
	pheno_matched$lesions[clust == k]
})
ndv = lapply(unique(clust), function(k) {
	pheno_matched$ndv[clust == k]
})
LDL = lapply(unique(clust), function(k) {
	pheno_matched$LDL[clust == k]
})
BMI = lapply(unique(clust), function(k) {
	pheno_matched$BMI[clust == k]
})
age = lapply(unique(clust), function(k) {
	pheno_matched$Age[clust == k]
})


clust_order = order(
	sapply(syntax_scores, mean, na.rm=T) +
	sapply(duke_scores, mean, na.rm=T) +
	sapply(lesions, mean, na.rm=T) * 10,
	decreasing=TRUE)

syntax_scores = syntax_scores[clust_order]
duke_scores = duke_scores[clust_order]
lesions = lesions[clust_order]
ndv = ndv[clust_order]
LDL = LDL[clust_order]
BMI = BMI[clust_order]
age = age[clust_order]

pdf("pheno/plots/patient_clusters/aor_mam_senescence.pdf", width=5, height=6)
par(mfrow=c(4, 1), mar=c(2, 4, 2, 4))
light_red = rgb(251, 180, 174, 60, maxColorValue=255)
red = rgb(228, 26, 28, 150, maxColorValue=255)
barplotErr(syntax_scores, ylim=c(0, 100),
	barcol=light_red,
	pt_col=red,
	main="Clusters based on senescence markers in AOR and MAM",
	ylab="SYNTAX")
barplotErr(duke_scores, ylim=c(0, 100),
	barcol=light_red,
	pt_col=red,
	ylab="DUKE")
barplotErr(lesions, ylim=c(0, 9),
	barcol=light_red,
	pt_col=red,
	ylab="Lesions")
barplotErr(LDL, ylim=c(0, 9), ylab="LDL")
dev.off()

# barplotErr(age, ylim=c(0, 100), ylab="Age")
# barplotErr(BMI, ylim=range(pheno$BMI, na.rm=TRUE), ylab="BMI")

# barplotErr(ndv, ylim=c(0, 3), ylab="ndv")

# Test all pairwise combinations of clusters
pairwise_ttests = pairwise.t.test(pheno_matched$syntax_score, clust,
	p.adjust="BH",
	pool.sd=TRUE)

sort(pairwise_ttests$p.value)


pheno_tests = lapply(colnames(pheno_matched), function(name) {
	pairwise_ttests = pairwise.t.test(pheno_matched[[name]], clust,
		p.adjust="BH",
		pool.sd=TRUE)
	sort(pairwise_ttests$p.value)
})
names(pheno_tests) = colnames(pheno_matched)

min_pheno_pval = sapply(pheno_tests, min, na.rm=TRUE)

sig_min_pheno_pval = min_pheno_pval[which(min_pheno_pval < 0.05)]
sig_min_pheno_pval = sig_min_pheno_pval[order(sig_min_pheno_pval, decreasing=TRUE)]



pdf("pheno/plots/patient_clusters/aor_mam_senescence_pheno_tests.pdf", width=5, height=5)
dotchart(-log10(sig_min_pheno_pval),
	main="Cluster-phenotype distinctions",
	cex.main=1.0,
	pch=16,
	xlab=expression("max -log"[10] * "p (pairwise t-tests, BH)"))
dev.off()






syntax_clust = data.frame(
	mean=sapply(syntax_scores, mean, na.rm=TRUE),
	sd=sapply(syntax_scores, sd, na.rm=TRUE),
	n=sapply(syntax_scores, function(x) sum(!is.na(x)))
)
syntax_clust




# Dendrograms 
dend = as.dendrogram(hc)

dend = set(dend, "leaves_pch", 16)

leaves_col = colorGradient(pheno_matched$syntax_score,
			gradlim=c(0, 100),
			colors=rev(brewer.pal(9, "Spectral")))
leaves_col[is.na(leaves_col)] = "grey"

dend = hang.dendrogram(dend)

dend = set(dend, "leaves_col", leaves_col)
dend = set(dend, "labels_cex", 0.3)

plot(dend)

# Run tSNE
# message("Running tSNE")
# embed = tsne(dmat, max_iter=2000, perplexity=5)

# tsnePlotSyntax(embed, pheno_matched)

# # message("Estimating distance matrix")
# dmat = Dist(t(mat[sel_transcripts, ]), method="euclidean", nbproc=6)

# # Run tSNE
# message("Running tSNE")
# embed = tsne(dmat, max_iter=2000, perplexity=20)

# tsnePlotSyntax(embed, pheno_matched)


# Filtered correlation with SYNTAX plots
pdf("pheno/plots/all_tissue/syntax_cor_all_tSNE.pdf", width=7, height=6)
embed = tsneSyntaxCmp(mat, pheno_matched, row_meta)
tsnePlotSyntax(embed, pheno_matched)
dev.off()

pdf("pheno/plots/all_tissue/syntax_cor_all_tSNE_null.pdf", width=7, height=6)
embed = tsneSyntaxCmp(mat[, sample(ncol(mat))], pheno_matched, row_meta)
tsnePlotSyntax(embed, pheno_matched)
dev.off()


# tSNE plots for each tissue
expr_mats_batch = expr_mats_batch[names(expr_mats_batch) != "COR"]

# Run tSNE for each expression matrix
embed = lapply(expr_mats_batch, function(mat) {
	dmat = dist(t(mat))
	res = tsne(dmat, max_iter=2000, perplexity=15)
	return(res)
})


# Selected phenotype for tSNE plots
# phenotype = "syntax_score"
phenotype = "BMI"
# phenotype = "ID"
# phenotype = "LDL"

# expr_mats_batch

pdf(paste0("pheno/plots/tSNE_", phenotype, ".pdf"), width=16)
par(mfrow=c(2, 5))
for (i in 1:length(expr_mats_batch)) {
	print(names(expr_mats_batch)[i])

	# i = 2
	mat = expr_mats_batch[[i]]
	res = embed[[i]]

	# Match phenotype data to selected gene expression matrix
	patient_ids = sapply(
		strsplit(colnames(mat), "_"),
		function(x) x[2]
	)

	pheno_matched = pheno[match(patient_ids, pheno$starnet.ID), ]
	covar_matched = covar[match(colnames(mat), covar$sample), ]

	# Pick color scheme
	if (phenotype == "syntax_score") {
		# # SYNTAX score
		# col = colorGradient(pheno_matched$syntax_score,
		# 	gradlim=c(0, 100),
		# 	# range(pheno$syntax_score, na.rm=T),
		# 	colors=rev(brewer.pal(9, "Spectral")))
		col = colorGradient(pheno_matched$syntax_score,
			gradlim=c(0, 100),
			# range(pheno$syntax_score, na.rm=T),
			colors=brewer.pal(9, "YlGnBu"))
	} else if (phenotype == "LDL") {
		# LDL 
		col = colorGradient(as.numeric(pheno_matched$LDL),
			gradlim=range(pheno$LDL),
			colors=rev(brewer.pal(9, "Spectral")))
	} else if (phenotype == "BMI") {
		# BMI
		col = colorGradient(as.numeric(pheno_matched$BMI),
			gradlim=range(pheno$BMI, na.rm=T),
			colors=rev(brewer.pal(9, "Spectral")))
	} else if (phenotype == "ID") {
		col = colorGradient(
			gradlim=range(as.numeric(pheno$starnet.ID), na.rm=T),
			# as.numeric(pheno_matched$starnet.ID),
			as.numeric(covar_matched$subject),
			colors=rev(brewer.pal(9, "Spectral")))
	} else {
		stop("Invalid phenotype")
	}

	# col=colorGradient(pheno_matched$DUKE, colors=rev(brewer.pal(9, "Spectral")))
	# col=colorGradient(pheno_matched$lesions, colors=rev(brewer.pal(9, "Spectral")))
	# col=colorGradient(as.numeric(pheno_matched$Smoking.Years), colors=rev(brewer.pal(9, "Spectral")))
	# col=colorGradient(as.numeric(covar_matched$read_length), colors=brewer.pal(9, "Spectral"))
	# col=colorGradient(as.numeric(pheno_matched$Age), colors=brewer.pal(9, "Spectral"))
	plot(res[,1], res[,2],
		col=col,
		pch=16,
		cex=1.5,
		main=names(expr_mats_batch)[i],
		xlab="", ylab=""
	)
	# text(res[,1], res[,2], labels=pheno_matched$starnet.ID, cex=0.5)
	# text(res[,1], res[,2], labels=pheno_matched$syntax_score, cex=0.5)
	# text(res[,1], res[,2], labels=pheno_matched$DUKE, cex=0.5)
}

if (phenotype == "syntax_score") {
	plotColorBar(
		colorGradient(seq(0, 100, length.out=100),
			colors=brewer.pal(9, "YlGnBu")),
		min=0, max=100, title="SYNTAX")
} else if (phenotype == "ID") {
	plotColorBar(
		colorGradient(
			seq(
				min(as.numeric(pheno$starnet.ID), na.rm=T),
				max(as.numeric(pheno$starnet.ID), na.rm=T),
				length.out=100),
		colors=rev(brewer.pal(9, "Spectral"))),
		min=0,
		max=max(as.numeric(pheno$starnet.ID), na.rm=T),
		# min=0,
		# max=100,
		title="STARNET ID")
}
dev.off()


# Plot color bars for selected tSNE plots
pdf("pheno/plots/legends.pdf", width=4.0)
par(mfrow=c(1, 2))
plotColorBar(
	colorGradient(seq(0, 100, length.out=100),
		colors=brewer.pal(9, "YlGnBu")),
	min=0, max=100, title="SYNTAX")

plotColorBar(
	colorGradient(
		seq(
			min(as.numeric(pheno$starnet.ID), na.rm=T),
			max(as.numeric(pheno$starnet.ID), na.rm=T),
			length.out=100),
	colors=rev(brewer.pal(9, "Spectral"))),
	min=0,
	max=max(as.numeric(pheno$starnet.ID), na.rm=T),
	# min=0,
	# max=100,
	title="STARNET ID")
dev.off()
