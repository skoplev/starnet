# First approach for analysing cross-tissue modules and endocrine factors
#
rm(list=ls())

library(WGCNA)
enableWGCNAThreads(nThreads=8)

library(data.table)
library(RColorBrewer)
library(gplots)
library(igraph)

library(compiler)
enableJIT(3)

data_dir = "/Users/sk/DataProjects/cross-tissue"  # root of data directory
setwd("/Users/sk/Google Drive/projects/cross-tissue")
source("src/parse.R")


# Load cross-tissue modules
# ------------------------------
between = new.env()
load(file.path(data_dir, "modules/between_within-cross-tissue.RData"),
	between,
	verbose=TRUE)

# Parse module data
between = parseModuleData(between)

# Main module table
mod_tab = fread("co-expression/tables/module_tab.csv")

# Endocrine factors
endocrine = fread("~/Google Drive/projects/STARNET-endocrine/data/endo_scores.csv")

# Load endocrine target enrichment
endocrine_enrich = fread("co-expression/tables/endocrine_tab.csv")
endocrine_enrich$endocrine_ids = paste(endocrine_enrich$endocrine_factor, endocrine_enrich$from_tissue, endocrine_enrich$target_tissue, sep="_")

# Pick endocrine factor with highest coverage (first found)
# endocrine_enrich = endocrine_enrich[match(unique(endocrine_enrich$endocrine_ids), endocrine_enrich$endocrine_ids), ]

# endocrine_enrich = endocrine_enrich[endocrine_enrich$module_coverage > 0.1, ]

# hist(endocrine_enrich$module_coverage, breaks=50)

# sum(endocrine_enrich$padj < 0.1)

# Plots of endocrine target enrichment stratified by cross-tissue/tissue-specific modules
# --------------------------------------------------------------------------
d_sec_score = list(
	"-"=endocrine_enrich$sec_score[!endocrine_enrich$endocrine_in_module],
	"+"=endocrine_enrich$sec_score[endocrine_enrich$endocrine_in_module]
)

d_module_coverage = list(
	"-"=endocrine_enrich$module_coverage[!endocrine_enrich$endocrine_in_module],
	"+"=endocrine_enrich$module_coverage[endocrine_enrich$endocrine_in_module]
)

names(d_sec_score) = paste(names(d_sec_score), sapply(d_sec_score, length))
names(d_module_coverage) = paste(names(d_module_coverage), sapply(d_module_coverage, length))

pdf("co-expression/plots/endocrine_module_enrichment.pdf", width=5, height=5)
par(mar=c(10.1, 4.1, 4.1, 2.1), mfrow=c(1, 2))
boxplot(d_sec_score, outline=FALSE,
	las=2,
	main=paste0("p=", format(wilcox.test(d_sec_score[[1]], d_sec_score[[2]])$p.value, digits=3)),
	col=rgb(0.95, 0.95, 0.95),
	ylab="Secretion score")

boxplot(d_module_coverage, outline=FALSE,
	las=2,
	main=paste0("p=", format(wilcox.test(d_module_coverage[[1]], d_module_coverage[[2]])$p.value, digits=3)),
	col=rgb(0.95, 0.95, 0.95),
	ylab="Module coverage")
dev.off()


# Endocrines found in cross-tissue modules
sub_endocrine_enrich = endocrine_enrich[endocrine_enrich$endocrine_in_module, ]
sub_endocrine_enrich[order(sub_endocrine_enrich$sec_score, decreasing=TRUE), ]

# Order by secretion score
sub_endocrine_enrich = sub_endocrine_enrich[order(sub_endocrine_enrich$sec_score, decreasing=TRUE), ]

# Count tissue combinations

tissues = levels(factor(sub_endocrine_enrich$from_tissue))

endocrine_counts = matrix(NA, ncol=length(tissues), nrow=length(tissues))
colnames(endocrine_counts ) = tissues
rownames(endocrine_counts ) = tissues

for (i in 1:length(tissues)) {
	for (j in 1:length(tissues)) {
		endocrine_counts[i, j] = sum(
			sub_endocrine_enrich$from_tissue == tissues[i] &
			sub_endocrine_enrich$target_tissue == tissues[j]
		)
	}
}

pdf("co-expression/plots/endocrine/endocrine_cross-tissue_heatmap.pdf", height=6, width=6)
heatmap.2(endocrine_counts, trace="none",
	col=colorRampPalette(brewer.pal(9, "YlOrRd"))(100),
	cellnote=endocrine_counts,
	notecol="black",
	mar=c(8, 8),
	ylab="From", xlab="To",
	key.title="",
	key.xlab="Endocrines",
	key.ylab=""
)
dev.off()

# adj_mat = endocrine_counts
adj_mat = log10(endocrine_counts)
adj_mat[is.infinite(adj_mat)] = 0

# adj_mat[adj_mat <= 1] = 0
g = graph_from_adjacency_matrix(adj_mat, weighted="endocrines")
write_graph(g, "co-expression/plots/endocrine_tissue_graph.gml", format="gml")

# plot(g)

# head(, 100)

unique(sub_endocrine_enrich$endocrine_ids)

# Top-20 endocrine scores
colors = brewer.pal(9, "Set1")[-6]
tissues = c("AOR", "BLOOD", "LIV", "MAM", "SKLM", "SF", "VAF")
# levels(factor(endocrine_enrich$from_tissue))  # vector of tissues 


pdf("co-expression/plots/endocrine/endocrine_cross_tissue_ranking.pdf", width=6)
idx = 1:20


pdf("co-expression/plots/endocrine/endocrine_cross_tissue_ranking_MAM_not_AOR.pdf", width=8)
idx = which(sub_endocrine_enrich$from_tissue != "AOR" & sub_endocrine_enrich$target_tissue == "MAM")


# idx = which(sub_endocrine_enrich$from_tissue == "VAF" & sub_endocrine_enrich$target_tissue == "MAM")

par(mfrow=c(2, 1))
bar_col = rgb(0.4, 0.4, 0.4)
barplot(sub_endocrine_enrich$sec_score[idx],
	names.arg=paste0(sub_endocrine_enrich$endocrine_factor[idx], " (", sub_endocrine_enrich$module[idx], ")"),
	ylab="Endocrine effect",
	las=2,
	col=bar_col,
	border=bar_col
	)

sub_endocrine_enrich[idx, ]

tissue_mat = sub_endocrine_enrich[idx, c("from_tissue", "target_tissue")]

tissue_mat$from_tissue = factor(tissue_mat$from_tissue, levels=tissues)
tissue_mat$target_tissue = factor(tissue_mat$target_tissue, levels=tissues)

plot(1:length(idx), rep(1, length(idx)), col=colors[as.integer(tissue_mat$from_tissue)],
	xlab="", ylab="",
	bty="n",
	cex=2.0,
	pch=15)
points(1:length(idx), rep(1.2, length(idx)), col=colors[as.integer(tissue_mat$target_tissue)],
	cex=2.0,
	pch=15)
dev.off()



# plot(density(d[[1]]))
# lines(density(d[[2]]), col="red")

unique(endocrine_enrich$endocrine_factor[endocrine_enrich$endocrine_in_module])
length(unique(endocrine_enrich$endocrine_factor[endocrine_enrich$endocrine_in_module]))

sub_endocrine_enrich[sub_endocrine_enrich$target_tissue == "AOR", ]
sub_endocrine_enrich[sub_endocrine_enrich$target_tissue == "MAM", ]
sub_endocrine_enrich[sub_endocrine_enrich$target_tissue == "BLOOD", ]


# Find modules of each endocrine factor
# ----------------------------------
endocrine_ids = paste(endocrine$from_tissue, endocrine$endocrine_factor, sep=":")
transcript_ids = paste(between$meta_genes$tissue, between$meta_genes$gene_symbol, sep=":")

endocrine$clust = between$clust[match(endocrine_ids, transcript_ids)]

# Add module purity to endocrine table
endocrine$mod_purity = mod_tab$purity[endocrine$clust]

# Cross-tissue?
endocrine$cross_tissue_10 = endocrine$mod_purity < 0.9
endocrine$cross_tissue_5 = endocrine$mod_purity < 0.95
endocrine$cross_tissue_1 = endocrine$mod_purity < 0.99

# Write example table for particular module
write.csv(endocrine[endocrine$clust == 98, ], "co-expression/plots/endocrines_mod98.csv")

sec_score_thresh = 2

sum(endocrine$sec_score > sec_score_thresh)
unique(endocrine$endocrine_factor[endocrine$sec_score > sec_score_thresh])
# unique(endocrine$from_tissue[endocrine$sec_score > sec_score_thresh])
# unique(endocrine$target_tissue[endocrine$sec_score > sec_score_thresh])


pdf("co-expression/plots/cross_tissue_endocrine.pdf", height=3.5, width=7)
par(mfcol=c(1, 2),
	mar=c(5.1, 4.1, 4.1, 2.1)  # default
)
hist(endocrine$sec_score,
	breaks=100,
	xlab="Endocrine secretion score",
	main="",
	col="black")
abline(v=2, col=brewer.pal(9, "Set1")[1], lty=3)
# sum(endocrine$sec_score > sec_score_thresh)

unique(endocrine$endocrine_factor)
unique(endocrine$endocrine_factor[endocrine$sec_score > sec_score_thresh])

mat = cbind(
	table(endocrine$cross_tissue_1[endocrine$sec_score > sec_score_thresh]),
	table(endocrine$cross_tissue_5[endocrine$sec_score > sec_score_thresh]),
	table(endocrine$cross_tissue_10[endocrine$sec_score > sec_score_thresh])
)
colnames(mat) = c("1", "5", "10")
mat = mat[2:1, ]  # flip rows


# Calculate cross-tissue percentages
mat[1,] / apply(mat, 2, sum)

par(mar=c(5.1, 4.1, 4.1, 8.1))
barplot(
	mat,
	ylab="Endocrines",
	xlab="Cross-tissue transcripts (%)"
)
dev.off()


# Load expression data
# -----------------------------------
emat_file = "STARNET/gene_exp_norm_reshape/expr_recast.RData"

load(file.path(data_dir, emat_file), verbose=TRUE)

emat = expr_recast[, 3:ncol(expr_recast)]
meta_genes = expr_recast[, 1:2]
meta_genes = as.data.frame(meta_genes)

rm(expr_recast)

# Test if the gene metadata is the same as for the cross-tissue modules
if (!all(meta_genes$transcript_id == between$meta_genes$transcript_id)) {
	stop("Transcript mismatch")
}

tissues = unique(meta_genes$tissue)

# t1 = "VAF"
# t2 = "LIV"

pairwise_stats = list()
for (i in 1:length(tissues)) {
	for (j in 1:length(tissues)) {
		if (j >= i) next
		gc()  # garbage collect

		message(tissues[i], " ", tissues[j])
		t1 = tissues[i]
		t2 = tissues[j]

		idx1 = which(meta_genes$tissue == t1)
		idx2 = which(meta_genes$tissue == t2)

		# Calculate cross-tissue correlations
		cmat = cor(
			t(emat[idx1, ]),
			t(emat[idx2, ]),
			use="pairwise.complete.obs")

		# Adjacency matrix
		cmat = abs(cmat)

		# cutoffs = seq(0, 1, length.out=11)
		cutoffs = seq(0, 1, length.out=31)  # adjacency cutoffs to test

		edges_stats = sapply(cutoffs, function(adj_cutoff) {
			message(adj_cutoff)

			# Get local edges with high correlations
			local_index = which(cmat > adj_cutoff, arr.ind=TRUE)

			# Translate edge indices to global module index
			edges = data.frame(
				m1=idx1[local_index[, 1]],
				m2=idx2[local_index[, 2]]
			)

			edges$clust1 = between$clust[edges$m1]
			edges$clust2 = between$clust[edges$m2]

			# Are the transcript pairs found in same module?
			edges$same_module = edges$clust1 == edges$clust2

			# n_edges = sapply(cross_tissue_edges, function(edges) sum(edges$same_module))

			# Collect stats
			stats = data.frame(
				cutoff=adj_cutoff,
				frac=sum(edges$same_module) / nrow(edges),
				edges=sum(edges$same_module))

			return(stats)
		})

		# Store
		edges_stats = as.data.frame(t(edges_stats))
		print(edges_stats)
		pairwise_stats[[paste(t1, t2, sep="_")]] = edges_stats
	}
}

save(pairwise_stats, file=file.path(data_dir, "R_workspaces/endocrineCrossTissueCoverage_pairwise_stats.RData"))



# Get all fractions, as matrix
all_frac = sapply(pairwise_stats, function(stats) {
	unlist(stats$frac)
})


# Calculate mean fraction
mean_frac = apply(all_frac, 1, mean, na.rm=TRUE)

# Calculate 95% confidence intervals
conf_int = apply(all_frac, 1, function(x) {
	if (all(is.na(x))) {
		return(c(NA, NA))
	} else {
		return(t.test(x)$conf.int)
	}
})
conf_int = t(conf_int)

conf_int[conf_int < 0] = 0

conf_int[31, 1] = 0
conf_int[31, 2] = 1



# Plot of mean cross-tissue coverage
pdf("co-expression/plots/cross-tissue-correlations.pdf", width=6, height=5)
mid_r = 0.5
tissue_col = brewer.pal(9, "Set1")[-6]
conf_col = brewer.pal(9, "Pastel1")[2]

plot(cutoffs, mean_frac, ylim=c(0, 1), type="n",
	xlab="Edge criteria (> |r|)",
	ylab="Cross-tissue interactions"
)

# 95% confidence interval area
polygon(c(cutoffs, rev(cutoffs)), c(conf_int[, 1], rev(conf_int[, 2])),
	col=conf_col,
	border=NA)


lines(cutoffs, mean_frac)

mid = which(cutoffs == mid_r)
mid_points = sapply(pairwise_stats, function(stats) stats$frac[mid])

t1 = sapply(strsplit(names(mid_points), "_"), function(x) x[1])
t2 = sapply(strsplit(names(mid_points), "_"), function(x) x[2])

abline(v=mid_r, col="grey", lty=3)

x_jit = jitter(rep(mid_r, length(pairwise_stats)), 3)
points(x_jit, mid_points,
	pch=16,
	col=tissue_col[match(t1, tissues)])

points(x_jit, mid_points,
	# pch=16,
	col=tissue_col[match(t2, tissues)])

text(x_jit, mid_points,
	labels=names(mid_points),
	pos=4,
	cex=0.3)

legend("topleft", legend=c("95% CI", "Mean", "Tissue A", "Tissue B"), pch=c(15, 3, 16, 1), col=c(conf_col, "black", "black", "black"))
dev.off()



pdf("co-expression/plots/cross-tissue-edges-in-modules.pdf", width=5, height=8)
par(mfrow=c(2, 1))
plot(df$cutoff, df$frac * 100,
	type="l",
	lwd=2,
	ylim=c(0, 100),
	main=paste0(t1, "-", t2),
	ylab="Cross-tissue %",
	xlab="> |r|")

plot(df$cutoff, log10(df$edges),
	main=paste0(t1, "-", t2),
	type="l",
	lwd=2,
	xlab="> |r|")
dev.off()


# i = 142917
# j = 63705

# between$meta_genes[i, ]
# between$meta_genes[j, ]

# between$clust[i]
# between$clust[j]


# Scatter plot comparing the module enrichment of CAD GWAS and secreted proteins
pdf("co-expression/plots/endocrine/CAD_secreted_enrichment.pdf", width=4, height=4.2)
purity_thresh = 0.99  # cross-tissue definition

sum(mod_tab$purity < purity_thresh) / length(mod_tab$purity)

x = -log10(mod_tab$CAD_qvalue)
y = -log10(mod_tab$secreted_protein_qvalue)

x[is.infinite(x)] = 16
y[is.infinite(y)] = 16

mod_colors = rep("black", length(x))
mod_colors[mod_tab$purity < purity_thresh] = brewer.pal(9, "Set1")[1]

sum(mod_tab$CAD_qvalue < 0.1 & mod_tab$secreted_protein_qvalue < 0.1)

n_cross_tissue = sum(mod_tab$CAD_qvalue < 0.1 &
	mod_tab$secreted_protein_qvalue < 0.1 &
	mod_tab$purity < purity_thresh)

n_tissue_specific = sum(mod_tab$CAD_qvalue < 0.1 &
	mod_tab$secreted_protein_qvalue < 0.1 &
	mod_tab$purity > purity_thresh)

plot(x, y,
	xlab="CAD enrichment (-log10 q)",
	ylab="Secreted enrichment (-log10 q)",
	main=paste0("p=", format(cor.test(x, y)$p.value, digits=3), " CT: ", n_cross_tissue, " TS: ", n_tissue_specific),
	pch=16,
	col=mod_colors
)

abline(v=-log10(0.1), col="grey", lty=2)
abline(h=-log10(0.1), col="grey", lty=2)

legend("topright",
	legend=c("Tissue-specific", "Cross-tissue"),
	pch=16,
	col=c("black", brewer.pal(9, "Set1")[1])
)

dev.off()



# sel = mod_tab$purity < 0.99

# # d = list(
# # 	"-"=-log10(mod_tab$secreted_protein_qvalue)[sel],
# # 	"+"=-log10(mod_tab$secreted_protein_qvalue)[!sel]
# # )

# d = list(
# 	"-"=-log10(mod_tab$CAD_qvalue)[!sel],
# 	"+"=-log10(mod_tab$CAD_qvalue)[sel]
# )


# # d = list(
# # 	"-"=-log10(mod_tab$CAD_pval)[sel],
# # 	"+"=-log10(mod_tab$CAD_pval)[!sel]
# # )


# # plot(density(d[[1]]))
# # lines(density(d[[2]]), col="red")

# d = lapply(d, function(x) {
# 	x[x > 10] = 10
# 	return(x)
# })

# sapply(d, length)

# boxplot(d)
# wilcox.test(d[[1]], d[[2]])



# sel = list(
# 	"Secreted"=mod_tab$secreted_protein_qvalue < 0.1,
# 	"CAD"=mod_tab$CAD_qvalue < 0.1
# )


