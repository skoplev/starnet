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

# endocrine = fread("~/Google Drive/projects/STARNET-endocrine/data/endo_scores_large.csv")

endocrine = fread("~/Google Drive/projects/STARNET-endocrine/data/endo_scores.csv")

# Main module table
mod_tab = fread("co-expression/tables/module_tab.csv")


# Endocrine target enrichment
endocrine_enrich = fread("co-expression/tables/endocrine_tab.csv")


endocrine_enrich$endocrine_ids = paste(endocrine_enrich$endocrine_factor, endocrine_enrich$from_tissue, endocrine_enrich$target_tissue, sep="_")

# Pick endocrine factor with highest coverage (first found)
# endocrine_enrich = endocrine_enrich[match(unique(endocrine_enrich$endocrine_ids), endocrine_enrich$endocrine_ids), ]

# endocrine_enrich = endocrine_enrich[endocrine_enrich$module_coverage > 0.1, ]

hist(endocrine_enrich$module_coverage, breaks=50)

# sum(endocrine_enrich$padj < 0.1)

d_sec_score = list(
	"-"=endocrine_enrich$sec_score[!endocrine_enrich$endocrine_in_module],
	"+"=endocrine_enrich$sec_score[endocrine_enrich$endocrine_in_module]
)

d_module_coverage = list(
	"-"=endocrine_enrich$module_coverage[!endocrine_enrich$endocrine_in_module],
	"+"=endocrine_enrich$module_coverage[endocrine_enrich$endocrine_in_module]
)


# sapply(d, length)
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


sub_endocrine_enrich = endocrine_enrich[endocrine_enrich$endocrine_in_module, ]

sub_endocrine_enrich[order(sub_endocrine_enrichr$sec_score, decreasing=TRUE), ]

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

n = 20

barplot(sub_endocrine_enrich$sec_score[1:n],
	names.arg=sub_endocrine_enrich$endocrine_factor[1:n],
	las=2
	)




# plot(density(d[[1]]))
# lines(density(d[[2]]), col="red")

unique(endocrine_enrich$endocrine_factor[endocrine_enrich$endocrine_in_module])
length(unique(endocrine_enrich$endocrine_factor[endocrine_enrich$endocrine_in_module]))

sub_endocrine_enrich[sub_endocrine_enrich$target_tissue == "AOR", ]
sub_endocrine_enrich[sub_endocrine_enrich$target_tissue == "MAM", ]
sub_endocrine_enrich[sub_endocrine_enrich$target_tissue == "BLOOD", ]

# Load cross-tissue modules
# ------------------------------
between = new.env()
load(file.path(data_dir, "modules/between_within-cross-tissue.RData"),
	between,
	verbose=TRUE)

# Parse module data
between = parseModuleData(between)


# Load module annotations
module_tab = fread("co-expression/tables/module_tab.csv")


# Find modules of each endocrine factor
# ----------------------------------
endocrine_ids = paste(endocrine$from_tissue, endocrine$endocrine_factor, sep=":")
transcript_ids = paste(between$meta_genes$tissue, between$meta_genes$gene_symbol, sep=":")

endocrine$clust = between$clust[match(endocrine_ids, transcript_ids)]

# Add module purity to endocrine table
endocrine$mod_purity = module_tab$purity[endocrine$clust]

# Cross-tissue?
endocrine$cross_tissue_10 = endocrine$mod_purity < 0.9
endocrine$cross_tissue_5 = endocrine$mod_purity < 0.95
endocrine$cross_tissue_1 = endocrine$mod_purity < 0.99

# Write example table for particular module
write.csv(endocrine[endocrine$clust == 98, ], "co-expression/plots/endocrines_mod98.csv")

sec_score_thresh = 2

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

df = pairwise_stats[[1]]

par(mfrow=c(3, ceiling(length(pairwise_stats)/3)))

for (i in 1:length(pairwise_stats)) {
	df = pairwise_stats[[i]]

	plot(df$cutoff, df$frac * 100,
		type="l",
		lwd=2,
		ylim=c(0, 100),
		main=paste0(t1, "-", t2),
		ylab="Cross-tissue %",
		xlab="> |r|")

}


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



sel = mod_tab$purity < 0.99

# d = list(
# 	"-"=-log10(mod_tab$secreted_protein_qvalue)[sel],
# 	"+"=-log10(mod_tab$secreted_protein_qvalue)[!sel]
# )

d = list(
	"-"=-log10(mod_tab$CAD_qvalue)[!sel],
	"+"=-log10(mod_tab$CAD_qvalue)[sel]
)


# d = list(
# 	"-"=-log10(mod_tab$CAD_pval)[sel],
# 	"+"=-log10(mod_tab$CAD_pval)[!sel]
# )


# plot(density(d[[1]]))
# lines(density(d[[2]]), col="red")

d = lapply(d, function(x) {
	x[x > 10] = 10
	return(x)
})

sapply(d, length)

boxplot(d)
wilcox.test(d[[1]], d[[2]])



sel = list(
	"Secreted"=mod_tab$secreted_protein_qvalue < 0.1,
	"CAD"=mod_tab$CAD_qvalue < 0.1
)


