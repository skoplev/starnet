rm(list=ls())

library(RColorBrewer)
library(gplots)

library(devtools)  # for installing heatmap.3
# Load heatmap.3
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")

library(compiler)
enableJIT(3)

data_dir = "/Users/sk/DataProjects/cross-tissue"  # root of data directory
setwd("/Users/sk/Google Drive/projects/cross-tissue")

countModuleTissueStat = function(modules) {

	out = list()

	out$tissue_counts = sapply(1:max(modules$clust), function(cl) {
		table(modules$meta_genes$tissue[modules$clust == cl])
	})

	out$purity = apply(out$tissue_counts, 2, max) / apply(out$tissue_counts, 2, sum)

	out$n_tissues = apply(out$tissue_counts, 2, function(col) {
		sum(col > 0)
	})

	out$size = apply(out$tissue_counts, 2, sum)

	return(out)
}

between_within = new.env()
load(file.path(data_dir, "modules/between_within-cross-tissue.RData"),
	between_within,
	verbose=TRUE)


complete = new.env()
load(file.path(data_dir, "modules/complete-cross-tissue.RData"),
	complete,
	verbose=TRUE)

sing = new.env()
load(file.path(data_dir, "modules/single-cross-tissue.RData"),
	sing,
	verbose=TRUE)

between_within$clust = as.numeric(factor(between_within$bwnet$colors))
complete$clust = as.numeric(factor(complete$bwnet$colors))
sing$clust = as.numeric(factor(sing$bwnet$colors))

sing_stats = countModuleTissueStat(sing)
comp_stats = countModuleTissueStat(complete)
bw_stats = countModuleTissueStat(between_within)

mod_stats = list("Between-within"=bw_stats,
	"Complete"=comp_stats,
	"Single"=sing_stats)

cols = brewer.pal(9, "Set1")[c(1:5, 7:9)]


pdf("co-expression/plotsModuleStats/module_size.pdf", width=4, height=5)
par(mfrow=c(3, 2))
barplot(
	log10(rev(sort(table(between_within$clust)))),
	col=cols[1], border=cols[1],
	main=paste0("Within-between, n=", length(unique(between_within$clust))),
	ylim=c(1, 4),
	xpd=FALSE,
	names.arg="",
	ylab=expression("log"[10] * " size"),
	xlab="Modules")

barplot(table(bw_stats$n_tissues),
	xlab="Tissues",
	ylab="Modules",
	main="Within-between",
	space=0,
	col=cols[1], border=cols[1])

barplot(
	log10(rev(sort(table(complete$clust)))),
	col=cols[2], border=cols[2],
	main=paste0("Complete, n=", length(unique(complete$clust))),
	ylim=c(1, 4),
	xpd=FALSE,
	names.arg="",
	ylab=expression("log"[10] * " size"),
	xlab="Modules")

barplot(table(comp_stats$n_tissues),
	xlab="Tissues",
	ylab="Modules",
	main="Complete",
	space=0,
	col=cols[2],
	border=cols[2])


barplot(
	log10(rev(sort(table(sing$clust)))),
	col=cols[3], border=cols[3],
	main=paste0("Single, n=", length(unique(sing$clust))),
	ylim=c(1, 4),
	xpd=FALSE,
	names.arg="",
	ylab=expression("log"[10] * " size"),
	xlab="Modules")

barplot(table(sing_stats$n_tissues),
	xlab="Tissues",
	ylab="Modules",
	main="Single",
	space=0,
	col=cols[3], border=cols[3])
dev.off()


# Plots counting number of cross tissue modules
# ------------------------------------------------
pdf("co-expression/plotsModuleStats/cross_tissue_module_counts.pdf", width=3, height=6)
cross_tissue_mod_05 = sapply(mod_stats, function(stat) {
	sum(stat$purity < 0.95)
})

cross_tissue_mod_005 = sapply(mod_stats, function(stat) {
	sum(stat$purity < 0.995)
})

par(mfrow=c(2, 1))
barplot(cross_tissue_mod_05, las=2, col=cols, main="Cross-tissue >5%", ylab="Modules")
barplot(cross_tissue_mod_005, las=2, col=cols, main="Cross-tissue >0.5%", ylab="Modules")
dev.off()


# Module counts based on purity
# ------------------------------------
pdf("co-expression/plotsModuleStats/cross_tissue_module_purity_size.pdf", width=4, height=4.5)
size_range = range(bw_stats$size, comp_stats$size, sing_stats$size)
size_range = c(0, 2000)
other_range = range(1 - bw_stats$purity, 1 - comp_stats$purity, 1 - sing_stats$purity)

plot(bw_stats$size, 1 - bw_stats$purity, col=cols[1],
	pch=16,
	xlab="Module size",
	ylab="Cross-tissue fraction",
	xlim=size_range, ylim=other_range)
points(comp_stats$size, 1 - comp_stats$purity,
	pch=16,
	col=cols[2])
points(sing_stats$size, 1 - sing_stats$purity,
	pch=16,
	col=cols[3])
abline(h=0.05, col="grey", lty=2, lwd=1.5)
legend("topright", legend=c("Between-within", "Complete", "Single"), pch=16, col=cols)
dev.off()


# Plot jaccard similarity of detected modules

# modules1 = between_within
# modules2 = complete

jaccard = function(a, b) length(intersect(a, b)) / length(union(a, b))

moduleSimil = function(modules1, modules2) {

	modules1$meta_genes$id = paste(modules1$meta_genes$tissue,
		modules1$meta_genes$transcript_id,
		sep="_")
	modules2$meta_genes$id = paste(modules2$meta_genes$tissue,
		modules2$meta_genes$transcript_id,
		sep="_")

	# if (!all(sapply(modules1$clust, is.integer))) {
	# 	stop("Cluster encoding is are not integers")
	# }

	# if (!all(sapply(modules2$clust, is.integer))) {
	# 	stop("Cluster encoding is are not integers")
	# }

	simil = matrix(NA,
		nrow=max(modules1$clust),
		ncol=max(modules2$clust))

	rownames(simil) = 1:max(modules1$clust)
	colnames(simil) = 1:max(modules2$clust)

	for (i in 1:max(modules1$clust)) {
		idx1 = modules1$clust == i
		gene_tissue_ids1 = modules1$meta_genes$id[idx1]

		for (j in 1:max(modules2$clust)) {
			idx2 = modules2$clust == j
			gene_tissue_ids2 = modules2$meta_genes$id[idx2]
			simil[i, j] = jaccard(gene_tissue_ids1, gene_tissue_ids2)
		}
	}
	return(simil)
}


simil_bw_comp = moduleSimil(between_within, complete)
simil_bw_sing = moduleSimil(between_within, sing)
simil_sing_comp = moduleSimil(sing, complete)


pdf("co-expression/plotsModuleStats/best_module_match.pdf", width=5, height=6)
par(mfrow=c(3, 2))
plot(density(apply(simil_bw_comp, 1, max)),
	lwd=2.0,
	main="within-between -> complete",
	col=cols[1],
	xlim=c(0, 1))  # within-between -> complete
plot(density(apply(simil_bw_comp, 2, max)),
	lwd=2.0,
	main="complete -> within-between",
	col=cols[2],
	xlim=c(0, 1))  # complete -> within-between

plot(density(apply(simil_bw_sing, 1, max)),
	lwd=2.0,
	main="within-between -> single",
	col=cols[1],
	xlim=c(0, 1))  # within-between -> single
plot(density(apply(simil_bw_sing, 2, max)),
	lwd=2.0,
	main="single -> within-between",
	col=cols[3],
	xlim=c(0, 1))  # single -> within-between

plot(density(apply(simil_sing_comp, 1, max)),
	lwd=2.0,
	main="single -> complete",
	col=cols[3],
	xlim=c(0, 1))  # single -> complete 
plot(density(apply(simil_sing_comp, 2, max)),
	lwd=2.0,
	main="complete -> single",
	col=cols[2],
	xlim=c(0, 1))  # complete -> single
dev.off()


# i = 4
# j = 9

# sort(unique(complete$clust))

# between_within$meta_genes[between_within$clust == i, ]
# complete$meta_genes[complete$clust == j, ]

# comp_stats$tissue_counts

colorByPurity = function(tissue_counts, cols) {
	module_cols = apply(tissue_counts, 2, function(counts) {
		tissue_present = counts / sum(counts) > 0.05  # purity threshold
		module_colors = rep("white", length(tissue_present))  # default
		module_colors[tissue_present] = cols[tissue_present]
		return(module_colors)
	})
	rownames(module_cols) = rownames(tissue_counts)
	return(module_cols)
}

comp_module_cols = colorByPurity(comp_stats$tissue_counts, cols)
bw_module_cols = colorByPurity(bw_stats$tissue_counts, cols)
sing_module_cols = colorByPurity(sing_stats$tissue_counts, cols)

# comp_module_cols = apply(comp_stats$tissue_counts, 2, function(counts) {
# 	tissue_present = counts / sum(counts) > 0.05
# 	module_colors = rep("white", length(tissue_present))  # default
# 	module_colors[tissue_present] = cols[tissue_present]
# 	return(module_colors)
# })

# bw_module_cols = apply(bw_stats$tissue_counts >= 1, 2, function(tissue_present) {
# 	module_colors = rep("white", length(tissue_present))  # default
# 	module_colors[tissue_present] = cols[tissue_present]
# 	return(module_colors)
# })

# sing_module_cols = apply(bw_stats$tissue_counts >= 1, 2, function(tissue_present) {
# 	module_colors = rep("white", length(tissue_present))  # default
# 	module_colors[tissue_present] = cols[tissue_present]
# 	return(module_colors)
# })



pdf("co-expression/plotsModuleStats/jaccard_bw_comp.pdf", width=10, height=10)
heatmap.3(simil_bw_comp,
	trace="none",
	ColSideColors=t(comp_module_cols),
	RowSideColors=bw_module_cols[nrow(bw_module_cols):1, ],
	ColSideColorsSize=3,
	RowSideColorsSize=3,
	labRow="", labCol="",
	main="Module overlap",
	xlab="Between-within", ylab="Complete",
	dendrogram="none",
	col=colorRampPalette(brewer.pal(9, "Greens"))(100),
	breaks=seq(0, 0.8, length.out=100 + 1)
	)
dev.off()

pdf("co-expression/plotsModuleStats/jaccard_bw_sing.pdf", width=10, height=10)
heatmap.3(simil_bw_sing,
	trace="none",
	ColSideColors=t(sing_module_cols),
	RowSideColors=bw_module_cols[nrow(bw_module_cols):1, ],
	ColSideColorsSize=2,
	RowSideColorsSize=2,
	labRow="", labCol="",
	main="Module overlap",
	xlab="Single", ylab="Between-within",
	dendrogram="none",
	col=colorRampPalette(brewer.pal(9, "Greens"))(100),
	breaks=seq(0, 0.8, length.out=100 + 1)
)
dev.off()

pdf("co-expression/plotsModuleStats/jaccard_sing_comp.pdf", width=10, height=10)
heatmap.3(simil_sing_comp,
	trace="none",
	ColSideColors=t(comp_module_cols),
	RowSideColors=sing_module_cols[nrow(sing_module_cols):1, ],
	ColSideColorsSize=2,
	RowSideColorsSize=2,
	labRow="", labCol="",
	main="Module overlap",
	xlab="Complete", ylab="Single",
	dendrogram="none",
	col=colorRampPalette(brewer.pal(9, "Greens"))(100),
	breaks=seq(0, 0.8, length.out=100 + 1)
)
dev.off()



