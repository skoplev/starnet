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


# Compare modules to Ariella's liver module
# --------------------------------------------------------
liver_module = read.table(file.path(data_dir, "ariella_modules/gold_ens_transcriptIDs"))
colnames(liver_module) = "ensembl"

# Identify clusters of 


# liver_module$ensembl %in% complete$meta_genes$ensembl

ensemblBase = function(ensembl_ids) {
	sapply(as.character(ensembl_ids), function(id) strsplit(id, "[.]")[[1]][1])
}

# ensemblBase(liver_module$ensembl)

# liver_module$ensembl %in% complete$meta_genes$ensembl

# ensemblBase(liver_module$ensembl) %in% ensemblBase(complete$meta_genes$ensembl)

pdf("co-expression/plotsModuleStats/liver_module_match.pdf", width=3.0, height=6)
par(mfrow=c(3, 1))

idx_between = which(between_within$meta_genes$tissue == "LIV" & 
	ensemblBase(between_within$meta_genes$ensembl) %in% ensemblBase(liver_module$ensembl))
between_within$clust[idx_between]
barplot(table(between_within$clust[idx_between]),
	main="Within-between",
	ylab="Liver module genes",
	xlab="Module",
	col=cols[1])

liver_modules_between = data.frame(
	sort(table(between_within$clust[idx_between]), decreasing=TRUE))
colnames(liver_modules_between)[1] = "clust"


idx_complete = which(complete$meta_genes$tissue == "LIV" &
	ensemblBase(complete$meta_genes$ensembl) %in% ensemblBase(liver_module$ensembl))
complete$clust[idx_complete]
barplot(table(complete$clust[idx_complete]),
	main="Complete",
	ylab="Liver module genes",
	xlab="Module",
	col=cols[2])

liver_modules_complete = data.frame(
	sort(table(complete$clust[idx_complete]), decreasing=TRUE))
colnames(liver_modules_complete)[1] = "clust"



idx_single = which(sing$meta_genes$tissue == "LIV" &
	ensemblBase(sing$meta_genes$ensembl) %in% ensemblBase(liver_module$ensembl))
sing$clust[idx_single]
barplot(table(sing$clust[idx_single]),
	main="Single",
	ylab="Liver module genes",
	xlab="Module",
	col=cols[3])

liver_modules_single = data.frame(
	sort(table(sing$clust[idx_single]), decreasing=TRUE))
colnames(liver_modules_single)[1] = "clust"

dev.off()


# Tissue distribution of identified modules


bw_stats$tissue_counts[, as.integer(as.character(liver_modules_between$clust))]
comp_stats$tissue_counts[, as.integer(as.character(liver_modules_complete$clust))]
sing_stats$tissue_counts[, as.integer(as.character(liver_modules_single$clust))]

# bw_stats$tissue_counts[, liver_modules_between$clust]


clust = as.integer(as.character(liver_modules_between$clust[1]))
liver_module_between = between_within$meta_genes[between_within$clust == clust ,]


clust = as.integer(as.character(liver_modules_complete$clust[1]))
liver_module_complete = complete$meta_genes[complete$clust == clust ,]

clust = as.integer(as.character(liver_modules_single$clust[1]))
liver_module_single = sing$meta_genes[sing$clust == clust ,]


liverModuleOverlap = function(module) {
	liver_overlap_n = sum(module$tissue == "LIV" & ensemblBase(module$ensembl) %in% ensemblBase(liver_module$ensembl))
	liver_missing_n = sum(module$tissue == "LIV" & !ensemblBase(module$ensembl) %in% ensemblBase(liver_module$ensembl))
	non_liver_n = sum(module$tissue != "LIV")

	counts = c(liver_overlap_n, liver_missing_n, non_liver_n)
	names(counts) = c("LIV overlap", "LIV new", "Non-LIV")
	return(counts)
}

# par(mfrow=c(3, 1))
# # bar_cols = brewer.pal(9, "Paired")[c(2, 2, 4)]
# # bar_cols = brewer.pal(9, "Set1")[c(4, 4, 5)]
# # bar_cols = gray.colors(3)[c(1, 1, 2)]
# bar_cols = c(gray.colors(3)[c(1, 1)], brewer.pal(9, "Set1")[5])
# barplot(as.matrix(liverModuleOverlap(liver_module_between)),
# 	ylab="Cross-tissue module genes",
# 	density=c(NA, 15, NA),
# 	col=bar_cols)
# barplot(as.matrix(liverModuleOverlap(liver_module_complete)),
# 	ylab="Cross-tissue module genes",
# 	density=c(NA, 15, NA),
# 	col=bar_cols)
# barplot(as.matrix(liverModuleOverlap(liver_module_single)),
# 	ylab="Cross-tissue module genes",
# 	density=c(NA, 15, NA),
# 	col=bar_cols)

pdf("co-expression/plotsModuleStats/liver_module_match_sel.pdf", width=3, height=4)
barplot(
	cbind(liverModuleOverlap(liver_module_between),
		liverModuleOverlap(liver_module_complete),
		liverModuleOverlap(liver_module_single)),
	names.arg=c(98, 120, 298),
	ylab="Cross-tissue module genes",
	density=c(NA, 15, NA),
	col=bar_cols)
dev.off()


liver_module_between[
	liver_module_between$tissue == "LIV" &
	!ensemblBase(liver_module_between$ensembl) %in% ensemblBase(liver_module$ensembl)
	, ]

getNewLivGenes = function(module) {
	new_genes = module[
		module$tissue == "LIV" &
		!ensemblBase(module$ensembl) %in% ensemblBase(liver_module$ensembl)
		, ]
	return(new_genes)
}

g1 = getNewLivGenes(liver_module_between)$gene_symbol
g2 = getNewLivGenes(liver_module_complete)$gene_symbol
g3 = getNewLivGenes(liver_module_single)$gene_symbol

unique(c(g1, g2, g3))

intersect(g1, intersect(g2, g3))


# Check specific genes
#------------------------------------------

# Returns list of transcripts in 
findModule = function(mod_env, gene) {
	module_ids = mod_env$clust[
		which(mod_env$meta_genes$gene_symbol == gene)]

	modules = lapply(module_ids, function(id) {
		mod_env$meta_genes[mod_env$clust == id, ]
	})

	names(modules) = module_ids

	return(modules)
}

gene = "KIAA1462"  # gene to look up

pdf("co-expression/plotsModuleSelected/KIAA1462/module_sizes.pdf", width=3)
par(mfrow=c(3, 1))
for (mod_env in c(between_within, complete, sing)) {
	modules = findModule(mod_env, gene)

	lapply(modules, dim)

	tissue_counts = sapply(modules, function(mod) {
		table(mod$tissue)
	})

	# Exclude modules containing more than 2000 genes
	include_modules = apply(tissue_counts, 2, sum) < 2000
	tissue_counts = tissue_counts[, include_modules]

	barplot(tissue_counts, names.arg=names(modules)[include_modules],
		main=mod_env$opts$method,
		ylab="Module size",
		border=NA,
		col=cols)
	legend("topleft", legend=rownames(tissue_counts),
		col=cols,
		pch=15, cex=0.6)
}
dev.off()


modules = findModule(between_within, gene)
# mod_name = 20
# mod_name = 142
# mod_name = 155
mod_name = 82
write.csv(modules[[which(names(modules) == mod_name)]],
	file=paste0("co-expression/plotsModuleSelected/KIAA1462/module_between_within_", mod_name, ".csv"),
	quote=FALSE)


modules = findModule(complete, gene)
mod_name = 104
write.csv(modules[[which(names(modules) == mod_name)]],
	file=paste0("co-expression/plotsModuleSelected/KIAA1462/module_complete_", mod_name, ".csv"),
	quote=FALSE)

# modules = findModule(complete, gene)
# modules = findModule(sing, gene)

