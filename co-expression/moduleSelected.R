options(java.parameters = "-Xmx8g")  # Max Java memory heap size

rm(list=ls())

library(RColorBrewer)
library(gplots)

library(devtools)  # for installing heatmap.3
# Load heatmap.3
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")

library(compiler)
enableJIT(3)

# library(bnlearn)
library(Mergeomics)
library(igraph)

# Note that installation of rcausal may require to run:
# $ sudo R CMD javareconf
# install.packages("rJava",type='source')
library(rcausal)


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

# Removes suffix of ensembl IDs
ensemblBase = function(ensembl_ids) {
	sapply(as.character(ensembl_ids), function(id) strsplit(id, "[.]")[[1]][1])
}


liverModuleOverlap = function(module, liver_module) {
	liver_overlap_n = sum(module$tissue == "LIV" & ensemblBase(module$ensembl) %in% ensemblBase(liver_module$ensembl))
	liver_missing_n = sum(module$tissue == "LIV" & !ensemblBase(module$ensembl) %in% ensemblBase(liver_module$ensembl))
	non_liver_n = sum(module$tissue != "LIV")

	counts = c(liver_overlap_n, liver_missing_n, non_liver_n)
	names(counts) = c("LIV overlap", "LIV new", "Non-LIV")
	return(counts)
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

# List-based reference to module statistics
mod_stats = list(
	"Between-within"=bw_stats,
	"Complete"=comp_stats,
	"Single"=sing_stats)

# Default color scheme
cols = brewer.pal(9, "Set1")[c(1:5, 7:9)]



# Load gene expession data
# --------------------------------
data_dir = "~/DataProjects/cross-tissue"
emat_file = "STARNET/gene_exp_norm_reshape/expr_recast.RData"

load(file.path(data_dir, emat_file), verbose=TRUE)

emat = expr_recast[, 3:ncol(expr_recast)]
meta_genes = expr_recast[, 1:2]
meta_genes = as.data.frame(meta_genes)

# tissue_transcript IDs
meta_genes$id = paste(meta_genes$tissue, meta_genes$transcript_id, sep="_")

rm(expr_recast)


# Compare modules to Ariella's liver module
# --------------------------------------------------------

# Load liver module gene ensembl IDs
liver_module = read.table(file.path(data_dir, "ariella_modules/gold_ens_transcriptIDs"))
colnames(liver_module) = "ensembl"

# Find genes in liver in each set of modules
idx_between = which(between_within$meta_genes$tissue == "LIV" & 
	ensemblBase(between_within$meta_genes$ensembl) %in% ensemblBase(liver_module$ensembl))
between_within$clust[idx_between]

idx_complete = which(complete$meta_genes$tissue == "LIV" &
	ensemblBase(complete$meta_genes$ensembl) %in% ensemblBase(liver_module$ensembl))
complete$clust[idx_complete]

idx_single = which(sing$meta_genes$tissue == "LIV" &
	ensemblBase(sing$meta_genes$ensembl) %in% ensemblBase(liver_module$ensembl))
sing$clust[idx_single]


pdf("co-expression/plotsModuleStats/liver_module_match.pdf", width=3.0, height=6)
par(mfrow=c(3, 1))

barplot(table(between_within$clust[idx_between]),
	main="Within-between",
	ylab="Liver module genes",
	xlab="Module",
	col=cols[1])

barplot(table(complete$clust[idx_complete]),
	main="Complete",
	ylab="Liver module genes",
	xlab="Module",
	col=cols[2])

barplot(table(sing$clust[idx_single]),
	main="Single",
	ylab="Liver module genes",
	xlab="Module",
	col=cols[3])

dev.off()


# Tissue distribution of genes from liver module, based on identified idx
# --------------------------------------------------------------
liver_modules_between = data.frame(
	sort(table(between_within$clust[idx_between]), decreasing=TRUE))
colnames(liver_modules_between)[1] = "clust"

liver_modules_complete = data.frame(
	sort(table(complete$clust[idx_complete]), decreasing=TRUE))
colnames(liver_modules_complete)[1] = "clust"

liver_modules_single = data.frame(
	sort(table(sing$clust[idx_single]), decreasing=TRUE))
colnames(liver_modules_single)[1] = "clust"


bw_stats$tissue_counts[, as.integer(as.character(liver_modules_between$clust))]
comp_stats$tissue_counts[, as.integer(as.character(liver_modules_complete$clust))]
sing_stats$tissue_counts[, as.integer(as.character(liver_modules_single$clust))]



clust = as.integer(as.character(liver_modules_between$clust[1]))
liver_module_between = between_within$meta_genes[between_within$clust == clust ,]

clust = as.integer(as.character(liver_modules_complete$clust[1]))
liver_module_complete = complete$meta_genes[complete$clust == clust ,]

clust = as.integer(as.character(liver_modules_single$clust[1]))
liver_module_single = sing$meta_genes[sing$clust == clust ,]


pdf("co-expression/plotsModuleStats/liver_module_match_sel.pdf", width=3, height=4)
bar_cols = c(gray.colors(3)[c(1, 1)], brewer.pal(9, "Set1")[5])
barplot(
	cbind(liverModuleOverlap(liver_module_between, liver_module),
		liverModuleOverlap(liver_module_complete, liver_module),
		liverModuleOverlap(liver_module_single, liver_module)),
	names.arg=c(98, 120, 298),
	ylab="Cross-tissue module genes",
	density=c(NA, 15, NA),
	col=bar_cols)
dev.off()


# print new genes identified in LIV
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


# Find modules associated with specific genes.
#-----------------------------------------------------------------

gene = "KIAA1462"  # gene to look up

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

# Create directory for selected gene
sel_gene_dir = file.path("co-expression/plotsModuleSelected", gene)
dir.create(sel_gene_dir)


# Plot tissue distribution for identified modules <2000 genes. Stacked barplot.
pdf(file.path(sel_gene_dir, "module_sizes.pdf"), width=2.0)
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

	# Reorder based on size
	tissue_counts = tissue_counts[, order(apply(tissue_counts, 2, sum), decreasing=TRUE)]

	barplot(tissue_counts,
		# names.arg=names(modules)[include_modules],
		names.arg=colnames(tissue_counts),
		main=mod_env$opts$method,
		ylab="Module size",
		border=NA,
		las=2,
		col=cols)
	abline(h=0)
	legend("topright",
		legend=rev(rownames(tissue_counts)),
		col=rev(cols[1:nrow(tissue_counts)]),
		pch=15, cex=0.6)
}
dev.off()


# Within-between modules,
modules = findModule(between_within, gene)

# Write modules to files
module_dir = file.path(sel_gene_dir, "between_within")
dir.create(module_dir)
for (mod_name in names(modules)) {
	write.csv(modules[[which(names(modules) == mod_name)]],
		file=file.path(module_dir, paste0("module_", mod_name, ".csv")),
		row.names=FALSE,
		quote=FALSE)
}



# Selected modules from other methods
modules = findModule(complete, gene)

module_dir = file.path(sel_gene_dir, "complete")
dir.create(module_dir)
for (mod_name in names(modules)) {
	write.csv(modules[[which(names(modules) == mod_name)]],
		file=file.path(module_dir, paste0("module_", mod_name, ".csv")),
		row.names=FALSE,
		quote=FALSE)
}

# mod_name = 104
# write.csv(modules[[which(names(modules) == mod_name)]],
# 	file=paste0("co-expression/plotsModuleSelected/KIAA1462/module_complete_", mod_name, ".csv"),
# 	quote=FALSE)


# Learn Bayesian networks for each module based on standardized expression data.
# Input: 
#		modules is a list of module (data.frames with columbs)
#   	emat is a tissue:transcript expression matrix with corresponding meta data
learnBayesNets = function(modules, meta_genes, emat) {
	if (nrow(meta_genes) != nrow(emat)) {
		stop("Expression matrix and transcript meta data mismatch.")
	}

	if (! "id" %in% colnames(meta_genes)) {
		stop("Transcript meta data does not contain id column")
	}

	bayes_nets = lapply(modules, function(mod) {
		# Tissue_transcipt IDs for module
		mod_ids = paste(mod$tissue, mod$transcript_id, sep="_")
		mod_gene_symbols = mod$gene_symbol

		# Print first 10 transcript IDs
		message("Module: ", paste(mod_ids[1:2], collapse=", "), "...(n=", length(mod_ids), ")")

		idx = match(mod_ids, meta_genes$id)

		# Get matrix
		submat = t(emat[idx, ])
		colnames(submat) = mod_ids

		# Standardize and impute missing to mean
		submat = scale(submat)
		submat[is.na(submat)] = 0.0
		submat = as.data.frame(submat)

		# bn = hc(submat)  # greedy hill climbing, BIC score
		bn = fges(submat, maxDegree=100)  # Fast-greedy equivalence search

		return(bn)
	})

	return(bayes_nets)
}



modules = findModule(between_within, gene)
module_dir = file.path(sel_gene_dir, "between_within")  # root module method dir

modules = findModule(complete, gene)
module_dir = file.path(sel_gene_dir, "complete")  # root module method dir

# Exclude modules with too many genes
modules = modules[
	sapply(modules, function(mod) nrow(mod)) < 2000
]

# lapply(modules, dim)


# Learn Bayesian network within each module
# WARNING: long runtime for networks with more than 1000 nodes
bayes_nets = learnBayesNets(modules, meta_genes, emat)

# Get all edge lists, for all Bayesian networks
# edge_list = sapply(bayes_nets, function(net) net$arcs)

edge_list = lapply(bayes_nets, function(net) {
	g = graph_from_graphnel(net$graphNEL)
	return(as_edgelist(g))
})

# Combine edge list
edge_list = Reduce(rbind, edge_list)

edge_list = cbind(edge_list, 1)
colnames(edge_list) = c("TAIL", "HEAD", "WEIGHT")
edge_list = as.data.frame(edge_list)

# write.csv(edge_list,
# 	file.path(module_dir, "BayesNet/all.csv"),
# 	row.names=FALSE,
# 	quote=FALSE
# )

dir.create(file.path(module_dir, "BayesNet"))
write.table(edge_list,
	sep="\t",
	file.path(module_dir, "BayesNet/all.tsv"),
	row.names=FALSE,
	quote=FALSE
)

# module table
mod_table = lapply(seq_along(modules), function(i) {
	mod = modules[[i]]  # module reference
	mod_name = names(modules)[i]  # name of module
	# get transcript IDs, tissue_mRNA
	mod_gene_ids = paste(mod$tissue, mod$transcript_id, sep="_")

	df = data.frame(MODULE=mod_name, NODE=mod_gene_ids)
	return(df)
})
mod_table = Reduce(rbind, mod_table)

write.table(mod_table,
	file.path(module_dir, "nodes.tsv"),
	sep="\t",
	row.names=FALSE,
	quote=FALSE
)

# Weighted key driver analysis (wKDA) of Bayesian networks within each cross-tissue module
kda_label = paste0(gene, "_key_driver")
job.kda = list()
job.kda$netfile = file.path(module_dir, "BayesNet", "all.tsv")
job.kda$modfile = file.path(module_dir, "nodes.tsv")
job.kda$label = kda_label
job.kda$folder = module_dir
job.kda$nperm = 20

job.kda = kda.configure(job.kda)
job.kda = kda.start(job.kda)
job.kda = kda.prepare(job.kda)
job.kda = kda.analyze(job.kda)
job.kda = kda.finish(job.kda)

# load results table
kda_results = read.table(
	file.path(module_dir, "kda", paste0(kda_label, ".results.txt")),
	header=TRUE
)

kda_results$gene_symbol = sapply(strsplit(as.character(kda_results$NODE), "_"), function(x) x[2])
kda_results$tissue = sapply(strsplit(as.character(kda_results$NODE), "_"), function(x) x[1])


# Selective plot of key drivers
# -----------------------------------------------
kda_sel = kda_results[kda_results$MODULE == 155, ]

kda_sel = kda_sel[kda_sel$FDR < 0.05, ]

write.table(kda_sel$gene_symbol,
	file=file.path(module_dir, "kda", "mod_155_genes.txt"),
	quote=FALSE, row.names=FALSE
)



n = 30

i = which(kda_sel$gene_symbol == gene)

x = kda_sel$FDR[c(1:n, i)]

names(x) = kda_sel$gene_symbol[c(1:n, i)]
# bar_col

pdf(file.path(module_dir, "kda", "mod_155_key_drivers.pdf"), width=3)
par(mar=c(4, 6, 4, 4))
barplot(-log10(rev(x)),
	las=2,
	main="Module 155 key drivers",
	xlab=expression("log"[10] * " p"),
	horiz=TRUE)

abline(v=0)
dev.off()

# Make subnetwork containing only identified key drivers
# -----------------------------------------------
alpha = 0.000001

alpha = 0.00000001

include_nodes = kda_results$NODE[kda_results$FDR < alpha]
# include_nodes = as.character(include_nodes)

sub_edge_list = edge_list[
	edge_list$TAIL %in% include_nodes |
	edge_list$HEAD %in% include_nodes
, ]

dim(sub_edge_list)


write.csv(sub_edge_list,
	file.path(module_dir, "BayesNet", paste0("key_driver", ".csv")),
	row.names=FALSE,
	quote=FALSE
)


# Make subnetwork involving target gene
# -----------------------------------------
from_gene = sapply(
	strsplit(as.character(edge_list$TAIL), "_"),
	function(x) x[2])

to_gene = sapply(
	strsplit(as.character(edge_list$HEAD), "_"),
	function(x) x[2])

sub_edge_list = edge_list[gene == from_gene | gene == to_gene, ]

# convert to igraph object
g = graph_from_edgelist(as.matrix(sub_edge_list[,1:2]))
V(g)$tissue = sapply(strsplit(names(V(g)), "_"), function(x) x[1])
V(g)$gene_symbol = sapply(strsplit(names(V(g)), "_"), function(x) x[2])
V(g)$module = as.character(mod_table$MODULE)[match(names(V(g)), mod_table$NODE)]

# Write edge list
write.csv(sub_edge_list,
	file.path(module_dir, "BayesNet", paste0("target.csv")),
	row.names=FALSE,
	quote=FALSE
)

write_graph(g, file=file.path(module_dir, "BayesNet", "target.gml"), format="gml")


# -------------------------------
modules = findModule(between_within, gene)
names(modules)
lapply(modules, nrow)


i = 3

mod_ids = paste(modules[[i]]$tissue, modules[[i]]$transcript_id, sep="_")
mod_gene_symbols = modules[[i]]$gene_symbol

idx = match(mod_ids, meta_genes$id)

# Get matrix
submat = t(emat[idx, ])
colnames(submat) = mod_ids

submat = scale(submat)
submat[is.na(submat)] = 0.0

submat = as.data.frame(submat)

# Learn Bayesian network structure
# bn = mmhc(submat)
# bn = hc(submat)  # greedy hill climbing, BIC

bn = fges(submat, maxDegree=100)  # Fast-greedy search

# bn$edges

# Convert to igraph object
g = graph_from_graphnel(bn$graphNEL)

degree_distribution(g)


edge_list = as_edgelist(g)

# edge list
edge_list = cbind(bn$arcs, 1)
colnames(edge_list)[3] = "weight"

# Get edge list
write.csv(edge_list,
	file.path(module_dir, "BayesNet", paste0(names(modules)[i], ".csv")),
	row.names=FALSE,
	quote=FALSE
)






# Correlation matrix
cmat = cor(t(emat[idx, ]), use="pairwise.complete.obs")
rownames(cmat) = mod_gene_symbols
colnames(cmat) = mod_gene_symbols

heatmap.2(cmat,
	trace="none",
	col=colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
	cexRow=0.25, cexCol=0.25
)

paste(meta_genes$tissue, meta

head(meta_genes)

