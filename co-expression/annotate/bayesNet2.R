rm(list=ls())

options(java.parameters = "-Xmx8g")  # Max Java memory heap size, for rcausal

library(withr)
library(devtools)

# Local install, don't update java etc...
with_libpaths(new="/Users/sk/GoogleDrive/projects/STARNET/cross-tissue/lib/R", install_github('bd2kccd/r-causal'))
with_libpaths(new="/Users/sk/GoogleDrive/projects/STARNET/cross-tissue/lib/R", library('rcausal'))


library(data.table)
library(igraph)
# library(rcausal)
library(bnlearn)
library(Mergeomics)
# library(plyr)

library(compiler)
enableJIT(3)

setwd("~/GoogleDrive/projects/STARNET/cross-tissue")
data_dir = "/Users/sk/DataProjects/cross-tissue"  # root of data directory

source("src/models/bayesNet.R")
source("src/parse.R")


# Load cross-tissue modules
# ------------------------------
between = new.env()
load(file.path(data_dir, "modules/between_within-cross-tissue.RData"),
	between,
	verbose=TRUE)

between = parseModuleData(between)


# Load normalized gene expression data, used for evaluating Bayesian networks within cross-tissue modules
# --------------------------------------------------------------
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


# Load TF definition from Lambert et al
# -------------------------------------------------------

tf_symbols = as.character(read.table("transcription-factors/lambert/TF_names_v_1.01.txt")$V1)


# Bayesian network for each cross-tissue module, Key driver analysis.
# -----------------------------------------------------------

bayes_dir = "co-expression/annotate/bayesNet2"

# Infer Bayeisan networks withing cross-tissue modules
bayes_nets = learnBayesNets.3(
	clust=between$clust,
	meta_genes=between$meta_genes,
	emat=emat,
	source_symbols=tf_symbols,
	max_size=3000)

edge_list = lapply(bayes_nets, function(net) {
	if (is.na(net)) {
		return(NA)
	} else {
		g = graph_from_graphnel(net$graphNEL)
		return(as_edgelist(g))
	}
})
names(edge_list) = 1:length(edge_list)  # cluster ids

# Remove entries for large modules without bayesian networks
edge_list = edge_list[!is.na(edge_list)]

# Combine edge list
edge_list = Reduce(rbind, edge_list)


edge_list = cbind(edge_list, 1)  # weight
colnames(edge_list) = c("TAIL", "HEAD", "WEIGHT")
edge_list = as.data.frame(edge_list)


# Write edge list to file, used for input to Weighted Key Driver analysis.
write.table(edge_list,
	file.path(bayes_dir, "all.tsv"),
	sep="\t",
	row.names=FALSE,
	quote=FALSE
)

# Make module table for KDA.
kda_mod_tab = data.frame(
	MODULE=between$clust,
	NODE=paste(between$meta_genes$tissue, between$meta_genes$transcript_id, sep="_")
)

write.table(kda_mod_tab,
	file.path(bayes_dir, "nodes.tsv"),
	sep="\t",
	row.names=FALSE,
	quote=FALSE
)

# Weighted key driver analysis (wKDA) of Bayesian networks within each cross-tissue module
# kda_label = paste0(gene, "_key_driver")
kda_label = "modules"
job.kda = list()
job.kda$netfile = file.path(bayes_dir, "all.tsv")
job.kda$modfile = file.path(bayes_dir, "nodes.tsv")
job.kda$label = kda_label
job.kda$folder = bayes_dir
job.kda$nperm = 20

job.kda = kda.configure(job.kda)
job.kda = kda.start(job.kda)
job.kda = kda.prepare(job.kda)
job.kda = kda.analyze(job.kda)
job.kda = kda.finish(job.kda)


# load results table
kda_results = read.table(
	file.path(bayes_dir, "kda", paste0(kda_label, ".results.txt")),
	header=TRUE
)
