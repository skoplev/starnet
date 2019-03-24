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
library(RColorBrewer)

library(compiler)
enableJIT(3)

setwd("~/GoogleDrive/projects/STARNET/cross-tissue")
data_dir = "/Users/sk/DataProjects/cross-tissue"  # root of data directory

source("src/models/bayesNet.R")
source("src/parse.R")
source("co-expression/annotate/lib/netwValidate.R")


# Load cross-tissue modules
# ------------------------------
mod_tab = fread("co-expression/tables/modules.csv")


# Load normalized gene expression data, used for evaluating Bayesian networks within cross-tissue modules
# --------------------------------------------------------------
emat_file = "STARNET/gene_exp_norm_reshape/expr_recast.RData"

load(file.path(data_dir, emat_file), verbose=TRUE)

expr = parseExprTable(expr_recast)

rm(expr_recast)

# Test if the gene metadata is the same as for the cross-tissue modules
if (!all(expr$meta_row$transcript_id == mod_tab$transcript_id)) {
	stop("Transcript mismatch")
}


# Load cis-trans-eQTL network
# ---------------------------------------
cis_eqtl_dir = "/Users/sk/DataProjects/STARNET/eQTL"
cis_eqtl_files = list.files(cis_eqtl_dir, pattern="*.cis.tbl")

trans_eqtl_dir = "/Users/sk/DataProjects/STARNET/eQTL/trans"
trans_eqtl_files = list.files(trans_eqtl_dir)

# Load cis-eQTL
cis_eqtl = lapply(
	cis_eqtl_files,
	function(file_name) {
		d = fread(file.path(cis_eqtl_dir, file_name))
		d = d[d$adj.p < 0.05, ]  # FDR < 5%
		d = d[order(d$p), ]
		return(d)
})
names(cis_eqtl) = sapply(strsplit(cis_eqtl_files, "[.]"), function(x) x[4])  # name by tissue code


# Load trans-eQTL
trans_eqtl = lapply(
	trans_eqtl_files,
	function(file_name) {
		d = fread(file.path(trans_eqtl_dir, file_name))
		colnames(d) = colnames(cis_eqtl[[1]])[-6]  # no adjusted p-values
		return(d)
})
names(trans_eqtl) = sapply(strsplit(trans_eqtl_files, "[.]"), function(x) x[4])


# Check if tissue names agree
stopifnot(all(names(cis_eqtl) == names(trans_eqtl)))


cis_trans_eqtl = lapply(1:length(cis_eqtl), function(i) {
	# join. x: cis, y: trans
	d = merge(cis_eqtl[[i]], trans_eqtl[[i]],
		by="SNP"
	)

	return(d)
})
names(cis_trans_eqtl) = names(cis_eqtl)

# Add tissue column
for (i in 1:length(cis_trans_eqtl)) {
	cis_trans_eqtl[[i]]$tissue = toupper(names(cis_trans_eqtl)[i])
}


# All cis-trans-eQTL, without macrophages and foam cells
cis_trans_eqtl_all = rbindlist(cis_trans_eqtl[!names(cis_trans_eqtl) %in% c("MP", "FC")])


# Load Marbach networks
# --------------------------
marbach_dir = "/Users/sk/DataBases/regulatorynetworks/FANTOM5_individual_networks/394_individual_networks"

netw_files_all = list.files(marbach_dir)
marbach_netw = lapply(netw_files_all, function(file_name) {
	message(file_name)
	rnetw = loadMarbachNetw(file.path(marbach_dir, file_name))

	# Filter for high-confidence interactions
	rnetw = rnetw[rnetw$weight > 0.1, ] 
	return(rnetw)
})

marbach_netw_all = rbindlist(marbach_netw)

# Sort by weight
marbach_netw_all = marbach_netw_all[order(marbach_netw_all$weight, decreasing=TRUE), ]

# Only unique assocations, highest weight across networks
marbach_netw_all = marbach_netw_all[match(unique(marbach_netw_all$edge_id), marbach_netw_all$edge_id), ]


# Bayesian networks
# --------------------------------------

bayes_dir = "co-expression/annotate/bayesNet3"

# Infer Bayeisan networks withing cross-tissue modules
bayes_nets = learnBayesNets.4(
	clust=mod_tab$clust,
	meta_genes=expr$meta_row,
	emat=expr$mat,
	ref_netw=marbach_netw_all,
	eqtl_netw=cis_trans_eqtl_all,
	max_size=3000)

# save(bayes_nets, file=file.path(bayes_dir, "bayes_nets_marbach01.RData"))


# Plots results
# ----------------------------------------

edge_stats = sapply(bayes_nets, function(x) {
	if (is.na(x)) {
		values = rep(NA, 3)
	} else {
		n1 = nrow(x$ref_netw_include)
		n2 = nrow(x$eqtl_netw_include)
		n3 = length(x$edges)
		values = c(n1, n2, n3 - n2 - n1)
		names(values) = c("Known TF-gene", "cis-trans-eQTL", "Inferred")
	}
	return(values)
})


pdf("co-expression/annotate/bayesNet3/plots/bayes_netw_edges_marbach01.pdf", width=5, height=4)
idx = order(apply(edge_stats, 2, sum), decreasing=TRUE)
barplot(
	edge_stats[, idx],
	xlab="Co-expression module",
	ylab="Gene-gene interactions",
	border=NA,
	space=-0.3,
	col=brewer.pal(9, "Set1")[c(5, 1, 2)],
	legend.text=TRUE)
dev.off()

# lapply(bayes_nets, function(x) {
# 	if (!is.na(x))
# 		x$eqtl_netw_include
# })
