rm(list=ls())

library(data.table)
library(caret)
library(gplots)
library(viridis)

library(compiler)
enableJIT(3)

# library(parallel)

setwd("~/GoogleDrive/projects/STARNET/cross-tissue")

# Load reference gene networks

# Containing unzipped network files
marbach_high_level_dir = "/Users/sk/DataBases/regulatorynetworks/Network_compendium/Tissue-specific_regulatory_networks_FANTOM5-v1/32_high-level_networks"
marbach_dir = "/Users/sk/DataBases/regulatorynetworks/FANTOM5_individual_networks/394_individual_networks"


loadMarbachNetw = function(file_name) {
	tab = fread(file_name)
	names(tab) = c("from", "to", "weight")

	# From-to edge IDs, for looking up if edges are present
	tab$edge_id = paste(tab$from, tab$to, sep="_")

	return(tab)
}

# netw_ref = lapply(file.path(marbach_high_level_dir, netw_files), function(file_name) {
# 	message(file_name)
# 	tab = loadMarbachNetw(file_name)
# 	return(tab)
# })
# names(netw_ref) = netw_files


# Load STARNET Bayesian networks
# ---------------------------
bnet = fread("co-expression/annotate/bayesNet/all.tsv")
# bnet_nodes = fread("co-expression/annotate/bayesNet/nodes.tsv")

# Some formatting
bnet$tail_symbol = sapply(strsplit(bnet$TAIL, "_"), function(x) x[2])
bnet$head_symbol = sapply(strsplit(bnet$HEAD, "_"), function(x) x[2])

bnet$edge_id = paste(bnet$tail_symbol, bnet$head_symbol, sep="_")


# Co-expression module table
# ------------------------------------
modules = fread("co-expression/tables/modules.csv")

modules$node_id = paste(modules$tissue, modules$transcript_id, sep="_")

# Check coverage of Bayesian network
# Note that not all modules have inferred networks due to being too large.
mean(modules$node_id %in% c(bnet$TAIL, bnet$HEAD))



# Sensitivity and precision analysis
cmpNetw = function(bnet, modules, mod_ids, rnetw, method="global") {
	# Test each STARNET subnetwork
	results = sapply(mod_ids, function(k) {
		idx = modules$clust == k
		node_ids = modules$node_id[idx]
		node_symbols = modules$gene_symbol[idx]

		if (method == "local") {
			sub_bnet = bnet[bnet$HEAD %in% node_ids & bnet$TAIL %in% node_ids, ] 
			sub_ref_netw = rnetw[rnetw$from %in% node_symbols & rnetw$to %in% node_symbols, ]
		} else if (method == "global") {
			sub_bnet = bnet[bnet$HEAD %in% node_ids | bnet$TAIL %in% node_ids, ]
			sub_ref_netw = rnetw[rnetw$from %in% node_symbols | rnetw$to %in% node_symbols, ]
		} else {
			stop("Unrecognized method: ", method)
		}

		sub_bnet$STARNET = TRUE
		sub_ref_netw$marbach = TRUE

		# Combine networks, keep all
		netw_merge = merge(
			sub_bnet[ , c("edge_id", "STARNET")],
			sub_ref_netw[, c("edge_id", "marbach", "weight")],  # V3 is the score
			by="edge_id", all=TRUE)

		netw_merge$STARNET[is.na(netw_merge$STARNET)] = FALSE  # not found in STARNET
		netw_merge$marbach[is.na(netw_merge$marbach)] = FALSE

		netw_merge$weight[is.na(netw_merge$weight)] = 0  # Marbach weight

		out = list()
		TP = sum(netw_merge$STARNET & netw_merge$marbach)  # true positive
		P = sum(netw_merge$marbach)  # positive
		FP = sum(netw_merge$STARNET & !netw_merge$marbach)

		# Sensitivity (validation rate)
		out$sensitivity =  TP / P
		out$precision = TP / (TP + FP)
		out$TP = TP

		return(out)
	})

	# Some formating
	results = t(results)
	rownames(results) = mod_ids

	results = data.matrix(as.data.frame(results))

	return(results)
}

# Convenient way of getting matrix from list of statistics
getStat = function(tests, stat) {
	i = which(colnames(netw_tests[[1]]) == stat)
	sensitivity = sapply(netw_tests, function(x) x[, i])
	sensitivity[is.na(sensitivity)] = 0

	sensitivity = data.matrix(data.frame(sensitivity))
	return(sensitivity)
}


# Get module IDs with inferred networks
mod_ids = unique(modules$clust[modules$node_id %in% c(bnet$TAIL, bnet$HEAD)])
mod_ids = sort(mod_ids)

netw_files = list.files(marbach_high_level_dir, pattern="*.txt")
netw_files = netw_files[netw_files != "_clusters.txt"]  # exclude

# Test high-level Marbach networks
# netw_tests = lapply(names(netw_ref), function(netw_name) {
netw_tests = lapply(netw_files, function(file_name) {
	message(netw_name)
	# rnetw = netw_ref[[netw_name]]

	rnetw = loadMarbachNetw(file.path(marbach_high_level_dir, file_name))

	# Loop over each reference network
	results = cmpNetw(bnet, modules, mod_ids, rnetw)

	return(results)
})
# names(netw_tests) = names(netw_ref)
names(netw_tests) = netw_files


# netw_file = "liver_adult.txt"

# Test all Marbach networks
# netw_files_all = list.files(marbach_dir)[1:6]
netw_files_all = list.files(marbach_dir)

netw_tests_all = lapply(netw_files_all, function(netw_file) {
	message(netw_file)
	rnetw = loadMarbachNetw(file.path(marbach_dir, netw_file))

	# Loop over each reference network
	results = cmpNetw(bnet, modules, mod_ids, rnetw)

	return(results)
})
names(netw_tests_all) = netw_files_all


# boxplot(results[, 1])

pt_size = results[, 3] * 0.1
x = results[, 1] * 100
y = results[, 2] * 100
plot(x, y,
	cex=pt_size,
	xlab="Sensitivity (%)", ylab="Specificity (%)")

text(x, y, legend=rownames(results), pos=1, cex=0.5)

boxplot(as.data.frame(results)$sensitivity)

boxplot(data.matrix(as.data.frame(results))[, 1] * 100)
boxplot(as.data.frame(results)[, 1] * 100)

[, 1]


sensitivity = getStat(netw_tests, "sensitivity")

apply(sensitivity, 2, max)


# data.frame(sensitivity)
# colnames(sensitivity) = names(netw_ref)

# sort(apply(sensitivity, 1, max))

heatmap.2(sensitivity * 100,
	trace="none",
	col=viridis(100)
)
