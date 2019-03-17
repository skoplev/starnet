rm(list=ls())

library(data.table)
library(caret)
library(gplots)
library(viridis)
library(RColorBrewer)

library(compiler)
enableJIT(3)

# library(parallel)

setwd("~/GoogleDrive/projects/STARNET/cross-tissue")

source("co-expression/annotate/lib/netwValidate.R")

# Load reference gene networks

# Containing unzipped network files
marbach_high_level_dir = "/Users/sk/DataBases/regulatorynetworks/Network_compendium/Tissue-specific_regulatory_networks_FANTOM5-v1/32_high-level_networks"
marbach_dir = "/Users/sk/DataBases/regulatorynetworks/FANTOM5_individual_networks/394_individual_networks"


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


# Get module IDs with inferred networks
mod_ids = unique(modules$clust[modules$node_id %in% c(bnet$TAIL, bnet$HEAD)])
mod_ids = sort(mod_ids)

netw_files = list.files(marbach_high_level_dir, pattern="*.txt")
netw_files = netw_files[netw_files != "_clusters.txt"]  # exclude

# Test high-level Marbach networks
# netw_tests = lapply(names(netw_ref), function(netw_name) {
netw_tests = lapply(netw_files, function(file_name) {
	message(file_name)
	# rnetw = netw_ref[[netw_name]]

	rnetw = loadMarbachNetw(file.path(marbach_high_level_dir, file_name))

	# Loop over each reference network
	results = cmpNetw(bnet, modules, mod_ids, rnetw)

	return(results)
})
names(netw_tests) = netw_files


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

# Write output to file
saveRDS(netw_tests, file="co-expression/annotate/netwValidation/bayesnet1MarbachHighLevel.rds")
saveRDS(netw_tests_all, file="co-expression/annotate/netwValidation/bayesnet1Marbach.rds")


# Visualize results
# -------------------------------------------------

sensitivity = getStat(netw_tests_all, "sensitivity")
precision = getStat(netw_tests_all, "precision")
TP = getStat(netw_tests_all, "TP")

plot(sensitivity * 100, precision * 100)


# Equal weight optimal match, get statistics
opt_idx = apply(sensitivity + precision, 1, which.max)

# Get matrix coordinates of optimal match
opt_mat_idx = cbind(
	1:length(opt_idx),
	opt_idx)

# Collect data frame of optimal validation statistics
netw_val = data.frame(
	module=rownames(sensitivity),
	network_match=colnames(sensitivity)[opt_idx],
	sensitivity=sensitivity[opt_mat_idx],
	precision=precision[opt_mat_idx],
	TP=TP[opt_mat_idx]
)
netw_val = netw_val[order(netw_val$precision, decreasing=TRUE), ]



# Plot from optimal precision-sensitivity networks
pdf("co-expression/annotate/netwValidation/plots/bnet1_maxPrecisionSensitivity.pdf", height=5, width=4.5)
pt_col = rep("black", nrow(netw_val))

colors = brewer.pal(9, "Set1")

highlight_netw = names(sort(table(netw_val$network_match), decreasing=TRUE)[2:4])
for (i in 1:length(highlight_netw)) {
	pt_col[netw_val$network_match == highlight_netw[i]] = brewer.pal(9, "Set1")[i]
}

plot(netw_val$sensitivity * 100, netw_val$precision * 100,
	main="BayesNet vs Marbach, best per module",
	xlab="Sensitivity (%)",
	ylab="Precision (%)",
	col=pt_col,
	pch=16)

legend("topright", legend=highlight_netw, col=colors, pch=16)
dev.off()




pdf("co-expression/annotate/netwValidation/plots/bnet1_precision_boxplots.pdf", height=4, width=14.0)
idx = order(apply(precision, 1, max), decreasing=TRUE)
# par(cex.axis=0.5)
boxplot(t(precision[idx, ]) * 100,
	frame=FALSE,
	cex.axis=0.5,
	las=2,
	pch=16,
	cex=0.5,
	ylab="Precision (%)")
dev.off()



heatmap.2(
	# sensitivity * 100,
	precision * 100,
	trace="none",
	col=viridis(100)
)


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

