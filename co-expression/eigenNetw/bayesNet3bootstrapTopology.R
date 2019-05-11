#Supernetwork bootstrap and topology analysis

options(java.parameters = "-Xmx8g")  # Max Java memory heap size, for rcausal

rm(list=ls())

data_dir = "~/DataProjects/cross-tissue"  # root of data directory
setwd("~/GoogleDrive/projects/STARNET/cross-tissue")

library(rcausal)
library(data.table)
library(igraph)
library(RColorBrewer)

source("src/parse.R")
source("src/base.R")

# Load cross-tissues modules
between = new.env()
load(file.path(data_dir, "modules/between_within-cross-tissue.RData"),
	between,
	verbose=TRUE)

between = parseModuleData(between)

# Load module overview table
module_tab = read.table("co-expression/tables/module_tab.csv",
	sep=",",
	header=TRUE
)

module_tab$class[module_tab$purity < 0.95] = "CT"  # cross-tissue
module_tab$class[module_tab$purity >= 0.95] = "TS"  # tissue-specific

pheno_pval = read.table("pheno/tables/pheno_pval.csv",
	sep=",",
	header=TRUE,
	check.names=FALSE
)

# Bootstrap Bayesian inference of eigengenes
# ------------------------------------
# Eigengenes
variables = cbind(
	between$bwnet$eigengenes
)

bootstrap_m = 1000
bn_bootstrap = lapply(1:bootstrap_m, function(k) {
	if (k %% 10 == 0) {
		message("iter ", k, " out of ", bootstrap_m)
	}
	# Sample with replacement from patient IDs
	idx = sample(nrow(variables), replace=TRUE)

	out = fges(variables[idx, ], maxDegree=100)

	return(out)
})

# save(bn_bootstrap, file="co-expression/eigenNetw/v3/bn_bootstrap_all.RData")
load("co-expression/eigenNetw/v3/bn_bootstrap_all.RData")

	
# Calculate average bootstrapped adjacency matrix
bn_bootstrap_adjmat = lapply(bn_bootstrap, function(bn) {
	# Character encoding of networks. Convert to adjacency matrix
	directed_edges = grep("-->", bn$edges)
	undirected_edges = grep("---", bn$edges)

	edge_list = do.call(rbind,
		strsplit(bn$edges[directed_edges], " --> ")
	)

	undirected_edge_list = do.call(rbind,
		strsplit(bn$edges[undirected_edges], " --- ")
	)

	# repeat inverse
	undirected_edge_list = rbind(undirected_edge_list, undirected_edge_list[, 2:1])

	# add to rest
	edge_list = rbind(edge_list, undirected_edge_list)

	edge_list = t(apply(edge_list, 1, as.integer))

	adj_mat = matrix(0, ncol(variables), ncol(variables))
	rownames(adj_mat) = colnames(variables)  # from
	colnames(adj_mat) = colnames(variables)  # to

	adj_mat[edge_list] = 1

	return(adj_mat)
})

boot_adjmat = Reduce("+", bn_bootstrap_adjmat) / length(bn_bootstrap_adjmat)
boot_p = 1 - boot_adjmat


# Assess significant supernetwork interactions
sig = boot_p < 0.2

snetw = which(sig, arr.ind=TRUE)

snetw = data.frame(snetw, check.rows=FALSE)
colnames(snetw) = c("from", "to")

snetw$p = boot_p[sig]

snetw = snetw[order(snetw$p), ]

# Add CT, TS status
snetw$from_class = module_tab$class[snetw$from]
snetw$to_class = module_tab$class[snetw$to]

snetw$from_to_class = paste(snetw$from_class, snetw$to_class, sep=" -> ")


# Plot of CT -> TS interactions
# ----------------------------------
# sum(snetw$from_to_class == "CT -> TS")

tissue_col = brewer.pal(9, "Set1")[-6][c(1, 5, 6, 4, 7, 3, 2)]
tissues = c("AOR", "MAM", "LIV", "VAF", "SF", "SKLM", "BLOOD")

# Fraction of tissues
tissue_counts = data.frame(module_tab)[, tissues]

tissue_frac = tissue_counts / apply(tissue_counts, 1, sum)
tissue_frac = as.matrix(tissue_frac)

module_tab$primary_tissue = colnames(tissue_frac)[apply(tissue_frac, 1, which.max)]

# Sort significant supernetwork edges by target tissue
snetw$to_tissue = factor(module_tab$primary_tissue[snetw$to], levels=tissues)

pdf("co-expression/eigenNetw/v3/plots/bootstrap_supernetw_by_tissue.pdf", height=6.0, width=4.0)
snetw = snetw[order(snetw$to_tissue, decreasing=TRUE), ]

idx = snetw$from_to_class == "CT -> TS"

# par(mfrow=c(2, 1), mar=c(3, 4, 2, 3))
par(mfrow=c(1, 2), mar=c(3, 4, 2, 3))
barplot(
	t(tissue_frac[snetw$from[idx], ]),
	main="From",
	col=tissue_col,
	names.arg=snetw$from[idx],
	las=2,
	legend.text=colnames(tissue_frac),
	space=0,
	horiz=TRUE,
	border=rgb(50, 50, 50, maxColorValue=255),
	cex.names=0.8
)

barplot(
	t(tissue_frac[snetw$to[idx], ]),
	main="To",
	col=tissue_col,
	names.arg=snetw$to[idx],
	las=2,
	legend.text=colnames(tissue_frac),
	space=0,
	border=rgb(50, 50, 50, maxColorValue=255),
	horiz=TRUE,
	cex.names=0.8
)
dev.off()



# Centrality of bootstrap eigennetwork
# ------------------------------------------
boot_p_bound = boot_p
boot_p_bound[boot_p_bound == 0] = 1/bootstrap_m

g = graph_from_adjacency_matrix(boot_p_bound,
	mode="directed",
	weighted=TRUE)

stopifnot(length(E(g)) == length(boot_p_bound))

# g = graph_from_adjacency_matrix((1 - boot_adjmat) < 0.2,
# 	mode="directed")
# 	# weighted=TRUE)

# # Eigenvector centrality
# centrality = eigen_centrality(g, directed=TRUE)

# stopifnot(all(V(g) == rownames(module_tab)))

# centrality_by_class = list(
# 	CT=centrality$vector[module_tab$class == "CT"],
# 	TS=centrality$vector[module_tab$class == "TS"]
# )

# boxplot(centrality_by_class, frame=FALSE)


# Weighted betweenness centrality
centrality = betweenness(g, directed=TRUE)

order(centrality, decreasing=TRUE)

centrality_by_class = list(
	TS=centrality[module_tab$class == "TS"],
	CT=centrality[module_tab$class == "CT"]
)


pdf("co-expression/eigenNetw/v3/plots/full_supernetwork_centrality.pdf", height=3.5, width=4.5)
par(mfrow=c(1, 2))
barplot(
	sort(table(snetw$from_to_class), decreasing=TRUE),
	las=2,
	ylab="Supernetwork interactions (P<0.2)",
	col=brewer.pal(8, "Paired")
)
abline(h=0)

colors = brewer.pal(9, "Set1")[c(3, 2)]

t_test = t.test(centrality_by_class$TS, centrality_by_class$CT)

boxplot(centrality_by_class,
	frame=FALSE,
	ylab="Network centrality (betweenness)",
	outline=FALSE,
	ylim=range(unlist(centrality_by_class)),
	main=paste0("P=", format(t_test$p.value, digits=3))
)

for (i in 1:length(centrality_by_class)) {
	points(
		jitter(rep(i, length(centrality_by_class[[i]])), amount=0.2),
		centrality_by_class[[i]],
		col=colors[i],
		pch=16,
		cex=0.5
	)
}
dev.off()
