options(java.parameters = "-Xmx8g")  # Max Java memory heap size, for rcausal

rm(list=ls())

library(rcausal)
library(data.table)
library(igraph)

data_dir = "/Users/sk/DataProjects/cross-tissue"  # root of data directory
setwd("/Users/sk/Google Drive/projects/cross-tissue")

source("src/parse.R")


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

variables = cbind(
	between$bwnet$eigengenes
)

bn = fges(variables, maxDegree=100)  # Fast-greedy equivalence search


g = graph_from_graphnel(bn$graphNEL)

epsilon = 10^-30

# Add featueres
# V(g)$log_P_cad = -log10(module_tab$cad_pval)
V(g)$mod_size_square = sqrt(module_tab$mod_size) + epsilon

V(g)$mod_size_square = sqrt(module_tab$mod_size) + epsilon

V(g)$frac_AOR = module_tab$AOR / module_tab$mod_size
V(g)$frac_BLOOD = module_tab$BLOOD / module_tab$mod_size
V(g)$frac_LIV = module_tab$LIV / module_tab$mod_size
V(g)$frac_MAM = module_tab$MAM / module_tab$mod_size
V(g)$frac_SKLM = module_tab$SKLM / module_tab$mod_size
V(g)$frac_SF = module_tab$SF / module_tab$mod_size
V(g)$frac_VAF = module_tab$VAF / module_tab$mod_size

V(g)$log_pval_CAD = -log10(module_tab$pval)  #  CAD enrichment
V(g)$log_pval_DUKE = -log10(module_tab$pval_DUKE)  #  DUKE correlation

V(g)$log_n_within_mod_endocrines = log10(module_tab$n_within_mod_endocrines)
V(g)$log_n_within_mod_endocrines[is.infinite(V(g)$log_n_within_mod_endocrines)] = 0

V(g)$log_n_key_driver_eQTL = log10(module_tab$n_key_driver_eQTL_genes)
V(g)$log_n_key_driver_eQTL[is.infinite(V(g)$log_n_key_driver_eQTL)] = 0

write_graph(g,
	file="co-expression/eigenNetw/v2/bayes_net.gml",
	format="gml")


