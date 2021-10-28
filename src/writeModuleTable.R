rm(list=ls())

data_dir = "/Users/sk/DataProjects/cross-tissue"  # root of data directory
setwd("/Users/sk/Google Drive/projects/STARNET/cross-tissue")
source("src/parse.R")

between = new.env()
load(file.path(data_dir, "modules/between_within-cross-tissue.RData"),
	between,
	verbose=TRUE)

# Parse module data
between = parseModuleData(between)


mod_tab = between$meta_genes

mod_tab$clust = between$clust
mod_tab$color = between$bwnet$colors

write.csv(mod_tab, "co-expression/tables/modules.csv", row.names=FALSE)
