# Writes matrix of eigengenes to file
rm(list=ls())

data_dir = "~/DataProjects/cross-tissue"  # root of data directory
setwd("~/Google Drive/projects/STARNET/cross-tissue")

source("src/parse.R")


between = new.env()
load(file.path(data_dir, "modules/between_within-cross-tissue.RData"),
	between,
	verbose=TRUE)


between = parseModuleData(between)

bwnet$eigengenes

eigengene_mat = between$bwnet$eigengenes
rownames(eigengene_mat) = between$patient_ids

write.csv(eigengene_mat,
	"co-expression/tables/eigengene_mat.csv",
	quote=FALSE
)