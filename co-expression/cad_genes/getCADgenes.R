rm(list=ls())

data_dir = "/Users/sk/DataProjects/cross-tissue"  # root of data directory
setwd("~/GoogleDrive/projects/STARNET/cross-tissue")

source("src/parse.R")

cad_genes = getCADGenes(data_dir)

# Nelson eQTL mapping, FDR <0.05
nelson_eqtl = fread("co-expression/cad_genes/Nelson_variants/nelson_eQTL_networks.csv")
nelson_genes = unique(nelson_eqtl$hgnc_symbol)

write(cad_genes, file="co-expression/cad_genes/cad_genes.txt")

length(cad_genes)
length(unique(nelson_genes, cad_genes))

write(unique(nelson_genes, cad_genes), file="co-expression/cad_genes/cad_genes_nelson_eQTL.txt")
