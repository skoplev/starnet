# Impute missing tissue data
rm(list=ls())

library(impute)

data_dir = "/Users/sk/DataProjects/cross-tissue"

# Loads expr_recast data frame with tissue-specific expression
load("/Users/sk/DataProjects/cross-tissue/STARNET/gene_exp_norm_reshape/expr_recast.RData")

# Get expression matrix from recasted data frame with all tissue 
mat = expr_recast[, 3:ncol(expr_recast)]
mat = data.matrix(mat)
row_meta = expr_recast[, 1:2]
colnames(mat) = colnames(expr_recast)[3:ncol(expr_recast)]

# Remove samples with more than 50% missing data
missing_frac = apply(mat, 2, function(col) {
		sum(is.na(col)) / length(col)
})
mat = mat[, missing_frac < 0.5]

# Missing data message, remining samples
message("missing data fraction: ", sum(is.na(mat))/(nrow(mat) * ncol(mat)))

# Impute missing data using k-nearest neighbors
mat = impute.knn(mat, k=20, colmax=0.5, maxp=10000)$data

save(mat, row_meta, file=file.path(data_dir, "STARNET/gene_exp_norm_batch_imp/all.RData"))
