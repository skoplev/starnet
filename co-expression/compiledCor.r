# experimental optimization of correlation calculations.
library(Rcpp)
library(Matrix)
sourceCpp("partCor/cor.cpp")

m = 20000  # features to test on
mat = expr_recast[sample(m), 3:ncol(expr_recast)]

# mat = expr_recast[, 3:ncol(expr_recast)]
# rownames(mat) = expr_recast$transcript_id
# rownames(mat) = paste(expr_recast$tissue, expr_recast$transcript_id, sep=";")

colnames(mat) = colnames(expr_recast)[3:ncol(expr_recast)]
mat = as.matrix(mat)
# mat[is.na(mat)] = 0
mat = t(mat)  # features in columns, obs in rows

min_cor = 0.3

# Chunk features for correlation calculations
chunk_size = 5000  # number of features in each chunk
chunks = split(1:ncol(mat), ceiling(1:ncol(mat) / chunk_size))  # feature chunks

# Combinations of feature chunks to be evaluated by the correlation
chunk_idx = c(
	# all diagonal chunk blocks
	lapply(1:length(chunks), function(n) c(n, n)),
	# lower triangular blocks, pairwise combinations, without order.
	combn(1:length(chunks), 2, simplify=FALSE)
)

flat_cors = mclapply(chunk_idx, function(mat_block) {
	message("Chunk: ", mat_block)
	i = mat_block[1]
	j = mat_block[2]

	if (i == j) {
		# Within-chunk correlation
		cmat = cor(
			mat[, chunks[[i]]],
			use="pairwise.complete.obs", nThreads=2)
	} else {
		# Cross-chunk correlation
		cmat = cor(
			mat[, chunks[[i]]],
			mat[, chunks[[j]]],
			use="pairwise.complete.obs", nThreads=2)
	}

	# Filter correlations below threshold, creating flat data stucture
	include = abs(cmat) > min_cor
	from_to = which(include, arr.ind=TRUE)

	cor_tab = cbind(
		which(include, arr.ind=TRUE),
		cmat[include]
	)
	colnames(cor_tab) = c("row", "col", "cor")

	# translate from chunk to global feature coordinates
	cor_tab[,1] = chunks[[i]][cor_tab[,1]]
	cor_tab[,2] = chunks[[j]][cor_tab[,2]]

	# Remove self edges
	cor_tab = cor_tab[cor_tab[,1] != cor_tab[,2],]

	return(cor_tab)
}, mc.cores=3)

flat_cors = do.call("rbind", flat_cors)
mat_ids = paste(flat_cors[,1], flat_cors[,2], sep="_")

# Forse triangular
swap = flat_cors[,1] > flat_cors[,2]
tmp = flat_cors[swap, 1]
flat_cors[swap, 1] = flat_cors[swap, 2]
flat_cors[swap, 2] = tmp

# flat_cors[,1] == flat_cors[,2] & flat_cors[,2] == flat_cors[,1]
sum(flat_cors[,1] == flat_cors[,2] & flat_cors[,2] == flat_cors[,1])
# sum(flat_cors[,1] == flat_cors[,2] & flat_cors[,2] == flat_cors[,1])

cor_mat = sparseMatrix(i=flat_cors[,1], j=flat_cors[,2], x=flat_cors[,3], symmetric=TRUE)
# max(cor_mat)
cor_mat = sparseMatrix(i=flat_cors[,1], j=flat_cors[,2], x=flat_cors[,3])
adj_mat = abs(cor_mat)^3
plot(density(as.vector(adj_mat)))

tom =

# cor_mat = sparseMatrix(i=flat_cors[,1], j=flat_cors[,2], x=flat_cors[,3])
# cor_mat = forceSymmetric(cor_mat)

# list(flat_cors=flat_cors, nodes=)




# part_cors = partCor(t(mat[1:n,]), 0.3)

# all_cors = cor(t(mat[idx,]), nThreads=4)

# all_cors[abs(all_cors) < min_cor] = 0.0

# include = abs(all_cors) > min_cor
# from_to = which(include, arr.ind=TRUE)
# weights = all_cors[include]

# network = cbind(from_to, weights)

sparse_cors = Matrix(all_cors, sparse=TRUE)

all_cors2 = cor(t(mat[idx + 100,]))
all_cors2[abs(all_cors2) < min_cor] = 0.0
sparse_cors2 = Matrix(all_cors2, sparse=TRUE)

merge(sparse_cors, sparse_cors2, by="row.names", all=TRUE)

comb_cors = all_cors + all_cors2

Matrix(all_cors, all_cors2)
