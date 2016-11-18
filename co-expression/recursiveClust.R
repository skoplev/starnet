# Experimental recursive clustering using CLARA;
# is a lot slower than k-means.

max_size = 40000


# Returns row idx for input matrix.
recursiveClust = function(dmat, max_size) {
	# Base case: cluster is small enough
	if (nrow(dmat) < max_size) {
		return(
			list(1:nrow(dmat))
		)
	}

	# Cluster rows into 2 clusters
	# clust = clara(dmat, 2, samples=50, correct.d=TRUE)
	clust = pam(dmat, 2)

	idx_clust1 = which(clust$cluster == 1)
	idx_clust2 = which(clust$cluster == 2)

	names(idx_clust1) = NULL
	names(idx_clust2) = NULL

	# Subdivide
	idx_sub1 = recursiveClust(dmat[idx_clust1, ], max_size)
	idx_sub2 = recursiveClust(dmat[idx_clust2, ], max_size)

	# Translate row ids from submatrices to input matrix
	idx1 = lapply(idx_sub1, function(idx) idx_clust1[idx])
	idx2 = lapply(idx_sub2, function(idx) idx_clust2[idx])

	return(c(idx1, idx2))
}


# Recursive CLARA clustering
n = 5000
max_size= 100
dmat = mat[1:n,]

clust = recursiveClust(dmat, max_size)
