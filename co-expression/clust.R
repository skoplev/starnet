# Clustering of combined tissue-specific expression matrices.

library(cluster)
# K-means
# n = 50000  # subset size
# k = 16
# dmat = mat[sample(n),]
# dmat = mat


# Standardize matrix data, overwrite for memory efficiency
mat = scale(t(mat))  # obs in rows
mat = abs(mat)  # absoute z-scores, to allow clustering of inversly correlated features
mat[is.na(mat)] = 0.0  # imputation to average
mat = t(mat)  # features in rows

# Run k-means in parallel for a range of k
# scaled_mat_sub = mat[,1:1000]
clust = mclapply(1:25, function(k) {
	kmeans(mat, k)
}, mc.cores=3)

# clust = kmeans(t(scaled_mat), k)
# table(clust$cluster)


par(
	mfrow=c(3, 1),
	mar=c(2, 4.2, 2, 2),
	oma=c(2, 0, 0, 0)
)
plot(
	sapply(clust, function(x) x$tot.withinss),
	type="l", lwd=1.5,
	ylab="Within sum of squares"
)

plot(
	sapply(clust, function(x) x$betweenss),
	type="l", lwd=1.5,
	ylab="Between sum of squares"
)

plot(
	sapply(clust, function(x) (max(table(x$cluster)))),
	type="l", lwd=1.5,
	ylab="Max cluster size"
)

k = 10
tissue_clust = sapply(1:k, function(i) {
	table(expr_recast$tissue[clust[[k]]$cluster == i])
})

tissue_col = brewer.pal(12, "Set3")

counts = dcast(melt(tissue_clust), L1~Var1, value.var="value")
counts = counts[,2:ncol(counts)]
counts[is.na(counts)] = 0
counts = counts[order(apply(counts, 1, sum), decreasing=TRUE),]  # reorder
counts = as.matrix(counts)
counts = t(counts)

barplot(counts, col=tissue_col, names.arg=rep("", k), xlab="Cluster")
legend("topright", legend=rownames(counts), pch=15, col=tissue_col)


# Plots for k-means clustering
# PCA of clusters
pca = prcomp(t(scaled_mat))


pts_col = rep("black", n)
# for (i in 1:length(clust)) {
# 	pts_col[clust[[i]]] = brewer.pal(9, "Set1")[i]
# }
pts_col = brewer.pal(9, "Set1")[clust$cluster]

plot(pca$x[,1], pca$x[,2], col=pts_col)
