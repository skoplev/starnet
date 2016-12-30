i = 17

plot(con_eval_pairs[[i]]$fitIndices$Power,
	con_eval_pairs[[i]]$fitIndices$SFT.R.sq,
	main=names(con_eval_pairs)[i],
	type="l"
)


# i = 2
par(mfrow=c(3, 3))

for (i in 1:length(thresh_eval)) {

	tissue = names(thresh_eval)[i]
	plot(thresh_eval[[i]]$fitIndices$Power,
		thresh_eval[[i]]$fitIndices$SFT.R.sq,
		main=names(thresh_eval)[i],
		type="l",
		lwd=2.0,
		ylim=c(0, 1)
	)

	# Find matching

	pairs_idx = which(apply(paired_tissue, 2, function(col) tissue %in% col))
	for (idx in pairs_idx) {
		lines(thresh_eval_pairs[[idx]]$fitIndices$Power,
			thresh_eval_pairs[[idx]]$fitIndices$SFT.R.sq,
			main=names(thresh_eval_pairs)[idx],
			type="l", col="grey"
		)
	}

}


sapply(thresh_eval_pairs, function(x) x$powerEstimate)

# Average best slope
deg_slopes = sapply(thresh_eval, function(x) x$fitIndices$slope)
mean_deg_slopes = apply(deg_slopes, 1, mean)

best_pow = which.min(abs(1 + mean_deg_slopes))
message("Best average beta for tissues: ", powers[best_pow])


dir.create("plots")

svg("plots/scale-free.svg", width=3.7, height=4)
colors = brewer.pal(9, "Set1")[c(1:5, 7:9)]
plot(powers, thresh_eval[[1]]$fitIndices$slope, type="n",
	main="Scale-free co-expression networks",
	ylim=c(-0.5, -1.5),
	xlab=expression(beta),
	ylab="log-log slope"
)

abline(h=-1.0, col="grey", lty=3, lwd=1.5)
abline(v=powers[best_pow], col="grey", lwd=1.5)
for (i in 1:length(thresh_eval)) {
	lines(powers, thresh_eval[[i]]$fitIndices$slope, col=colors[i], lwd=1.2)
	points(powers, thresh_eval[[i]]$fitIndices$slope, pch=16, col=colors[i], cex=0.4)
}
legend("bottomright", legend=names(thresh_eval), col=colors, pch=16, cex=0.5)
dev.off()

plot(powers, apply(deg_slopes, 1, mean))

# par(mfrow())
# lapply(thresh_eval, function(sft) {
# 	plot(powers, -sft$fitIndices$slope)
# 	abline(h=1)
# })

# # cex1=0.9
# plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
# 	xlab="Soft Threshold (power)",
# 	ylab="Scale Free Topology Model Fit,signed R^2",
# 	# type="n",
# 	main = paste("Scale independence"),
# 	pch=16
# )

# text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
# labels=powers,cex=cex1,col="red");