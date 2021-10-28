rm(list=ls())

library(data.table)
library(RColorBrewer)

setwd("~/GoogleDrive/projects/STARNET/cross-tissue")



mod_tab = fread("co-expression/tables/module_tab.csv")



pdf("co-expression/determinePower/plots/purity_histogram.pdf", height=4, width=7)
hist((1 - mod_tab$purity) * 100,
	breaks=100,
	xlab="Percent cross-tissue transcripts, per module",
	main="",
	col=brewer.pal(9, "Set1")[2])
abline(v=5,
	# col="grey",
	col=brewer.pal(9, "Set1")[1],
	lty=2)
dev.off()

# x = seq(from=0, to=1-min(mod_tab$purity), length.out=200)

# frac_CT = sapply(x, function(cutoff) {
# 	mean(mod_tab$purity <= (1 - cutoff))
# })

# plot(x * 100, frac_CT, type="l", bty="l")
# abline(v=5, col="grey", lty=2)
