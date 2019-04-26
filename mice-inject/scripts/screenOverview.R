rm(list=ls())

setwd("~/GoogleDrive/projects/STARNET/cross-tissue")

library(RColorBrewer)

n = 16  # total candidates
k = 10  # 
m = 4

candidates = c(m, k-m, n-k)

pdf("mice-inject/plots/validation_rate.pdf", width=1.7, height=4)
barplot(as.matrix(candidates), 
	# besides=TRUE,
	col=rev(brewer.pal(3, "Greens")),
	ylab="Endocrine candidates"
)
dev.off()