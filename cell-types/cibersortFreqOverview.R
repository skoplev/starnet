library(data.table)
library(gplots)

data_dir = "/Users/sk/DataProjects/cross-tissue"  # root of data directory
setwd("/Users/sk/Google Drive/projects/cross-tissue")

source("co-expression/base.R")
source("parse/io.R")

# Load CIBERSORT frequency data
freq = loadCibersortFreq(file.path(data_dir, "CIBERSORT/out_freq"))

all_freq = do.call(rbind, freq)


# Get tissue of combined frequency matrix
tissue = sapply(
	strsplit(all_freq$tissue_patient_id, "_"),
	function(x) x[1]
)
tissue = factor(tissue)


# Heatmap of all tissue samples
# -----------------------------------------------------
colors = c(brewer.pal(9, "Set1"), brewer.pal(8, "Set2"))
sample_col = colors[as.integer(tissue)]

# sd_thresh = 0.005
sd_thresh = 0.01
# sd_thresh = 0.015
# sd_thresh = 0.02

mat = all_freq[, 2:(ncol(all_freq) - 3), with=FALSE]
mat = as.matrix(mat)


mat = mat[, apply(mat, 2, sd) > sd_thresh]
mat = scale(mat)

mat = t(mat)

pdf("cell-types/plots/cell_type_freq_heatmap2.pdf")
heatmap.2(mat,
	trace="none",

	# col=colorRampPalette(brewer.pal(9, "Greens"))(100),
	# breaks=seq(0, 0.15, length.out=101),  # cap of coloring 

	col=colorRampPalette(rev(brewer.pal(9, "Spectral")))(100),
	breaks=seq(-4, 4, length.out=101),  # cap of coloring 

	ColSideColors=sample_col,
	# cexRow=0.4,
	# cexRow=0.7,
	margins=c(8, 18),
	xlab="Tissue samples",
	labCol="",
	key.xlab="Frequency z-score",
	key.title=NA
)
legend("topright", col=colors, pch=15, legend=levels(tissue), cex=0.5)
dev.off()
