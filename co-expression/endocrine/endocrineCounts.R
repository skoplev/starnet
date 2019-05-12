# Counts and tables of endocrine factors identified
#
rm(list=ls())

library(data.table)
library(WGCNA)
library(RColorBrewer)
library(gplots)

setwd("~/GoogleDrive/projects/STARNET/cross-tissue")

source("src/permuteTest.R")
source("src/parse.R")

data_dir = "/Users/sk/DataProjects/cross-tissue"  # root of data directory

endo = fread("co-expression/tables/CT_endocrines_TS_interactions.csv")


source_tissue_counts = sort(table(endo$tissue), decreasing=TRUE)

target_tissue_counts = table(endo$target_tissue_primary)

# Order by number of source endocrine candidates
target_tissue_counts = target_tissue_counts[match(names(source_tissue_counts), names(target_tissue_counts))]

pdf("co-expression/endocrine/plots/endocrine_candidates_793_margins.pdf", width=5, height=6)
par(mfrow=c(2, 1))
barplot(source_tissue_counts, las=3, main="Source")
barplot(target_tissue_counts, las=2, main="Target")
dev.off()

tissues = names(source_tissue_counts)

pdf("co-expression/endocrine/plots/endocrine_candidates_793_count_matrix.pdf", width=6, height=6)
endo_count = matrix(NA, ncol=length(tissues), nrow=length(tissues))
for (i in 1:length(tissues)) {
	for (j in 1:length(tissues)) {
		endo_count[i, j] = sum(endo$tissue == tissues[i] & endo$target_tissue_primary == tissues[j])
	}
}
colnames(endo_count) = tissues
rownames(endo_count) = tissues


heatmap.2(endo_count,
	Rowv=FALSE,
	Colv=FALSE,
	trace="none",
	col=colorRampPalette(c("white", brewer.pal(9, "Greens")))(100),
	cellnote=endo_count,
	notecol="black"
)
dev.off()


# endocrine factors from module 78
# from VAF
sum(endo$clust == 78 & endo$tissue == "VAF")
# from SF
sum(endo$clust == 78 & endo$tissue == "SF")

# targeting module 98
sum(endo$clust == 78 & endo$tissue == "VAF" & endo$target_clust == 98) 
sum(endo$clust == 78 & endo$tissue == "SF" & endo$target_clust == 98)

