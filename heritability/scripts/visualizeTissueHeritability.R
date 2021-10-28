rm(list=ls())

library(data.table)
library(RColorBrewer)

setwd("~/GoogleDrive/projects/STARNET/cross-tissue")


# Load tissue heritability estimates
heri = fread("heritability/from_johan/tissue/Tissue_H2 (3)20200107.csv")

# Remove empty rows
heri = heri[heri[["Module ID"]] != "", ]

# Comparative plots of heritability estimates
# ------------------------------------

# Extract features to render
features = c("CombinedModules2k", "CombinedModules1k", "CombinedModules500", "All tissue-specific", "All cross-tissue")

h2_CAD = heri$pct_h2_CAD[match(features, heri[["Module ID"]])]
names(h2_CAD) = features


# Prepend
h2_CAD = c(
	c("Lead SNPs"=21.2, "eQTL SNPs"=59.79),  # Nelson_2017
	h2_CAD
)



pdf("heritability/plots/heritability_overview.pdf", width=3, height=4)
par(mar=c(10, 4, 2, 2))

bar_col = c(
	brewer.pal(9, "Set1")[9],
	rep(brewer.pal(9, "Set1")[1], 4),
	rep(brewer.pal(9, "Set1")[2], 2)
)

barplot(h2_CAD,
	ylab="H2 CAD (%)",
	border=NA,
	# space=0.66,
	space=1.33,
	col=bar_col,
	las=2)

abline(h=h2_CAD[1]
	# col=rgb(100, 100, 100, maxColorValue=255),
	# lty=2
)

abline(h=0)
dev.off()


h2_CAD_stacked = rbind(
	rep(h2_CAD[1], length(h2_CAD) - 1),
	h2_CAD[-1]
)
colnames(h2_CAD_stacked) = names(h2_CAD)[-1]

pdf("heritability/plots/heritability_overview_stacked.pdf", width=2.8, height=4)

bar_col=c("white", brewer.pal(8, "Pastel1")[1])

par(mar=c(10, 4, 2, 2))
barplot(h2_CAD_stacked,
	ylab="H2 CAD (%)",
	# ylim=c(0, 100),
	las=2,
	col=bar_col
)
abline(h=0)

dev.off()


# Pie charts
# ---------------------------
h2_CAD_pie = rbind(
	h2_CAD_stacked,
	100 - apply(h2_CAD_stacked, 2, sum)
)
rownames(h2_CAD_pie) = c("Lead SNPs", "STARNET eQTLs", "Epistasis")



pdf("heritability/plots/heritability_pies.pdf")
par(mfrow=c(2, 3))

for (i in 1:ncol(h2_CAD_pie)) {
	pie(h2_CAD_pie[, i], main=colnames(h2_CAD_pie)[i],
		col=brewer.pal(9, "Set1")
		)
}
dev.off()


# Module size exclusion tiration of heritability
pdf("heritability/plots/heritability_module_titration.pdf", width=2.6, height=4)
par(mar=c(10, 4, 2, 2))
barplot(h2_CAD[2:5], 
	ylab="CAD H2 (%)",
	las=2,
	# col=rev(brewer.pal(6, "Greens"))[-1],
	col=rev(brewer.pal(7, "Blues")),
	border=NA)
abline(h=0)
dev.off()


pdf("heritability/plots/heritability_CT_TS_barplot.pdf", width=1.8, height=4)
par(mar=c(10, 4, 2, 2))
barplot(
	h2_CAD[6:7],
	ylab="CAD H2 (%)",
	las=2,
	col=brewer.pal(9, "Set1")[2],
	border=NA
)
abline(h=0)
dev.off()


all_idx = 20:26
ct_idx = 13:19
ts_idx = 1:7

heri[["Module ID"]][all_idx]
heri[["Module ID"]][ct_idx]
heri[["Module ID"]][ts_idx]


h2_by_tissue = rbind(
	heri$pct_h2_CAD[all_idx],
	heri$pct_h2_CAD[ts_idx],
	heri$pct_h2_CAD[ct_idx]
)

rownames(h2_by_tissue) = c("All", "Cross-tissue", "Tissue-specific")

colnames(h2_by_tissue) = sapply(strsplit(heri[["Module ID"]][all_idx], " "), function(x) x[2])


pdf("heritability/plots/heritability_by_tissue_barplot.pdf", height=4)
barplot(h2_by_tissue,
	beside=TRUE,
	# las=2,
	ylab="CAD H2 (%)",
	density=c(0, 50, 30),
	angle=c(0, 0, 45)
	# col=brewer.pal(9, "Set1")
)
dev.off()



# par(mfrow=c(3, 1))

# for (i in 1:3) {
# 	barplot(h2_by_tissue[i, ],
# 		main=rownames(h2_by_tissue)[i],
# 		col=brewer.pal(9, "Set1"),
# 		border=NA
# 	)
# 	abline(h=0)
# }


pdf("heritability/plots/heritability_by_tissue_dotplot.pdf", width=3, height=4)
tissue_colors = brewer.pal(9, "Set1")[c(1, 5, 7, 3, 4, 8, 2)]

plot(0, 0, type="n",
	xlim=c(0.5, 3.5),
	ylim=c(0, max(h2_by_tissue)),
	bty="l",
	ylab="CAD H2 (%)",
	xlab="All, TS, CT"
)


for (j in 1:ncol(h2_by_tissue)) {
	lines(1:3, h2_by_tissue[, j],
		lwd=1.5,
		col=tissue_colors[j])
}

for (i in 1:3) {
	points(
		rep(i, ncol(h2_by_tissue)),
		h2_by_tissue[i, ],
		pch=21,
		cex=1.2,
		bg=tissue_colors)
}

legend("topright",
	legend=colnames(h2_by_tissue),
	pt.bg=tissue_colors,
	pch=21
)
dev.off()