rm(list=ls())

library(data.table)
library(biomaRt)
library(fmsb)
library(RColorBrewer)

setwd("~/Google Drive/projects/STARNET/cross-tissue")

source("src/base.R")


# Load biomaRt database for annotating the genes
ensembl = useEnsembl(biomart="ensembl",
	dataset="hsapiens_gene_ensembl", version=89)

gene_map = getBM(
	attributes=c("ensembl_gene_id", "hgnc_symbol", "gene_biotype"),
	mart=ensembl)


# Load co-expression modules
modules = fread("co-expression/tables/modules.csv")

modules$gene_biotype = gene_map$gene_biotype[match(modules$gene_symbol, gene_map$hgnc_symbol)]



# Size of module 78
sum(modules$clust == 78)  #module size
sum(modules$clust == 78 & modules$gene_biotype == "protein_coding", na.rm=TRUE)
table(modules$tissue[modules$clust == 78])
table(modules$tissue[modules$clust == 78 & modules$gene_biotype == "protein_coding"])

table(modules$gene_biotype[modules$clust == 78])


# Get tissue-specific signatures of modules 78 and 98
table(modules$tissue[modules$clust == 78])
modules[modules$clust == 78 & modules$tissue == "MAM", ]


table(modules$tissue[modules$clust == 98])
modules[modules$clust == 78 & modules$tissue == "MAM", ]


# Module 98
# -------------------------------------------------
write.table(
	modules$gene_symbol[modules$clust == 98 & modules$tissue == "LIV"],
	"co-expression/module78to98/signatures/mod98_LIV.txt",
	row.names=FALSE, col.names=FALSE, quote=FALSE
)




# All, module 78
# --------------------------------------
# MAM
write.table(
	modules$gene_symbol[modules$clust == 78 & modules$tissue == "MAM"],
	"co-expression/module78to98/signatures/mod78_MAM.txt",
	row.names=FALSE, col.names=FALSE, quote=FALSE
)

# LIV
write.table(
	modules$gene_symbol[modules$clust == 78 & modules$tissue == "LIV"],
	"co-expression/module78to98/signatures/mod78_LIV.txt",
	row.names=FALSE, col.names=FALSE, quote=FALSE
)

# SF
write.table(
	modules$gene_symbol[modules$clust == 78 & modules$tissue == "SF"],
	"co-expression/module78to98/signatures/mod78_SF.txt",
	row.names=FALSE, col.names=FALSE, quote=FALSE
)

# VAF
write.table(
	modules$gene_symbol[modules$clust == 78 & modules$tissue == "VAF"],
	"co-expression/module78to98/signatures/mod78_VAF.txt",
	row.names=FALSE, col.names=FALSE, quote=FALSE
)


# Protein-coding mRNA only
# --------------------------------
# LIV
write.table(
	na.omit(modules$gene_symbol[modules$clust == 78 &modules$gene_biotype == "protein_coding" & modules$tissue == "LIV"]),
	"co-expression/module78to98/signatures/mod78_protein-coding_LIV.txt",
	row.names=FALSE, col.names=FALSE, quote=FALSE
)

# SF
write.table(
	na.omit(modules$gene_symbol[modules$clust == 78 &modules$gene_biotype == "protein_coding" & modules$tissue == "SF"]),
	"co-expression/module78to98/signatures/mod78_protein-coding_SF.txt",
	row.names=FALSE, col.names=FALSE, quote=FALSE
)

# VAF
write.table(
	na.omit(modules$gene_symbol[modules$clust == 78 &modules$gene_biotype == "protein_coding" & modules$tissue == "VAF"]),

	"co-expression/module78to98/signatures/mod78_protein-coding_VAF.txt",
	row.names=FALSE, col.names=FALSE, quote=FALSE
)




# Load gene set enrichment tables
enrich = list()
enrich$mod78vaf = fread("co-expression/module78to98/signatures/enrichment/PANTHER_GO_mod78_VAF_all.txt", skip=10)
enrich$mod78sf = fread("co-expression/module78to98/signatures/enrichment/PANTHER_GO_mod78_SF_all.txt", skip=10)
enrich$mod78liv = fread("co-expression/module78to98/signatures/enrichment/PANTHER_GO_mod78_LIV_all.txt", skip=10)
enrich$mod98liv = fread("co-expression/module78to98/signatures/enrichment/PANTHER_GO_mod98_LIV_all.txt", skip=10)

# overrepresented only
enrich = lapply(enrich, function(d) {
	d[d[["upload_1 (fold Enrichment)"]] > 1, ]
})


go_terms = lapply(enrich, function(d) {
	d = d[d[["upload_1 (fold Enrichment)"]] > 2.0, ]
	terms = d[["GO biological process complete"]][1:5]
	return(terms)
})

go_terms = unique(unlist(go_terms))


fdr_mat = sapply(enrich, function(d) {
	d[["upload_1 (FDR)"]][match(go_terms, d[["GO biological process complete"]])]
})
rownames(fdr_mat) = sapply(strsplit(go_terms, "[(]"), function(x) x[1])

fdr_mat[fdr_mat > 1] = 1  # probabilities!!

min_pval = 1e-16

# cap pvalues
fdr_mat[fdr_mat < min_pval] = min_pval





colors = brewer.pal(9, "Set1")[c(4, 8, 7, 7)]

pdf("co-expression/module78to98/plots/GO_radar_89to98_2.pdf", width=12, height=4)
par(mfrow=c(1, 4))

for (i in 1:ncol(fdr_mat)) {

	data = data.frame(
		max=16,
		min=0,
		-log10(fdr_mat[, i]))
	data = t(data)


	# if (i == 1) {
	# 	axis_labels = colnames(data)
	# } else {
	# 	axis_labels = NA
	# }

	axis_labels = 1:ncol(data)


	radarchart(data.frame(data),
		# Quantitative axis label
		axistype=1,  # axis label, 1:cente
		axislabcol=rgb(100, 100, 100, maxColorValue=255),

		caxislabels=paste0("10e-", seq(0, 16, length.out=5)),
		calcex=0.8,  # axis tick marks font size

		# Axis names
		# vlcex=0.8,  # axis name font size
		vlcex=0.9,  # axis name font size
		vlabels=axis_labels,

		# Grid
		cglwd=0.7,  # grid line width
		cglcol=rgb(150, 150, 150, maxColorValue=255),  # grid color

		# Polygon
		pfcol=addAlpha(colors[i], 0.5),

		title=colnames(fdr_mat)[i],  # main
	)
}
dev.off()

pdf("co-expression/module78to98/plots/GO_radar_89to98_2_legend.pdf")
plot(0, 0, type="n")
legend("topleft", legend=paste0(1:nrow(fdr_mat), ". ", rownames(fdr_mat)), bty="n")
dev.off()


# Create data: note in High school for Jonathan:
data=as.data.frame(matrix( sample( 2:20 , 10 , replace=T) , ncol=10))
colnames(data)=c("math" , "english" , "biology" , "music" , "R-coding", "data-viz" , "french" , "physic", "statistic", "sport" )

# To use the fmsb package, I have to add 2 lines to the dataframe: the max and min of each topic to show on the plot!
data=rbind(rep(20,10) , rep(0,10) , data)

# The default radar chart proposed by the library:
radarchart(data)

