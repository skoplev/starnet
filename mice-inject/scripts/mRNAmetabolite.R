rm(list=ls())

library(data.table)
library(RColorBrewer)
library(reshape2)
library(WGCNA)
library(corrplot)

setwd("~/GoogleDrive/projects/STARNET/cross-tissue")

# Load gene expression
# ---------------------------------------
count_norm = fread("mice-inject/data/count_norm_batch.csv")

# gene_symbols = make.names(count_norm$V1, unique=TRUE)
gene_symbols = count_norm$V1

count_norm = data.matrix(count_norm[, -1])

rownames(count_norm) = gene_symbols

# Load RNA expression metadata
samples = fread("mice-inject/data/samples_meta.csv")
samples$condition = factor(samples$condition)

stopifnot(colnames(count_norm) == samples$id)

# Load metabolite data
# -----------------------------

metabolite_dir = "mice-inject/data/metabolite"
metabolite_files = list.files(metabolite_dir)

metabolite = lapply(metabolite_files, function(file_name) {
	fread(file.path(metabolite_dir, file_name))
})
names(metabolite) = sapply(strsplit(metabolite_files, "[.]"), function(x) x[1])


# (Re) construct sample IDs

# Day 1
metabolite$LBP$sample_id = paste0(metabolite$LBP$Protein, "_1")

# Day 2
metabolite$EPDR1$sample_id = paste0(metabolite$EPDR1$Protein, "_2")
metabolite$FCN2$sample_id = paste0(metabolite$FCN2$Protein, "_2")
metabolite$FSTL3$sample_id = paste0(metabolite$FSTL3$Protein, "_2")

# Append .n suffix
metabolite = lapply(metabolite, function(tab) {
	tab$exp_group = paste(tab$Protein, tab$Metabolite)
	# add incremental .n suffix to IDs per group
	for (group in unique(tab$exp_group)) {
		idx = tab$exp_group == group

		prop_ids = tab$sample_id[idx]

		prop_ids = paste0(prop_ids, ".", 1:length(prop_ids))

		tab$sample_id[idx] = prop_ids
	}
	return(tab)
})

# Exclude FSTL Veh (these might have beend added erroneously)
metabolite$FSTL3 = metabolite$FSTL3[metabolite$FSTL3$Protein != "Veh", ]


# Recast metabolite data
# ------------------------------------------
# dcast(metabolite[[1]], sample_id~Metabolite, value.var="Value")
# dcast(metabolite[[2]], sample_id~Metabolite, value.var="Value")

metab_mat = lapply(metabolite, function(tab) {
	dcast(tab, sample_id~Metabolite, value.var="Value")
})

metab_mat = rbindlist(metab_mat)

sample_ids = metab_mat$sample_id
metab_mat = data.matrix(metab_mat[, -1])
rownames(metab_mat) = sample_ids

# Match to gene expression
metab_mat_match = metab_mat[match(samples$id, rownames(metab_mat)), ]





genes = c(
	"Fdft1",
	"Nsdhl",
	"Tdrkh",
	"Dhcr7",
	"Acat3",
	"Fam213a",
	"Rdh11",

	"Cyp51",
	"Sqle",
	"Pgd",
	"Msmo1",
	"Hmgcs1",

	"Pcsk9",
	"Lss",

	"Acaca"
	# "Dhcr24"
)

idx = match(genes, rownames(count_norm))

cor_tests = corAndPvalue(x=t(count_norm[idx, ]), y=metab_mat_match[, -1])


sum(cor_tests$p < 0.05)

t(cor_tests$p)

pdf("mice-inject/plots/mRNA_metabolite_corrplot.pdf", width=6, height=3)
# pdf("mice-inject/plots/mRNA_metabolite_corrplot.pdf", width=3, height=6)
corrplot(
	t(cor_tests$cor),
	p.mat=t(cor_tests$p),
	# cor_tests$cor,
	# p.mat=cor_tests$p,
	method="ellipse",
	col=colorRampPalette(rev(brewer.pal(9, "RdBu")))(200),
	outline=TRUE,
	addgrid.col=NA,
	tl.col="black",
	insig="label_sig",
	sig.level=c(0.001, 0.01, 0.05),
	pch.cex=1.0
)
dev.off()

# heatmap.2(cor_tests$cor,
# 	mar=c(20, 10),
# 	trace="none",
# 	col=colorRampPalette(rev(brewer.pal(9, "RdBu")))
# )



scatterPlot = function(x, y, ...) {
	cor_test = cor.test(x, y)

	plot(x, y,
		main=paste0(
			"P=", format(cor_test$p.value, digits=3),
			", r=", format(cor_test$estimate, digits=3)
		),
		bty="n",
		pch=21,
		bg=pts_col,
		...
	)
}

colors = c(brewer.pal(9, "Set1")[c(1, 2, 4, 3)], "white")
# levels(samples$condition)
pts_col = colors[as.integer(factor(samples$condition))]


# i = which(rownames(count_norm) == "Pcsk9")

# i = which(rownames(count_norm) == "Dhcr7")

i = which(rownames(count_norm) == "Dhcr24")
j = 4  # Plasma TC

x = count_norm[i, ]
y = metab_mat_match[, j]

scatterPlot(x, y,
	xlab=paste0(rownames(count_norm)[i], "(adj. counts)"),
	ylab=ylab
)



pdf("mice-inject/plots/scatter_mRNA_metabolites.pdf", height=3.5, width=6)
j = 5  # liver TC

par(mfrow=c(1, 2))
i = which(rownames(count_norm) == "Rdh11")
ylab = gsub("X 100", "", colnames(metab_mat_match)[j])

x = count_norm[i, ]
y = metab_mat_match[, j] / 100

scatterPlot(x, y,
	xlab=paste0(rownames(count_norm)[i], "(adj. counts)"),
	ylab=ylab
)


i = which(rownames(count_norm) == "Hmgcs1")

x = count_norm[i, ]
y = metab_mat_match[, j] / 100

scatterPlot(x, y,
	xlab=paste0(rownames(count_norm)[i], "(adj. counts)"),
	ylab=ylab
)

legend("topright",
	legend=levels(samples$condition),
	pch=21,
	pt.bg=colors)
dev.off()
