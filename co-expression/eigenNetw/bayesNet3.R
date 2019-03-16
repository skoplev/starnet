#
options(java.parameters = "-Xmx8g")  # Max Java memory heap size, for rcausal

rm(list=ls())

data_dir = "~/DataProjects/cross-tissue"  # root of data directory
setwd("~/GoogleDrive/projects/STARNET/cross-tissue")

library(rcausal)
library(RColorBrewer)
library(gplots)
library(igraph)
library(reshape)

source("src/parse.R")
source("src/base.R")

# Load cross-tissues modules
between = new.env()
load(file.path(data_dir, "modules/between_within-cross-tissue.RData"),
	between,
	verbose=TRUE)

between = parseModuleData(between)

# Load module overview table
module_tab = read.table("co-expression/tables/module_tab.csv",
	sep=",",
	header=TRUE
)

pheno_pval = read.table("pheno/tables/pheno_pval.csv",
	sep=",",
	header=TRUE,
	check.names=FALSE
)


# CAD modules criteria
# ---------------------------------------------------------------

# pheno_padj = data.frame(
# 	SYNTAX=p.adjust(pheno_pval$syntax_score, method="BH"),
# 	DUKE=p.adjust(pheno_pval$DUKE, method="BH"),
# 	Case_Ctrl_DEG=p.adjust(pheno_pval$case_control_DEG, method="BH")
# )

features = c("syntax_score", "DUKE", "case_control_DEG")
pheno_padj = matrix(
	p.adjust(
		data.matrix(pheno_pval[, features]),
		method="BH"),
	ncol=3)
colnames(pheno_padj) = c("SYNTAX", "DUKE", "Case-Ctrl DEG")

fdr = 0.01



# Bootstrapped network inferrence over eigengenes
# ----------------------------------------------------

variables = cbind(
	between$bwnet$eigengenes[, cad_modules]
)

mat = data.matrix(variables)

sd_cutoff = quantile(apply(mat, 1, sd), probs=0.75)[1]


# Plot standard deviation distribution and cutoff
pdf("co-expression/eigenNetw/v3/plots/CAD_eigengene_variability.pdf", height=4, width=6)
patient_sd = apply(mat, 1, sd)
patient_sd_density = density(patient_sd, from=0)
plot(patient_sd_density,
	xlab=expression("Std. dev. (" * sigma * ")"),
	main="CAD eigengene variability",
	zero.line=FALSE,
	bty="n"
)
polygon(
	c(
		0,
		patient_sd_density$x,
		max(patient_sd_density$x)
	),
	c(
		0,
		patient_sd_density$y,
		0
	),
	col=brewer.pal(9, "Pastel1")[2]
)

abline(v=sd_cutoff,
	lwd=1.5,
	col=brewer.pal(9, "Set1")[1])

legend("topright", legend=c("75th percentile"),
	pch="|",
	col=brewer.pal(9, "Set1")[1],
	# lwd=1.5,
	bty="n")
dev.off()


# Plots heatmap of eigengenes for high-variance patients
pdf("co-expression/eigenNetw/v3/plots/CAD_eigengene_heatmap.pdf", width=8)
mat = mat[apply(mat, 1, sd) > sd_cutoff, ]

heatmap.2(mat,
	trace="none",
	col=colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
	breaks=seq(-0.2, 0.2, length.out=101),
	xlab="CAD eigengenes",
	ylab="High-variance patients (75th percentile)",
	key.title="",
	key.xlab="Eigengene",
	tracecol="black",
	margins=c(4, 16),
	cexRow=0.3
)
dev.off()


#
bn = fges(variables, maxDegree=100)  # Fast-greedy equivalence search


bootstrap_m = 1000
bn_bootstrap = lapply(1:bootstrap_m, function(k) {
	if (k %% 10 == 0) {
		message("iter ", k, " out of ", bootstrap_m)
	}
	# Sample with replacement from patient IDs
	idx = sample(nrow(variables), replace=TRUE)

	out = fges(variables[idx, ], maxDegree=100)

	return(out)
})


# Calculate average bootstrapped adjacency matrix
bn_bootstrap_adjmat = lapply(bn_bootstrap, function(bn) {
	adj_mat = as(bn$graphNEL, "matrix")
	return(adj_mat)
})

boot_adjmat = Reduce("+", bn_bootstrap_adjmat) / length(bn_bootstrap_adjmat)
boot_p = 1 - boot_adjmat

library(magicaxis)

# plot(sort(boot_p)[1:50])

pdf("co-expression/eigenNetw/v3/plots/supernetwork_bootstrap.pdf", width=5)
n_edges = 50
sig_colors = c("white", brewer.pal(4, "YlGn"))
sig_levels = c(0.05, 0.1, 0.2, 0.5)



# To edge list
boot_p_tab = melt(boot_p)
colnames(boot_p_tab) = c("from", "to", "p")

# Order by significance
boot_p_tab = boot_p_tab[order(boot_p_tab$p), ]

p = boot_p_tab$p

pts_col = rep(sig_colors[1], n_edges)
for (i in 1:length(sig_levels)) {
	sig = rev(sig_levels)[i]
	pts_col[p < sig] = sig_colors[i + 1]
}

plot(-log10(p)[1:n_edges],
	type="n",
	bty="l",
	ylab=expression("-log"[10] * " p (bootstrap)"),
	xlab="Supernetwork edge rank"
)

sapply(sig_levels, function(sig) {
	abline(h=-log10(sig), col="grey", lty=2)
	text(45, -log10(sig),
		pos=3,
		labels=paste("p < ", sig),
		col="grey"
	)
})

points(-log10(p)[1:n_edges],
	bg=pts_col,
	pch=21
)

text(1:k, -log10(p[1:k]),
	pos=4,
	labels=paste(boot_p_tab$from, boot_p_tab$to, sep=" -> ")[1:k])

dev.off()



# Parse GWAS gene enrichment results
fdr = 0.1
GWAS_enrich = list()
GWAS_enrich$CAD = p.adjust(module_tab$CAD_pval, method="BH") < fdr
GWAS_enrich$T2D = p.adjust(module_tab$Type.2.diabetes_pval, method="BH") < fdr
GWAS_enrich$LDL = p.adjust(module_tab$LDL.cholesterol_pval, method="BH") < fdr
GWAS_enrich$HDL = p.adjust(module_tab$HDL.cholesterol_pval, method="BH") < fdr
GWAS_enrich$Total.Cholesterol = p.adjust(module_tab$Cholesterol..total_pval, method="BH") < fdr
GWAS_enrich$BMI = p.adjust(module_tab$Body.mass.index_pval , method="BH") < fdr
GWAS_enrich$WHRadjBMI = p.adjust(module_tab$Waist.to.hip.ratio.adjusted.for.body.mass.index_pval , method="BH") < fdr
GWAS_enrich$Fasting.Glucose = p.adjust(module_tab$Fasting.plasma.glucose_pval, method="BH") < fdr
GWAS_enrich$Blood.pressure = p.adjust(module_tab$Blood.pressure_pval, method="BH") < fdr
GWAS_enrich$Triglycerides = p.adjust(module_tab$Triglycerides_pval, method="BH") < fdr

gwas_bool = Reduce(rbind, GWAS_enrich)
gwas_bool = rbind(apply(gwas_bool, 2, function(x) !any(x)), gwas_bool)  # add dummy variable for first row if empty


# Pooled GWAS results
GWAS_enrich_comb = list()
GWAS_enrich_comb$CAD = GWAS_enrich$CAD
GWAS_enrich_comb$Blood.lipids =
	GWAS_enrich$LDL |
	GWAS_enrich$HDL |
	GWAS_enrich$Total.Cholesterol |
	GWAS_enrich$Triglycerides

GWAS_enrich_comb$Glucose.Metabolism =
	GWAS_enrich$T2D |
	GWAS_enrich$Fasting.Glucose

GWAS_enrich_comb$Obesity = GWAS_enrich$BMI |
	GWAS_enrich$WHRadjBMI

GWAS_enrich_comb$Blood.pressure = GWAS_enrich$Blood.pressure


gwas_comb_bool = Reduce(rbind, GWAS_enrich_comb)
gwas_comb_bool = rbind(apply(gwas_comb_bool, 2, function(x) !any(x)), gwas_comb_bool)  # add dummy variable for first row if empty



# Supernetwork visualization
# ----------------------------------------


# boot_cuttoff = 0.05
# boot_cuttoff = 0.1
boot_cuttoff = 0.2
# boot_cuttoff = 0.5
# boot_cuttoff = 0.7
# boot_cuttoff = 0.8

boot_p_tab_sel = boot_p_tab[boot_p_tab$p < boot_cuttoff, ]

# String encoding of node IDs
boot_p_tab_sel$from = as.character(boot_p_tab_sel$from)
boot_p_tab_sel$to = as.character(boot_p_tab_sel$to)

g = graph_from_edgelist(
	as.matrix(boot_p_tab_sel)[, 1:2]
)

# Hierarcical layout
lay = layout_with_sugiyama(g,
	maxiter=1000
)

# Scale layout
xscale = 1.2
yscale = 2.0
lay$layout[, 1] =  xscale * lay$layout[, 1]
lay$layout[, 2] =  yscale * lay$layout[, 2]


# Wrapper function ofr igraph network plot
plotNetw = function(g, layout, ...) {
	plot(g,
		# main="Tissue",
		layout=lay$layout,
		edge.width=1.2,
		edge.arrow.size=0.5,
		# edge.color=rgb(50, 50, 50, maxColorValue=255),
		edge.color="black",
		vertex.label.family="Helvetica",
		# vertex.label.color="black",
		# vertex.label.color="grey",
		vertex.label.color=rgb(150, 150, 150, maxColorValue=255),
		vertex.label.degree=pi/2,  # labels below vertex
		vertex.label.dist=10.0,
		vertex.size=100,
		rescale=FALSE,
		xlim=range(lay$layout[, 1]),
		ylim=range(lay$layout[, 2]),
		...
	)
}


# pdf("co-expression/eigenNetw/v3/plots/netw_tissue.pdf", width=12, height=12)
pdf("co-expression/eigenNetw/v3/plots/netw_tissue.pdf", width=6, height=18)
par(mfrow=c(3, 1))

# Tissue composition pie charts
# ---------------------
tissues = colnames(module_tab)[4:10]
pie_col = brewer.pal(9, "Set1")[c(1, 2, 7, 5, 3, 8, 4)]

V(g)$pie.color = rep(list(pie_col),
	length(V(g))
)
V(g)$shape = rep("pie", length(V(g)))


# module_tab
mod_idx = match(V(g)$name, module_tab[, 1])

V(g)$pie = as.list(
	as.data.frame(
		t(module_tab[mod_idx, 4:10])
	)
)


plotNetw(g, layout=lay, main="Tissue")

legend("bottomleft",
	legend=tissues,
	# pt.bg=pie_col,
	col=pie_col,
	bty="n",
	pch=15
)


# GWAS gene enrichment
# -------------------------------------------------
# pie_col = c("white", brewer.pal(8, "Set1"))
# pie_col[c(3, 5, 6)] =  pie_col[c(7, 9, 3)]

pie_col = c("white", brewer.pal(8, "Dark2")[c(2, 6, 1, 4, 3)])

V(g)$shape = rep("pie", length(V(g)))
V(g)$pie = as.list(data.frame(gwas_comb_bool[, mod_idx] + 0))  # pie fractions

V(g)$pie.color = rep(
	list(pie_col),
	length(V(g)))

# Plots
first_pie = sapply(V(g)$pie, function(x) {
	return(match(1, x))
})
V(g)$color = pie_col[first_pie]

plotNetw(g, layout=lay, main="GWAS")

legend("bottomleft",
	legend=c("None", names(GWAS_enrich_comb)),
	pch=15,
	col=pie_col,
	bty="n"
)


# Secreted protein enrichment
# -------------------------------------------------

vert_col = list()

# Enrichment for endocrine factors
fdr = 0.1
vert_col$secreted = rep("white", length(V(g)))
vert_col$secreted[
	p.adjust(module_tab$secreted_protein_pval, method="BH") < fdr
] = brewer.pal(9, "Set1")[2]

V(g)$color = vert_col$secreted[mod_idx]

V(g)$shape = rep("circle", length(V(g)))
plotNetw(g, layout=lay, main="Secreted proteins")

dev.off()


# Store workspace for recovery of exact results
save.image(file="co-expression/eigenNetw/Rworkspace/bayesNet3.RData")


# Plot of previous eigengene supernetwork
# ----------------------------------------------

# Load previous network into environment variable
# This contains the global bayesian network and layout previously computed
# see 'bayesNet2.R'
netw2 = environment()
load("co-expression/eigenNetw/Rworkspace/bayesNet2.RData", netw2,
	verbose=TRUE)


hideVertices = function(g, vertices) {
	# Vertices
	# V(g)$shape = rep("circle", length(V(g)))
	V(g)$shape[vertices] = "none"

	# V(g)$label = names(V(g))
	V(g)$label[vertices] = NA

	return(g)
}

showAllVertices = function(g, shape="circle") {
	V(g)$shape = rep(shape, length(V(g)))
	V(g)$label = names(V(g))

	E(g)$lty = rep(1, length(E(g)))
	return(g)
}


plotNetw = function(g, name, mod_sel, lay, shape="circle", ...) {
	# if (mod_sel == "CT") {
	# 	g = showAllVertices(g, shape)
	# 	g = hideVertices(g, module_tab$purity > 0.95)
	# } else if (mod_sel == "TS"){
	# 	g = showAllVertices(g, shape)
	# 	g = hideVertices(g, module_tab$purity <= 0.95)
	# } else if (mod_sel == "all") {
	# 	g = showAllVertices(g, shape=shape)
	# # } else if (mod_sel == "allpie") {
	# # 	g = showAllVertices(g, shape="pie")
	# }

	# Hide vertices without any edges
	g = hideVertices(g, degree(g) == 0)

	pdf(paste0("co-expression/eigenNetw/v3/igraph/by_", name, "_", mod_sel, ".pdf"),
		width=20, height=20)
	plot(g,
		# edge.arrow.size=0.5,
		edge.width=3,
		edge.arrow.size=0.5,
		# edge.color=rgb(0.9, 0.9, 0.9),
		vertex.label.cex=0.7,
		vertex.label.family="Helvetica",
		vertex.label.color="black",
		layout=lay,
		...
	)
	dev.off()
}

# Color outgoing edges based on vertex color, unless white
colorOutEdge = function(g) {
	E(g)$color = rep(rgb(0.9, 0.9, 0.9), length(E(g)))  # reset edge colors

	from = get.edgelist(g)[, 1]
	cols = V(g)$color[as.integer(from)]

	E(g)$color[cols != "white"] = cols[cols != "white"]
	return(g)
}

resetEdgeColor = function(g) {
	E(g)$color = rep(rgb(0.9, 0.9, 0.9), length(E(g)))
	return(g)
}


V(netw2$g)$shape = rep("circle", length(V(netw2$g)))
V(netw2$g)$shape[module_tab$purity < 0.95] = "square"

 # = rep("circle", length(V(netw2$g)))

# Reset node color
V(netw2$g)$color = rep("white", length(V(netw2$g)))


fdr = 0.01
V(netw2$g)$color[data.frame(pheno_padj)$SYNTAX < fdr] = brewer.pal(9, "Pastel1")[1]
V(netw2$g)$color[data.frame(pheno_padj)$DUKE < fdr] = brewer.pal(9, "Pastel1")[2]
V(netw2$g)$color[data.frame(pheno_padj)$Case.Ctrl.DEG < fdr] = brewer.pal(9, "Pastel1")[3]

# All modules
V(netw2$g)$color[cad_modules] = rgb(50, 50, 50, maxColorValue=255)
# V(netw2$g)$color[cad_modules] = rgb(100, 100, 100, maxColorValue=255)
# V(netw2$g)$color[cad_modules] = rgb(0, 0, 0, maxColorValue=255)

netw2$g = colorOutEdge(netw2$g)

plotNetw(netw2$g, name="tissue", mod_sel="all", lay=netw2$lay, vertex.label=NA)
# plotNetw(netw2$g, name="tissue", mod_sel="CT", lay=netw2$lay)
# plotNetw(netw2$g, name="tissue", mod_sel="TS", lay=netw2$lay)

