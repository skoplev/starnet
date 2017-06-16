options(java.parameters = "-Xmx8g")  # Max Java memory heap size, for rcausal

rm(list=ls())

library(rcausal)
library(data.table)
library(igraph)
library(RColorBrewer)
library(squash)

data_dir = "/Users/sk/DataProjects/cross-tissue"  # root of data directory
setwd("/Users/sk/Google Drive/projects/cross-tissue")

source("src/parse.R")


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

variables = cbind(
	between$bwnet$eigengenes
)

bn = fges(variables, maxDegree=100)  # Fast-greedy equivalence search

g = graph_from_graphnel(bn$graphNEL)

epsilon = 10^-30  # small number

# Add featueres
# V(g)$log_P_cad = -log10(module_tab$cad_pval)
V(g)$mod_size_square = sqrt(module_tab$mod_size) + epsilon

V(g)$mod_size_square = sqrt(module_tab$mod_size) + epsilon

V(g)$purity = module_tab$purity
V(g)$cross_tissue = as.character(module_tab$purity < 0.95)

V(g)$frac_AOR = module_tab$AOR / module_tab$mod_size
V(g)$frac_BLOOD = module_tab$BLOOD / module_tab$mod_size
V(g)$frac_LIV = module_tab$LIV / module_tab$mod_size
V(g)$frac_MAM = module_tab$MAM / module_tab$mod_size
V(g)$frac_SKLM = module_tab$SKLM / module_tab$mod_size
V(g)$frac_SF = module_tab$SF / module_tab$mod_size
V(g)$frac_VAF = module_tab$VAF / module_tab$mod_size

# Primary and secondary tissues
tissues = c("AOR", "BLOOD", "LIV", "MAM", "SKLM", "SF", "VAF")
V(g)$primary_tissue = tissues[
	apply(module_tab[, tissues], 1, which.max)
]

V(g)$secondary_tissue = tissues[
	apply(module_tab[, tissues], 1, function(vals) {
		i = order(vals, decreasing=TRUE)[2]
		if (vals[i] == 0) {
			return(NA)
		} else {
			return(i)
		}
	})
]

# Secondary tissue for cross-tissue modules only
V(g)$secondary_tissue_cross_tissue = V(g)$secondary_tissue
V(g)$secondary_tissue_cross_tissue[module_tab$purity > 0.95] = NA

V(g)$log_pval_CAD2 = -log10(module_tab$CAD_pval)  #  CAD enrichment
V(g)$log_pval_DUKE = -log10(module_tab$pval_DUKE)  #  DUKE correlation
V(g)$log_pval_BMI = -log10(module_tab$pval_BMI)
V(g)$log_pval_LDL = -log10(module_tab$pval_LDL)
V(g)$log_pval_chol = -log10(module_tab$pval_p_chol)


V(g)$log_n_within_mod_endocrines = log10(module_tab$n_within_mod_endocrines)
V(g)$log_n_within_mod_endocrines[is.infinite(V(g)$log_n_within_mod_endocrines)] = 0

write_graph(g,
	file="co-expression/eigenNetw/v2/bayes_net.gml",
	format="gml")


# Bayesian network of cross-tissue only
# ---------------------------------------------

# idx = module_tab$purity < 0.95

idx = module_tab$CAD_qval < 0.1 |
	module_tab$pval_syntax_score < 0.05 |
	module_tab$pval_DUKE < 0.05 |
	module_tab$pval_ndv < 0.05 |
	module_tab$pval_lesions < 0.05
	# module_tab$pval_LDL < 0.05 |
	# module_tab$pval_p_chol < 0.05


# idx = module_tab$purity < 0.95
variables = cbind(
	between$bwnet$eigengenes[idx]
)

bn2 = fges(variables, maxDegree=100)  # Fast-greedy equivalence search
module_tab_sub = module_tab[idx, ]

g = graph_from_graphnel(bn2$graphNEL)

epsilon = 10^-30  # small number
# epsilon = 10^-6  # small number

# Add featueres
V(g)$mod_size_square = sqrt(module_tab_sub$mod_size) + epsilon
V(g)$mod_size_square = sqrt(module_tab_sub$mod_size) + epsilon

V(g)$purity = module_tab_sub$purity
V(g)$cross_tissue = as.character(module_tab_sub$purity < 0.95)

V(g)$frac_AOR = module_tab_sub$AOR / module_tab_sub$mod_size
V(g)$frac_BLOOD = module_tab_sub$BLOOD / module_tab_sub$mod_size
V(g)$frac_LIV = module_tab_sub$LIV / module_tab_sub$mod_size
V(g)$frac_MAM = module_tab_sub$MAM / module_tab_sub$mod_size
V(g)$frac_SKLM = module_tab_sub$SKLM / module_tab_sub$mod_size
V(g)$frac_SF = module_tab_sub$SF / module_tab_sub$mod_size
V(g)$frac_VAF = module_tab_sub$VAF / module_tab_sub$mod_size

# Primary and secondary tissues
tissues = c("AOR", "BLOOD", "LIV", "MAM", "SKLM", "SF", "VAF")

V(g)$primary_tissue = tissues[
	apply(module_tab_sub[, tissues], 1, which.max)
]

V(g)$secondary_tissue = tissues[
	apply(module_tab_sub[, tissues], 1, function(vals) {
		i = order(vals, decreasing=TRUE)[2]
		if (vals[i] == 0) {
			return(NA)
		} else {
			return(i)
		}
	})
]

# Secondary tissue for cross-tissue modules only
V(g)$secondary_tissue_cross_tissue = V(g)$secondary_tissue
V(g)$secondary_tissue_cross_tissue[module_tab_sub$purity > 0.95] = NA

V(g)$log_pval_CAD2 = -log10(module_tab_sub$CAD_pval)  #  CAD enrichment
V(g)$log_pval_DUKE = -log10(module_tab_sub$pval_DUKE)  #  DUKE correlation

V(g)$log_n_within_mod_endocrines = log10(module_tab_sub$n_within_mod_endocrines)
V(g)$log_n_within_mod_endocrines[is.infinite(V(g)$log_n_within_mod_endocrines)] = 0

write_graph(g,
	file="co-expression/eigenNetw/v2/bayes_net2.gml",
	format="gml")



# Plot graphs
# ------------------------

# Blends two vectors of colors
blendColors = function(c1, c2) {
	if (length(c1) != length(c2)) {
		stop("Color vector mismatch")
	}

	sapply(1:length(c1), function(k) {
		rgb_col_blend = (col2rgb(c1[k]) + col2rgb(c2[k])) / 2

		# rgb_col_blend = (col2rgb(c1[k]) + col2rgb(c2[k]))
		# rgb_col_blend[rgb_col_blend > 255] = 255

		# Convert to R color
		col = rgb(rgb_col_blend[1], rgb_col_blend[2], rgb_col_blend[3],
			maxColorValue=255)

		return(col)
	})
}


tissues = c("AOR", "BLOOD", "SKLM", "VAF", "MAM", "LIV",  "SF")
tissue_cols = brewer.pal(9, "Set1")[-6]


# lay = layout_with_kk(g)
# lay = layout_with_dh(g)  # David-Harel layout
# lay = layout_with_gem(g)
lay = layout_with_fr(g, niter=20000)  # Furctherman-Reingold


# Size of nodes
V(g)$size = V(g)$mod_size_square / 15 + 1.5



# V(g)$shape[module_tab$purity > 0.95] = "none"  # hide tissue-specific

hideVertices = function(g, vertices) {
	# Vertices
	V(g)$shape = rep("circle", length(V(g)))
	V(g)$shape[vertices] = "none"

	V(g)$label = names(V(g))
	V(g)$label[vertices] = NA

	# # Edges
	# E(g)$lty = rep(1, length(E(g)))
	# E(g)$lty[
	# 	vertices[as.integer(get.edgelist(g)[, 1])]
	# ] = 0

	return(g)
}

showAllVertices = function(g) {
	V(g)$shape = rep("circle", length(V(g)))
	V(g)$label = 1:length(V(g))
	V(g)$label = names(V(g))

	E(g)$lty = rep(1, length(E(g)))
	return(g)
}

plotNetw = function(g, col_feature, mod_sel) {
	if (mod_sel == "CT") {
		g = showAllVertices(g)
		g = hideVertices(g, module_tab$purity > 0.95)
	} else if (mod_sel == "TS"){
		g = showAllVertices(g)
		g = hideVertices(g, module_tab$purity <= 0.95)
	} else if (mod_sel == "all") {
		g = showAllVertices(g)
	}

	# Hide vertices without any edges
	g = hideVertices(g, degree(g) == 0)

	pdf(paste0("co-expression/eigenNetw/v2/igraph/by_", col_feature, "_", mod_sel, ".pdf"),
		width=20, height=20)
	plot(g,
		# edge.arrow.size=0.5,
		edge.width=3,
		edge.arrow.size=0.5,
		# edge.color=rgb(0.9, 0.9, 0.9),
		vertex.label.cex=0.7,
		vertex.label.family="Helvetica",
		vertex.label.color="black",
		layout=lay
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




# Vertex colors
vert_col = list()
vert_col$tissue = tissue_cols[
	as.integer(factor(V(g)$primary_tissue, levels=tissues))
]

# Color based on cross-tissue
vert_col$cross_tissue = rep("white", length(V(g)))
# vert_col$cross_tissue[module_tab$purity < 0.95] = rgb(153, 153, 153, maxColorValue=255)
vert_col$cross_tissue[module_tab$purity < 0.95] = brewer.pal(8, "Dark2")[1]

vert_col$tissue_specific = rep("white", length(V(g)))
# vert_col$tissue_specific[module_tab$purity >= 0.95] = rgb(153, 153, 153, maxColorValue=255)
vert_col$tissue_specific[module_tab$purity >= 0.95] = brewer.pal(8, "Dark2")[6]

vert_col$cad = rep("white", length(V(g)))
vert_col$cad[module_tab$CAD_qval < 0.1] = brewer.pal(9, "Set1")[1]  # Red

# x = -log10(module_tab$CAD_qval)
# x[x>5] = 5
# vert_col$cad = cmap(x,
# 	# Color map
# 	makecmap(x,
# 		n=100,
# 		colFn=colorRampPalette(c("white", brewer.pal(9, "Reds")))
# 	)
# )

vert_col$DUKE = rep("white", length(V(g)))
vert_col$DUKE[module_tab$pval_DUKE < 0.05] = brewer.pal(9, "Set1")[3]  # Green


# BLend colors
# vert_col$cad_cross_tissue = blendColors(vert_col$cad, vert_col$cross_tissue)
# vert_col$DUKE_cross_tissue = blendColors(vert_col$DUKE, vert_col$cross_tissue)

# vert_col$cad_cross_tissue = vert_col$cad
# vert_col$DUKE_cross_tissue = vert_col$DUKE



# # Edge colors based on source tissue
# E(g)$color = tissue_cols[
# 	as.integer(
# 		factor(
# 			V(g)$primary_tissue[as.integer(get.edgelist(g)[, 1])],
# 			levels=tissues
# 		)
# 	)
# ]

g = colorOutEdge(g)

# Some counts
sum(degree(g) > 0 & module_tab$purity < 0.95)  # CT shown
sum(degree(g) > 0 & module_tab$purity >= 0.95)  # TS shown


plotNetw(g, "tissue", "all")

# Reset edge colors
E(g)$color = rep(rgb(0.9, 0.9, 0.9), length(E(g)))

# Edge colors based on source tissue
E(g)$color = tissue_cols[
	as.integer(
		factor(
			V(g)$primary_tissue[as.integer(get.edgelist(g)[, 1])],
			levels=tissues
		)
	)
]



# Color by CAD enrichment
V(g)$color = vert_col[["cross_tissue"]]
g = colorOutEdge(g)
plotNetw(g, "cross_tissue", "all")

V(g)$color = vert_col[["tissue_specific"]]
g = colorOutEdge(g)
plotNetw(g, "tissue_specific", "all")

V(g)$color = vert_col[["DUKE"]]
plotNetw(g, "DUKE", "CT")
plotNetw(g, "DUKE", "TS")
plotNetw(g, "DUKE", "all")

V(g)$color = vert_col[["cad"]]
plotNetw(g, "cad", "CT")
plotNetw(g, "cad", "TS")
plotNetw(g, "cad", "all")




# Hierarcical graph
# ------------------------
# tissues
tissue_cols = brewer.pal(9, "Set1")[-6]



# Tissue node colors
V(g)$color = tissue_cols[
	as.integer(factor(V(g)$primary_tissue, levels=tissues))
]

# Color based on cross-tissue
V(g)$color = rep("white", length(V(g)))
V(g)$color[module_tab$purity < 0.95] = brewer.pal(9, "Set1")[2]  # Blue


# Node colors based on additional features
# x = V(g)$log_pval_CAD2
# x = V(g)$log_pval_DUKE
# x = V(g)$log_pval_BMI
x = V(g)$log_pval_LDL
# x = V(g)$log_pval_LDL
x[x>5] = 5

# colormap
map = makecmap(x, n=100,
	colFn=colorRampPalette(brewer.pal(9, "Blues"))
) 

V(g)$color = cmap(x, map)


# Size of nodes
V(g)$size = V(g)$mod_size_square / 30 + 1.0

# Edge colors based on source tissue
E(g)$color = tissue_cols[
	as.integer(
		factor(
			V(g)$primary_tissue[as.integer(get.edgelist(g)[, 1])],
			levels=tissues
		)
	)
]









lay = layout_with_sugiyama(g,
	# hgap=30,
	# vgap=30,
	# maxiter=10000
	# maxiter=1,
	attributes="all"
)




# lay$layout

# layout = lay$layout

# moveOverlapHorizontal = function(layout, sizes) {
# 	# which modules overlap?

# 	# Calculate distance matrix
# 	dmat = dist(layout)
# 	# dmat = as.matrix(dmat)

# 	which(dmat == min(dmat))
# 	which.min(dmat)
# }

# sort(hub_score(g)$vector)

# # betweenness(g)
# # degree(g)
# # hierarchy(g)


plot(lay$layout[, 2], -log10(module_tab$CAD_pval))
cor.test(lay$layout[, 2], -log10(module_tab$CAD_pval))


plot(lay$layout[, 2], -log10(module_tab$secreted_protein_pval))
cor.test(lay$layout[, 2], -log10(module_tab$secreted_protein_pval))


layer = max(lay$layout[, 2]) - lay$layout[, 2] + 1


par(mfrow=c(4, 2))
plot(layer, -log10(module_tab$pval_DUKE),
	col=tissue_cols[
		as.integer(factor(V(g)$primary_tissue, levels=tissues))
	],
	pch=16
)
abline(h=-log10(0.1), col="grey", lty=2)

# density()


plot(layer, -log10(module_tab$CAD_qval),
	col=tissue_cols[
		as.integer(factor(V(g)$primary_tissue, levels=tissues))
	],
	pch=16
)
abline(h=-log10(0.1), col="grey", lty=2)

plot(layer, -log10(module_tab$secreted_protein_qval),
	col=tissue_cols[
		as.integer(factor(V(g)$primary_tissue, levels=tissues))
	],
	pch=16
)
abline(h=-log10(0.1), col="grey", lty=2)


plot(layer, 1 - module_tab$purity,
	col=tissue_cols[
		as.integer(factor(V(g)$primary_tissue, levels=tissues))
	],
	pch=16
)
abline(h=0.05, col="grey", lty=2)


x = -log10(module_tab$pval_BMI)
x[x>10] = 10
plot(layer, x,
	col=tissue_cols[
		as.integer(factor(V(g)$primary_tissue, levels=tissues))
	],
	ylab="BMI",
	pch=16
)
abline(h=-log10(0.1), col="grey", lty=2)


x = -log10(module_tab$pval_LDL)
x[x>10] = 10
plot(layer, x,
	col=tissue_cols[
		as.integer(factor(V(g)$primary_tissue, levels=tissues))
	],
	ylab="LDL",
	pch=16
)
abline(h=-log10(0.1), col="grey", lty=2)






cor.test(lay$layout[, 2], -log10(module_tab$pval_DUKE))

# plot(lay$layout[, 2], -log10(module_tab$pval_syntax_score))
# cor.test(lay$layout[, 2], -log10(module_tab$pval_syntax_score))


# plot(degree(g, mode="out"), -log10(module_tab$CAD_pval))
# cor.test(degree(g, mode="out"), -log10(module_tab$CAD_pval))


# cor.test(betweenness(g), -log10(module_tab$CAD_pval))

# plot(betweenness(g), module_tab$purity)

# cor.test(betweenness(g), module_tab$purity)


pdf("co-expression/eigenNetw/v2/igraph/h1.pdf", width=20, height=20)
plot.igraph(
	g,
	# lay$extd_graph
	# vertex.size=1.0,
	edge.arrow.size=0.3,
	edge.color=rgb(0.9, 0.9, 0.9),
	vertex.label.cex=0.4,
	vertex.label.color="black",
	asp=0.5,  # aspect ratio
	layout=lay$layout
)
dev.off()


# pdf("co-expression/eigenNetw/v2/igraph/h1.pdf", width=20, height=20)
# plot.igraph(
# 	# g,
# 	lay$extd_graph,
# 	# vertex.size=1.0,
# 	edge.arrow.size=0.3,
# 	edge.color=rgb(0.9, 0.9, 0.9),
# 	vertex.label.cex=0.4,
# 	vertex.label.color="black",
# 	asp=0.5  # aspect ratio
# 	# layout=lay$layout
# )
# dev.off()


library(gplots)

heatmap.2(as.matrix(as_adj(g)), col=c("white", "black"),
	trace="none")
