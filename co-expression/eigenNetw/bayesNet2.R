options(java.parameters = "-Xmx8g")  # Max Java memory heap size, for rcausal

rm(list=ls())

library(rcausal)
library(data.table)
library(igraph)
library(RColorBrewer)
library(squash)
library(qvalue)

data_dir = "/Users/sk/DataProjects/cross-tissue"  # root of data directory
setwd("~/Google Drive/projects/STARNET/cross-tissue")

source("src/parse.R")
source("src/base.R")


# Order of tissues
tissues = c("AOR", "BLOOD", "SKLM", "VAF", "MAM", "LIV",  "SF")
tissue_cols = brewer.pal(9, "Set1")[-6]


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

# Load significant endocrine-eigengene correlations
endocrines = read.csv("co-expression/tables/CT_endocrines_TS_interactions.csv")


variables = cbind(
	between$bwnet$eigengenes
)

bn = fges(variables, maxDegree=100)  # Fast-greedy equivalence search

g = graph.from.graphNEL(bn$graphNEL)

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
# tissues = c("AOR", "BLOOD", "LIV", "MAM", "SKLM", "SF", "VAF")
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

g = igraph.from.graphNEL(bn2$graphNEL)

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
# tissues = c("AOR", "BLOOD", "LIV", "MAM", "SKLM", "SF", "VAF")

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

# V(g)$shape[module_tab$purity > 0.95] = "none"  # hide tissue-specific

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

plotNetw = function(g, name, mod_sel, shape="circle", ...) {
	if (mod_sel == "CT") {
		g = showAllVertices(g, shape)
		g = hideVertices(g, module_tab$purity > 0.95)
	} else if (mod_sel == "TS"){
		g = showAllVertices(g, shape)
		g = hideVertices(g, module_tab$purity <= 0.95)
	} else if (mod_sel == "all") {
		g = showAllVertices(g, shape=shape)
	# } else if (mod_sel == "allpie") {
	# 	g = showAllVertices(g, shape="pie")
	}

	# Hide vertices without any edges
	g = hideVertices(g, degree(g) == 0)

	pdf(paste0("co-expression/eigenNetw/v2/igraph/by_", name, "_", mod_sel, ".pdf"),
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



# lay = layout_with_kk(g)
# lay = layout_with_dh(g)  # David-Harel layout
# lay = layout_with_gem(g)
lay = layout_with_fr(g, niter=20000)  # Furctherman-Reingold


# Size of nodes
V(g)$size = V(g)$mod_size_square / 15 + 1.5



# Vertex colors
vert_col = list()
vert_col$tissue = tissue_cols[
	as.integer(factor(V(g)$primary_tissue, levels=tissues))
]

# Color based on cross-tissue
vert_col$cross_tissue = rep("white", length(V(g)))
# vert_col$cross_tissue[module_tab$purity < 0.95] = rgb(153, 153, 153, maxColorValue=255)
vert_col$cross_tissue[module_tab$purity < 0.95] = brewer.pal(8, "Dark2")[5]

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

angio_pval =data.frame(
	DUKE=module_tab$pval_DUKE,
	SYNTAX=module_tab$pval_syntax_score,
	ndv=module_tab$pval_ndv,
	ndv=module_tab$pval_lesions
)

# apply(angio_pval, 1, min) < 0.95
vert_col$angio_scores = rep("white", length(V(g)))
vert_col$angio_scores[
	apply(angio_pval, 1, min) < 0.05
] = brewer.pal(9, "Set1")[1]


# Enrichment for endocrine factors
vert_col$secreted = rep("white", length(V(g)))
vert_col$secreted[
	p.adjust(module_tab$secreted_protein_pval, method="BH") < 0.1
] = brewer.pal(9, "Set1")[2]


# p.adjust(module_tab$pval_DUKE, method="BH")
# sum(p.adjust(module_tab$pval_DUKE, method="BH")< 0.2)
# sum(p.adjust(module_tab$pval_syntax_score, method="BH")< 0.2)
# sum(p.adjust(module_tab$pval_ndv, method="BH")< 0.2)

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

tissue_specific = module_tab$purity >= 0.95
cross_tissue = module_tab$purity < 0.95


# Tissues, cross-tissue, tissue-specific
g = colorOutEdge(g)
plotNetw(g, "tissue", "all")

# plotNetw(g, "test", "all")


# Reset edge colors
E(g)$color = rep(rgb(0.9, 0.9, 0.9), length(E(g)))

V(g)$color = vert_col[["cross_tissue"]]
g = colorOutEdge(g)
plotNetw(g, "cross_tissue", "all")

V(g)$color = vert_col[["tissue_specific"]]
g = colorOutEdge(g)
plotNetw(g, "tissue_specific", "all")


# GWAS enrichment

# GWAS enrichment collection
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


# GWAS coverage for each module 
gwas = fread("~/DataBases/GWAS/gwas_catalog_v1.0-associations_e88_r2017-05-29.tsv")

traits = c(
	"Type 2 diabetes",
	# "Lipid metabolism phenotypes",
	"LDL cholesterol",
	"HDL cholesterol",
	"Cholesterol, total",
	"Body mass index",
	"Waist-to-hip ratio adjusted for body mass index",
	"Fasting plasma glucose",
	"Blood pressure",
	"Triglycerides"
)

gwas_genes = list()
gwas_genes = lapply(traits, function(trait) getGWAS(gwas, trait))
names(gwas_genes) = traits

gwas_genes = c(list(CAD=getCADGenes(data_dir)), gwas_genes)

# check trait names manually
data.frame(n1=names(GWAS_enrich), n2=names(gwas_genes))

# Calculate overlap in enriched modules
frac_found = sapply(1:length(GWAS_enrich), function(i) {
	# Get gene symbols of all enriched modules
	sig_gene_symbols_CT = between$meta_genes$gene_symbol[
		between$clust %in% which(GWAS_enrich[[i]] & cross_tissue)
	]

	sig_gene_symbols_TS = between$meta_genes$gene_symbol[
		between$clust %in% which(GWAS_enrich[[i]] & tissue_specific)
	]

	found_gwas_genes_CT = gwas_genes[[i]] %in% sig_gene_symbols_CT
	found_gwas_genes_TS = gwas_genes[[i]] %in% sig_gene_symbols_TS

	return(c(
		mean(found_gwas_genes_CT & !found_gwas_genes_TS),
		mean(found_gwas_genes_TS & found_gwas_genes_CT),
		mean(found_gwas_genes_TS & ! found_gwas_genes_CT)
	))
})

rownames(frac_found) = c("CT", "CT+TS", "TS")
colnames(frac_found) = paste0(names(GWAS_enrich), " (n=", sapply(gwas_genes, length), ")")


# Order by total fraction found
frac_found = frac_found[, order(apply(frac_found, 2, sum), decreasing=TRUE)]

# bar_col = brewer.pal(8, "Set2")
bar_col = gray.colors(3)

pdf("co-expression/plots/CAD_in_enriched_modules.pdf", width=5, height=6.5)
par(mar=c(14, 4.1, 4.1, 2.1))
barplot(frac_found * 100, las=2,
	ylab="Loci in enriched RGNs (%)",
	xlab="GWAS trait loci",
	col=bar_col
)

legend("topright", legend=c("CT", "CT+TS", "TS"), pch=22, pt.bg=bar_col)
dev.off()

# Table of CAD genes in enriched modules
i = 1  # CAD
sig_gene_symbols_CT = between$meta_genes$gene_symbol[
	between$clust %in% which(GWAS_enrich[[i]] & cross_tissue)
]
gwas_genes[[i]][gwas_genes[[i]] %in% sig_gene_symbols_CT]

idx = 
	(between$clust %in% which(GWAS_enrich[[i]] & cross_tissue)) &   # transripts in CAD enriched CT modules
	(between$meta_genes$gene_symbol %in% gwas_genes[[i]])  # which are CAD genes

cad_enrich_CT_genes = cbind(
	between$meta_genes[idx, ],
	list(module=between$clust[idx])
)
write.csv(cad_enrich_CT_genes, file="co-expression/tables/CAD_enriched_CT.csv")


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



# pie_col = c("white", brewer.pal(9, "Set1"))
# pie_col = c("white", brewer.pal(8, "Dark2")[-8], brewer.pal)
# pie_col = c("white", brewer.pal(9, "Set1")[c(-6, -9)], brewer.pal(8, "Dark2"))
# pie_col = c("white", brewer.pal(9, "Set1")[1], brewer.pal(8, "Set2")[-8], brewer.pal(9, "Accent"))


# GWAS enrichment network plot (combined categories)
pie_col = c("white", brewer.pal(8, "Set1"))
pie_col[c(3, 5, 6)] =  pie_col[c(7, 9, 3)]
V(g)$shape = rep("pie", length(V(g)))
V(g)$pie = as.list(data.frame(gwas_comb_bool + 0))  # pie fractions

V(g)$pie.color = rep(
	list(pie_col),
	length(V(g)))

# Plots
first_pie = sapply(V(g)$pie, function(x) {
	return(match(1, x))
})
V(g)$color = pie_col[first_pie]
g = colorOutEdge(g)

g = resetEdgeColor(g)
V(g)$color = vert_col[["angio_scores"]]
edge_endocrine_support = apply(get.edgelist(g), 1, paste, collapse="_") %in% paste(endocrines$clust, endocrines$target_clust, sep="_")
E(g)$color[edge_endocrine_support] = "black"

plotNetw(g, "GWAS_comb2", "all", shape="pie")

V(g)$color = rep("white", length(V(g)))
V(g)$color[cross_tissue] = pie_col[first_pie][cross_tissue]
g = colorOutEdge(g)
plotNetw(g, "GWAS_comb", "CT", shape="pie")

V(g)$color = rep("white", length(V(g)))
V(g)$color[tissue_specific] = pie_col[first_pie][tissue_specific]
g = colorOutEdge(g)
plotNetw(g, "GWAS_comb", "TS", shape="pie")

# Plot GWAS legend
pdf("co-expression/eigenNetw/v2/igraph/GWAS_comb_legend.pdf")
plot(0, 0, type="n")
legend("topright", legend=c("None", names(GWAS_enrich_comb)), pch=15, col=pie_col)
dev.off()


# GWAS enrichment network plot (all GWAS separately)
pie_col = c("white", brewer.pal(12, "Set3")[c(-2, -9)])
pie_col[c(2, 4)] = pie_col[c(4, 2)]
V(g)$shape = rep("pie", length(V(g)))
V(g)$pie = as.list(data.frame(gwas_bool + 0))  # pie fractions

V(g)$pie.color = rep(
	list(pie_col),
	length(V(g)))

# Plots
first_pie = sapply(V(g)$pie, function(x) {
	return(match(1, x))
})
V(g)$color = pie_col[first_pie]
g = colorOutEdge(g)
plotNetw(g, "GWAS", "all", shape="pie")

V(g)$color = rep("white", length(V(g)))
V(g)$color[cross_tissue] = pie_col[first_pie][cross_tissue]
g = colorOutEdge(g)
plotNetw(g, "GWAS", "CT", shape="pie")

V(g)$color = rep("white", length(V(g)))
V(g)$color[tissue_specific] = pie_col[first_pie][tissue_specific]
g = colorOutEdge(g)
plotNetw(g, "GWAS", "TS", shape="pie")

# Plot GWAS legend
pdf("co-expression/eigenNetw/v2/igraph/GWAS_legend.pdf")
plot(0, 0, type="n")
legend("topright", legend=c("None", names(GWAS_enrich)), pch=15, col=pie_col)
dev.off()



# Combined (min) DUKE, SYTNAX, ndv, and lesions correlations
V(g)$color = vert_col[["angio_scores"]]
g = colorOutEdge(g)
plotNetw(g, "angio_scores", "all")

V(g)$color = rep("white", length(V(g)))
V(g)$color[cross_tissue] = vert_col[["angio_scores"]][cross_tissue]
g = colorOutEdge(g)
plotNetw(g, "angio_scores", "CT")

V(g)$color = rep("white", length(V(g)))
V(g)$color[tissue_specific] = vert_col[["angio_scores"]][tissue_specific]
g = colorOutEdge(g)
plotNetw(g, "angio_scores", "TS")


V(g)$color = vert_col[["secreted"]]
g = colorOutEdge(g)
plotNetw(g, "secreted", "all")

V(g)$color = rep("white", length(V(g)))
V(g)$color[cross_tissue] = vert_col[["secreted"]][cross_tissue]
g = colorOutEdge(g)
plotNetw(g, "secreted", "CT")

V(g)$color = rep("white", length(V(g)))
V(g)$color[tissue_specific] = vert_col[["secreted"]][tissue_specific]
g = colorOutEdge(g)
plotNetw(g, "secreted", "TS")



V(g)$color = vert_col[["DUKE"]]
plotNetw(g, "DUKE", "CT")
plotNetw(g, "DUKE", "TS")
plotNetw(g, "DUKE", "all")

# Color by CAD enrichment
V(g)$color = vert_col[["cad"]]
plotNetw(g, "cad", "CT")
plotNetw(g, "cad", "TS")
plotNetw(g, "cad", "all")



# Plots with different edge colors


# head(endocrines)
# endocrines[endocrines$target_clust == 98, ]
# endocrines[endocrines$target_tissue_primary == "AOR", ]
# endocrines[endocrines$target_tissue_primary == "MAM", ]

# length(table(endocrines$target_clust))
# unique(endocrines$target_clust)

g = resetEdgeColor(g)
V(g)$color = vert_col[["angio_scores"]]
edge_endocrine_support = apply(get.edgelist(g), 1, paste, collapse="_") %in% paste(endocrines$clust, endocrines$target_clust, sep="_")
E(g)$color[edge_endocrine_support] = "black"

plotNetw(g, "endocrine", "TS")
plotNetw(g, "endocrine", "CT")
plotNetw(g, "endocrine", "all")

# unique(c(endocrines$clust, endocrines$target_clust))


get.edgelist(g)[, 1] 


# Some counts
sum(module_tab$purity < 0.95)  # Cross-tissue
sum(cross_tissue)  # Cross-tissue
mean(cross_tissue)
sum(module_tab$purity >= 0.95)  # Tissue-specific
sum(tissue_specific)  # Tissue-specific
mean(tissue_specific)


gwas_sig = apply(Reduce(rbind, GWAS_enrich), 2, any)
sum(gwas_sig)
mean(gwas_sig)

sum(gwas_sig[cross_tissue])
mean(gwas_sig[cross_tissue])

sum(gwas_sig[tissue_specific])
mean(gwas_sig[tissue_specific])


angio_sig = apply(angio_pval, 1, min) < 0.05
sum(angio_sig)
mean(angio_sig)

sum(angio_sig[cross_tissue])
mean(angio_sig[cross_tissue])

sum(angio_sig[tissue_specific])
mean(angio_sig[tissue_specific])

secreted_sig = p.adjust(module_tab$secreted_protein_pval, method="BH") < 0.1
sum(secreted_sig)
mean(secreted_sig)

sum(secreted_sig[cross_tissue])
mean(secreted_sig[cross_tissue])

sum(secreted_sig[tissue_specific])
mean(secreted_sig[tissue_specific])


# Select modules
comb_sel = (gwas_sig & angio_sig & secreted_sig)
# comb_sel = gwas_sig + angio_sig + secreted_sig >= 2
which(comb_sel)


sum(comb_sel)
mean(comb_sel)

sum(comb_sel[cross_tissue])
mean(comb_sel[cross_tissue])

sum(comb_sel[tissue_specific])
mean(comb_sel[tissue_specific])


vert_col$criteria = rep("white", nrow(module_tab))
# vert_col$criteria[comb_sel] = "black"
vert_col$criteria[comb_sel] = brewer.pal(9, "Set1")[1]
# "black"

V(g)$color = vert_col$criteria
g = colorOutEdge(g)
plotNetw(g, "criteria", "all")

V(g)$color = rep("white", length(V(g)))
V(g)$color[cross_tissue] = vert_col$criteria[cross_tissue]
g = colorOutEdge(g)
plotNetw(g, "criteria", "CT")

V(g)$color = rep("white", length(V(g)))
V(g)$color[tissue_specific] = vert_col$criteria[tissue_specific]
g = colorOutEdge(g)
plotNetw(g, "criteria", "TS")


# which(module_tab$CAD_qvalue < 0.1 & module_tab$n_within_mod_endocrines > 0)
# which(p.adjust(module_tab$CAD_pval, method="BH") < 0.1 & module_tab$n_within_mod_endocrines > 0)
# which(apply(angio_pval, 1, min) < 0.01 & module_tab$n_within_mod_endocrines > 0)

# nrow(module_tab)
# sum(!gwas_sig)


# Hierarcical graph
# ------------------------
# tissues
# tissue_cols = brewer.pal(9, "Set1")[-6]



# Tissue node colors
V(g)$color = tissue_cols[
	as.integer(factor(V(g)$primary_tissue, levels=tissues))
]

V(g)$shape = rep("circle", length(V(g)))

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


# plot(lay$layout[, 2], -log10(module_tab$CAD_pval))
# cor.test(lay$layout[, 2], -log10(module_tab$CAD_pval))


# plot(lay$layout[, 2], -log10(module_tab$secreted_protein_pval))
# cor.test(lay$layout[, 2], -log10(module_tab$secreted_protein_pval))


# Hierarcical layout layer
layer = max(lay$layout[, 2]) - lay$layout[, 2] + 1

# cross_tissue

# ct_col = rep("grey", length(cross_tissue))
ct_col = rep(rgb(0.5,0.5,0.5), length(cross_tissue))
ct_col[cross_tissue] = "black"
line_col = rgb(0.2, 0.2, 0.2)

tis_col = tissue_cols[
	as.integer(factor(V(g)$primary_tissue, levels=tissues))
]

pdf("co-expression/eigenNetw/v2/igraph/tissue_legend.pdf")
plot(0, 0, type="n")
legend("topright", legend=tissues, col=tissue_cols, pch=16)
dev.off()

pdf("co-expression/eigenNetw/v2/igraph/h2_features.pdf")
par(mfrow=c(4, 2))

# Dim non-significant associations
col = tis_col
# idx = module_tab$pval_DUKE > 0.05
# col[idx] = addAlpha(col[idx], 0.5)
plot(layer, -log10(module_tab$pval_DUKE),
	col=cross_tissue,  # outline
	bg=col,
	pch=21,
	main="DUKE",
	ylab="-log10 p",
	xlab="Sugiyama layer"
)
abline(h=-log10(0.05), col=line_col, lty=2)


col = tis_col
# idx = module_tab$pval_syntax_score > 0.05
# col[idx] = addAlpha(col[idx], 0.5)
plot(layer, -log10(module_tab$pval_syntax_score),
	col=cross_tissue,  # outline
	bg=col,
	pch=21,
	main="SYNTAX",
	ylab="-log10 p",
	xlab="Sugiyama layer"
)
abline(h=-log10(0.05), col=line_col, lty=2)

col = tis_col
# idx = module_tab$pval_lesions > 0.05
# col[idx] = addAlpha(col[idx], 0.5)
plot(layer, -log10(module_tab$pval_lesions),
	col=cross_tissue,  # outline
	bg=col,
	pch=21,
	main="Lesions",
	ylab="-log10 p",
	xlab="Sugiyama layer"
)
abline(h=-log10(0.05), col=line_col, lty=2)

col = tis_col
# idx = module_tab$pval_ndv > 0.05
# col[idx] = addAlpha(col[idx], 0.5)
plot(layer, -log10(module_tab$pval_ndv),
	col=cross_tissue,  # outline
	bg=col,
	pch=21,
	main="ndv",
	ylab="-log10 p",
	xlab="Sugiyama layer"
)
abline(h=-log10(0.05), col=line_col, lty=2)


# CAD enrichment
x = p.adjust(module_tab$CAD_pval, method="BH")
col = tis_col
# idx = x > 0.1
# col[idx] = addAlpha(col[idx], 0.5)
plot(layer, -log10(x),
	col=cross_tissue,
	bg=col,
	pch=21,
	main="CAD enrichment",
	ylab="-log10 p (BH)",
	xlab="Sugiyama layer"
)
abline(h=-log10(0.1), col=line_col, lty=2)


x = -log10(p.adjust(module_tab$secreted_protein_pval, method="BH"))
x[x>16] = 16
col = tis_col
# idx = x > 0.1
# col[idx] = addAlpha(col[idx], 0.5)
plot(layer, x,
	bg=col,
	col=cross_tissue,
	pch=21,
	main="Secreted protein enrichment ",
	ylab="-log10 p (BH)",
	xlab="Sugiyama layer"
)
abline(h=-log10(0.1), col=line_col, lty=2)


# plot(layer, 1 - module_tab$purity,
# 	col=tissue_cols[
# 		as.integer(factor(V(g)$primary_tissue, levels=tissues))
# 	],
# 	pch=16
# )
# abline(h=0.05, col=line_col, lty=2)


x = -log10(p.adjust(module_tab$pval_BMI, method="BH"))
x[x>16] = 16
col = tis_col
plot(layer, x,
	bg=col,
	col=cross_tissue,
	main="BMI",
	ylab="-log10 p (BH)",
	xlab="Sugiyama layer",
	pch=21
)
abline(h=-log10(0.1), col=line_col, lty=2)


x = -log10(p.adjust(module_tab$pval_LDL, method="BH"))
x[x>16] = 16
plot(layer, x,
	bg=tis_col,
	col=cross_tissue,
	main="LDL",
	ylab="-log10 p (BH)",
	xlab="Sugiyama layer",
	pch=21
)
abline(h=-log10(0.1), col=line_col, lty=2)
dev.off()



cor.test(lay$layout[, 2], -log10(module_tab$pval_DUKE))

# plot(lay$layout[, 2], -log10(module_tab$pval_syntax_score))
# cor.test(lay$layout[, 2], -log10(module_tab$pval_syntax_score))


# plot(degree(g, mode="out"), -log10(module_tab$CAD_pval))
# cor.test(degree(g, mode="out"), -log10(module_tab$CAD_pval))


# cor.test(betweenness(g), -log10(module_tab$CAD_pval))

# plot(betweenness(g), module_tab$purity)

# cor.test(betweenness(g), module_tab$purity)


pdf("co-expression/eigenNetw/v2/igraph/h2.pdf", width=20, height=20)
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
