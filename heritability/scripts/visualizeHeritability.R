rm(list=ls())

library(data.table)
library(RColorBrewer)
library(magicaxis)
library(gplots)
library(magicaxis)


setwd("~/GoogleDrive/projects/STARNET/cross-tissue")

# Load data
mod_tab = fread("co-expression/tables/module_tab.csv")

# heri = fread("heritability/from_johan/STARNETmoduleH2_LD2GWAS302_20180920.csv")  # v1
heri = fread("heritability/from_johan/v2/STARNETmoduleH2_LD2GWAS302_20190318.csv")  # v2
# Note this table has an error in module ID 221. Doesn't affect results.

heri = heri[
	heri[["Module ID"]] != "" & heri[["Module ID"]] != "all",
]
heri[["Module ID"]] = as.integer(heri[["Module ID"]])

heri$CT = heri[["Module ID"]] %in% which(mod_tab$purity < 0.95)


# n_key_drivers = fread("co-expression/annotate/grn_vamsi_eqtl/kd_per_module.csv")
n_key_drivers = fread("co-expression/annotate/grn_vamsi_eqtl/kd_per_module_directed.csv")
colnames(n_key_drivers)[1] = "Module ID"

# Heritability per number of key drivers
heri = merge(heri, n_key_drivers, by="Module ID", all.x=TRUE)

heri$perKD_h2_CAD = heri$module_H2_CAD / heri$n_key_drivers
heri$perKD_h2_CAD[heri$n_key_drivers == 0] = NA  # no divisions by zeros.

head(heri[order(heri$n_key_drivers), ], 60)

# Manual QC of networks with <= 10 KDs
exclude_networks = c(
	13, 217, 218, 78, 162, 167,  # Networks also excluded in tissueHeritabilityTests.R
	30, 166, 58, 68  # blood networks
)
heri_kd = heri[!heri[["Module ID"]] %in% exclude_networks, ]

heri_kd = heri_kd[heri_kd$n_key_drivers != 0, ]



pheno_pval = read.table("pheno/tables/pheno_pval.csv",
	sep=",",
	header=TRUE,
	check.names=FALSE
)


# Define CAD-associated modules. Copied from bayesNet3.R
# ------------------------------------------
features = c("syntax_score", "DUKE", "case_control_DEG")
pheno_padj = matrix(
	p.adjust(
		data.matrix(pheno_pval[, features]),
		method="BH"),
	ncol=3)
colnames(pheno_padj) = c("SYNTAX", "DUKE", "Case-Ctrl DEG")


fdr = 0.01

cad_modules_idx = apply(pheno_padj < fdr,
	1,
	function(x) sum(x) >= 2)

cad_modules = which(cad_modules_idx)
length(cad_modules)


# Modules only
# heri = heri[!heri[["Module ID"]] %in% c("n/a", "", "all"), ]


plotHeriRanking = function(values, modules, mod_sel, sel_col=brewer.pal(9, "Set1")[1], ...) {
	par(xpd=FALSE)
	names(values) = modules

	values = sort(values, decreasing=TRUE)

	# t test
	print(t.test(values[names(values) %in% mod_sel], values[!names(values) %in% mod_sel]))

	# bar_col = rep("white", length(values))
	bar_col = rep(rgb(230, 230, 230, maxColorValue=255), length(values))
	bar_col[names(values) %in% mod_sel] = sel_col

	# border_col = rep("grey", nrow(heri))
	# border_col[heri[["Module ID"]][idx] %in% cad_modules] = brewer.pal(9, "Set1")[1]  # red

	par(lwd=0.3)
	barplot(values,
		col=bar_col,
		border=rgb(100, 100, 100, maxColorValue=255),
		space=0,
		las=2,
		cex.names=0.3,
		...
	)

	par(lwd=1)
	abline(h=0)

	par(xpd=TRUE)
	idx = names(values) %in% mod_sel
	points(which(idx) - 0.5, values[idx],
		cex=1.0,
		pch=21,
		bg=sel_col
	)

}

# pdf("heritability/plots/heritability_ranking.pdf", width=6.0)
pdf("heritability/plots/heritability_ranking_v2.pdf", width=6.0)
# par(mfrow=c(2, 1), lwd=0.3)
par(mfrow=c(2, 1))
plotHeriRanking(heri$module_H2_CAD, heri[["Module ID"]], cad_modules,
	ylab="CAD H2 (%)",
	xlab="Co-expression module")
plotHeriRanking(heri$perSNP_h2_CAD, heri[["Module ID"]], cad_modules,
	ylab="CAD H2 / eQTL SNP (%)",
	xlab=paste0("Co-expression modules (n=", sum(!is.na(heri$perSNP_h2_CAD)), ")")
)

# plotHeriRanking(heri$per_gene_h2_CAD, heri[["Module ID"]], cad_modules,
# 	ylab="CAD H2 / gene (%)",
# 	xlab="Co-expression modules")

dev.off()



pdf("heritability/plots/modSize_heritability_v2.pdf", width=2.9, height=3.4)
magplot(
	heri$genes,
	heri$module_H2_CAD,
	xlab="Module size (genes)",
	ylab="CAD H2 (%)",
	log="xy",
	pch=21,
	bg="grey",
	cex=0.7,
	bty="n")
dev.off()





# Cross-tissue vs tissue-specific comparisons
# -----------------------------------------------------

ct_modules = which(mod_tab$purity < 0.95)


boxplotTest = function(values, colors, ...) {

	t_test = t.test(values[[1]], values[[2]])

	boxplot(values,
		frame=FALSE,
		col=rgb(210, 210, 210, maxColorValue=255),
		border=rgb(150, 150, 150, maxColorValue=255),
		outline=FALSE,
		ylim=range(values, na.rm=TRUE),
		main=paste0("P=", format(t_test$p.value, digits=3)),
		...
	)

	for (i in 1:length(values)) {
		points(
			jitter(rep(i, length(values[[i]])), amount=0.2),
			values[[i]],
			pch=16,
			col=colors[i],
			cex=0.8
		)
	}
}


# pdf("heritability/plots/heritability_cross_tissue.pdf", width=5.5, height=2.5)
# pdf("heritability/plots/heritability_cross_tissue_v2.pdf", width=7.0, height=2.5)
pdf("heritability/plots/heritability_cross_tissue_v3.pdf", width=7.0, height=2.5)  # directed KDA, for # key drivers

par(mfrow=c(1, 5))

# colors = brewer.pal(9, "Set1")
colors = c("black", brewer.pal(9, "Set1"))

values = list(
	TS=heri$module_H2_CAD[!heri[["Module ID"]] %in% ct_modules],
	CT=heri$module_H2_CAD[heri[["Module ID"]] %in% ct_modules]
)
boxplotTest(values, colors, ylab="CAD H2 (%)")

values = list(
	TS=heri$perSNP_h2_CAD[!heri[["Module ID"]] %in% ct_modules],
	CT=heri$perSNP_h2_CAD[heri[["Module ID"]] %in% ct_modules]
)
boxplotTest(values, colors, ylab="CAD H2 / eQTL SNPs (%)")


values = list(
	TS=heri$per_eQTLS_h2_CAD[!heri[["Module ID"]] %in% ct_modules],
	CT=heri$per_eQTLS_h2_CAD[heri[["Module ID"]] %in% ct_modules]
)
boxplotTest(values, colors, ylab="CAD H2 / eQTL (%)")


values = list(
	TS=heri$per_gene_h2_CAD[!heri[["Module ID"]] %in% ct_modules],
	CT=heri$per_gene_h2_CAD[heri[["Module ID"]] %in% ct_modules]
)
boxplotTest(values, colors, ylab="CAD H2 / gene (%)")

# Using QC excluded table
values = list(
	TS=heri_kd$perKD_h2_CAD[!heri_kd[["Module ID"]] %in% ct_modules],
	CT=heri_kd$perKD_h2_CAD[heri_kd[["Module ID"]] %in% ct_modules]
)
boxplotTest(values, colors, ylab="CAD H2 / KD (%)")


dev.off()


plot(heri_kd$perKD_h2_CAD, heri_kd$genes)
cor.test(heri_kd$perKD_h2_CAD, heri_kd$genes)



# Regression analaysis of KD-normalized heritability
# Is there an effecto of CT vs TS after adjusting for module size?
# -----------------------------------------------
fit = lm(perKD_h2_CAD ~ CT, data=heri_kd)
summary(fit)

fit = lm(perKD_h2_CAD ~ CT + genes, data=heri_kd)
summary(fit)


fit = lm(perKD_h2_CAD ~ CT + genes + n_key_drivers, data=heri_kd)
summary(fit)




# Density plots of number of key drivers, TS vs CT
# Modified from qqPlotModules.R
# --------------------------------

densPlot = function(x, y, cex_legend=1.0, ...) {
	densy = density(y)
	densx = density(x)

	wilcox_test = wilcox.test(x, y)

	plot(0,
		type="n",
		xlim=range(densx$x, densy$x),
		ylim=range(densx$y, densy$y),
		bty="l",
		ylab="Density",
		...
	)
	colors = c(
		rgb(228, 26, 28, 150, maxColorValue=255),
		rgb(55, 126, 184, 150, maxColorValue=255)
	)

	polygon(densx, col=colors[1])
	polygon(densy, col=colors[2])

	legend("topright",
		legend=c("Cross-tissue", "Tissue-specific",
			paste0("p=", format(wilcox_test$p.value, digits=4))
		),
		cex=cex_legend,
		fill=c(colors, "white"))
}

# pdf("heritability/plots/n_kd_density.pdf", height=2.9, width=4)
pdf("heritability/plots/n_kd_density_v2.pdf", height=2.9, width=4)  # directed
x = heri_kd$n_key_drivers[heri_kd$CT]
y = heri_kd$n_key_drivers[!heri_kd$CT]
densPlot(
	log10(x),
	log10(y),
	# x,
	# y,
	cex_legend=0.7,
	xlab=expression("Key drivers (log"[10] * " n)"))
dev.off()


