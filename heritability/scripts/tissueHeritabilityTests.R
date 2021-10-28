rm(list=ls())

library(data.table)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
# library(gridExtra)
library(cowplot)

setwd("~/GoogleDrive/projects/STARNET/cross-tissue")

heri = fread("heritability/from_johan/tissue/Report_hsq_res_snpmod224_20190318_LZ_tissues20200107.csv")

n_key_drivers = fread("co-expression/annotate/grn_vamsi_eqtl/kd_per_module.csv")
colnames(n_key_drivers)[1] = "Module ID"


# Remove emtpy rows and columns
heri = heri[
	heri[["Module ID"]] != "" & heri[["Module ID"]] != "All",
	1:26]
heri[["Module ID"]] = as.integer(heri[["Module ID"]])

# Heritability per number of key drivers
heri = merge(heri, n_key_drivers, by="Module ID", all.x=TRUE)




heri$perKD_h2_CAD = heri$pct_h2_CAD / heri$n_key_drivers
heri$perKD_h2_CAD[heri$n_key_drivers == 0] = NA  # no divisions by zeros.
# heri$perKD_h2_CAD[heri$n_key_drivers <= 1] = NA
# heri$perKD_h2_CAD[heri$n_key_drivers <= 2] = NA


# heri$perKD_h2_CAD_log = log10(heri$perKD_h2_CAD + 1)
# heri[order(heri$perKD_h2_CAD, decreasing=T), ]



# Investigate modules with low numbers of key drivers
head(heri[order(heri$n_key_drivers), ], 60)

# Exclude BLOOD modules from subsequent comparisons
heri = heri[heri$primary_tissue != "BLOOD", ]
dim(heri)

# Set tissue colors
# heri$primary_tissue = factor(heri$primary_tissue, levels=c("AOR", "MAM", "LIV", "SKLM", "VAF", "SF", "BLOOD"))
heri$primary_tissue = factor(heri$primary_tissue, levels=c("AOR", "MAM", "LIV", "SKLM", "VAF", "SF"))

tissue_colors = brewer.pal(9, "Set1")[c(1, 5, 7, 3, 4, 8, 2)]



# QC for key drivers, manual assessment for networks <= 10 key drivers
exclude_networks = c(13, 217, 218, 78, 162, 167)
heri_kd = heri[!heri[["Module ID"]] %in% exclude_networks, ]
dim(heri)



# Syntactically correct names
colnames(heri) = make.names(colnames(heri))


# Interpret as numeric values
heri$perSNP_h2_CAD..ld02. = as.numeric(heri$perSNP_h2_CAD..ld02.)
heri$perSNP_h2_CAD..ld05. = as.numeric(heri$perSNP_h2_CAD..ld05.)
heri$Module.H2.CAD.gene = as.numeric(heri$Module.H2.CAD.gene)

# tissue_color = brewer.pal(9, "Set1")[1:7]



# -----------------


# feature = "pct_h2_CAD"
# feature = "perSNP_h2_CAD..ld02."
# feature = "perSNP_h2_CAD..ld05."
# feature = "Module.H2.CAD.gene"

# Returns ggplot object
plotHeri = function(heri_in, subset, feature, max_value=NA) {
	heri_sub = subset(heri_in, type==subset)

	var.equal = T  # t-test

	if (is.na(max_value)) {
		max_value = max(heri_in[[feature]], na.rm=TRUE)
	}

	# Pairwise comparisons
	all_comparisons = as.list(data.frame(combn(
		length(levels(heri_sub$primary_tissue)),  # 6 (or 7 including BLOOD)
		2)))

	# Test 
	all_pairwise_tests = lapply(all_comparisons, function(pair) {
		idx_x = as.integer(heri_sub$primary_tissue) == pair[1]
		idx_y = as.integer(heri_sub$primary_tissue) == pair[2]

		t.test(heri_sub[[feature]][idx_x], heri_sub[[feature]][idx_y], var.equal=var.equal)
	})

	test_pvals = sapply(all_pairwise_tests, function(x) x$p.value)

	sig_comparisons = all_comparisons[test_pvals < 0.05]

	message(subset)
	message(feature)
	message(sig_comparisons)

	p = ggplot(
		heri_sub,
		aes(x=primary_tissue,
			y=!!sym(feature),  # string to symbol
			color=primary_tissue)) +
		geom_boxplot(outlier.shape=NA) +
		coord_cartesian(ylim=c(0, max_value)) +  # same y-scale for tissue- and cross-tissue
		ggtitle(subset) +
		geom_jitter(shape=16, position=position_jitter(0.2)) +
		stat_compare_means(
			# comparisons=all_comparisons,
			comparisons=sig_comparisons,
			method="t.test", method.args=list(var.equal=var.equal)
		) +
		scale_color_manual(values=tissue_colors) + 
		theme_classic() +
		theme(axis.text.x=element_text(angle=90, hjust=1))  # vertical tissue labels

	return(p)
}



# pdf("heritability/plots/heritability_by_module.pdf")
# pdf("heritability/plots/heritability_by_module_v2.pdf", height=9)
# pdf("heritability/plots/heritability_by_module_v2.pdf")
# pdf("heritability/plots/heritability_by_module_v3.pdf")
# pdf("heritability/plots/heritability_by_module_v4.pdf")
pdf("heritability/plots/heritability_by_module_v5.pdf", width=6.5)

plot_grid(
	plotHeri(heri_in=heri, subset="Cross-tissue", feature="pct_h2_CAD"),
	plotHeri(heri_in=heri, subset="Tissue-specific", feature="pct_h2_CAD", max_value=7.5),  # 5
	plotHeri(heri_in=heri, subset="Cross-tissue", feature="perSNP_h2_CAD..ld02."),
	plotHeri(heri_in=heri, subset="Tissue-specific", feature="perSNP_h2_CAD..ld02."),
	# plotHeri(heri_in=heri, subset="Cross-tissue", feature="Module.H2.CAD.gene"),
	# plotHeri(heri_in=heri, subset="Tissue-specific", feature="Module.H2.CAD.gene"),
	plotHeri(heri_in=heri_kd, subset="Cross-tissue", feature="perKD_h2_CAD"),
	plotHeri(heri_in=heri_kd, subset="Tissue-specific", feature="perKD_h2_CAD", max_value=0.2),  # 0.4

	# plotHeri(heri=heri, subset="Cross-tissue", feature="perKD_h2_CAD_log"),
	# plotHeri(heri=heri, subset="Tissue-specific", feature="perKD_h2_CAD_log"),

	ncol=2,
	align="v")
dev.off()


t.test(
	heri$perKD_h2_CAD[heri$type == "Tissue-specific" & heri$primary_tissue == "AOR"],
	heri$perKD_h2_CAD[heri$type == "Tissue-specific" & heri$primary_tissue == "SKLM"],
	var.equal=T
)
