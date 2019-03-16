rm(list=ls())
library(RColorBrewer)
library(data.table)

data_dir = "/Users/sk/DataProjects/cross-tissue"  # root of data directory
setwd("/Users/sk/Google Drive/projects/cross-tissue")

# Load module overview table
module_tab = read.table("co-expression/tables/module_tab.csv",
	sep=",",
	header=TRUE
)

endocrine_tab = read.table("co-expression/tables/endocrine_tab.csv",
	sep=",",
	header=TRUE
)


# Ordered tissue composition based on correlation p-values 
# -----------------------------------------------------
n = 30
# tissues = c("AOR", "BLOOD", "LIV", "MAM", "SKLM", "SF", "VAF")

tissues = c("AOR", "BLOOD", "SKLM", "VAF", "MAM", "LIV",  "SF")
# tissue_cols = brewer.pal(9, "Set1")[-6]


# feature = "cad_qvalue"
# feature = "pval_DUKE"
# feature = "pval_syntax_score"
# feature = "pval_LDL"

features = c(
	# "cad_qvalue",
	"CAD_pval",
	"pval_DUKE",
	"pval_syntax_score",
	"pval_LDL")

# par(mfcol=c(2, 1), mar=c(3, 4, 3, 2))
pdf("co-expression/annotate/plots/tissue_composition2.pdf",
	width=14, height=4)
par(mfcol=c(2, length(features)), mar=c(1, 4, 3, 2))

for (feature in features) {
	module_tab_order = module_tab[order(module_tab[feature]), ]
	count_mat = t(as.matrix(module_tab_order[1:n, tissues]))
	barplot(count_mat, col=brewer.pal(9, "Set1")[-6],
		ylim=c(0, min(2000, max(count_mat))),
		ylab="Module size",
		main=feature,
		border=NA,
		las=2
	)

	logp = -log10(module_tab_order[1:n, feature])
	names(logp) = rownames(module_tab_order)[1:n]

	# barplot(logp,
	# 	ylab=paste("-log10 ", feature),
	# 	las=2,
	# 	ylim=c(0, max(logp)))

	plot(logp,
		ylab=paste("-log10 ", feature),
		bty="n",
		xaxt="n",
		pch=16,
		ylim=c(0, max(logp)))

	abline(h=-log10(0.1), lty=2,
		col=brewer.pal(9, "Set1")[1])
}

legend("topright",
	legend=tissues,
	col=brewer.pal(9, "Set1")[-6],
	pch=15
)
dev.off()



# Within-module endocrine factors of CAD and DUKE modules
include_modules = which(module_tab$cad_pval < 0.05 | module_tab$pval_DUKE < 0.05)

# include_modules = which(module_tab$cad_pval < 0.05 | module_tab$pval_DUKE < 0.05 | module_tab$pval_syntax_score < 0.05)

# sub_endocrine_tab = endocrine_tab[endocrine_tab$endocrine_in_module == TRUE, ]

sub_endocrine_tab = endocrine_tab[endocrine_tab$endocrine_in_module == TRUE & endocrine_tab$module %in% include_modules, ]

write.table(sub_endocrine_tab,
	"co-expression/annotate/plots/endocrine_DUKE_CAD_modules.csv",
	row.names=FALSE,
	sep=",")

# CAD genes that are also eQTLs



# barplot(
# 	sort(table(sub_endocrine_tab$module), decreasing=TRUE),
# 	main="CAD and DUKE-associated modules (p < 0.1)",
# 	ylab="Endocrine factors",
# 	xlab="Module",
# 	las=2)
