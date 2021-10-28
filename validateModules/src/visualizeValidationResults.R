rm(list=ls())

library(data.table)
library(pheatmap)
library(viridis)
library(RColorBrewer)

setwd("~/GoogleDrive/projects/STARNET/cross-tissue")



mod_tab = fread("co-expression/tables/module_tab.csv")
colnames(mod_tab)[1] = "module"

ct_modules = mod_tab$module[mod_tab$purity < 0.95]
ts_modules = mod_tab$module[mod_tab$purity >= 0.95]

# Load NetRep results
blocks = c(1, 2, 3, 4, 5)
d = lapply(blocks, function(i) {
	readRDS(paste0("validateModules/results/network_preservation_block", i, ".rds"))
})


# names(d[[1]])
# sapply(d, function(x) x$p.values)

pvals = rbindlist(
	lapply(d, function(x) {
		tab = x$p.value
		df = data.frame(tab)
		df$module = rownames(tab)
		return(df)
	})
)



# pval_mat = pvals[, -8]
pval_mat = pvals[, c(-2, -5, -7, -8)]
pval_mat = t(pval_mat)
colnames(pval_mat) = pvals$module


pheatmap(
	-log10(pval_mat),
	# color = colorRampPalette(rev(brewer.pal(9, "Spectral")))(100),
	color = colorRampPalette(brewer.pal(9, "YlGnBu"))(20),
	cluster_rows=FALSE,
	fontsize_col=8
)

pval_mat_ts = pval_mat[, colnames(pval_mat) %in% ts_modules]
pval_mat_ct = pval_mat[, colnames(pval_mat) %in% ct_modules]


pdf("validateModules/plots/NetRep_tissue-specific.pdf", height=2, width=12)
pheatmap(
	-log10(pval_mat_ts),
	# color = colorRampPalette(rev(brewer.pal(9, "Spectral")))(100),
	main="Tissue-specific",
	color = colorRampPalette(brewer.pal(9, "YlGnBu"))(20),
	cluster_rows=FALSE,
	border_color="white",
	fontsize_col=8
)
dev.off()

pdf("validateModules/plots/NetRep_cross-tissue.pdf", height=2, width=12)
pheatmap(
	-log10(pval_mat_ct),
	# color = colorRampPalette(rev(brewer.pal(9, "Spectral")))(100),
	main="Cross-tissue",
	color = colorRampPalette(brewer.pal(9, "YlGnBu"))(20),
	cluster_rows=FALSE,
	border_color="white",
	fontsize_col=8
)
dev.off()
