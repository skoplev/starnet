# Make minimun module annotation table containing cross-tissue and tissue-specific annotation

rm(list=ls())

library(data.table)

setwd("/Users/sk/GoogleDrive/projects/STARNET/cross-tissue")

mod_tab = fread("co-expression/tables/module_tab.csv")
mod_tab = data.frame(mod_tab)

# Add primary tissue column
tissues = c("AOR", "BLOOD", "LIV", "MAM", "SKLM", "SF", "VAF")
mod_tab$primary_tissue = tissues[
	apply(mod_tab[, tissues], 1, which.max)
]

# mod_tab_small = mod_tab[, c(1, 2, 3)]

mod_tab_small = mod_tab[, c("V1", "mod_size", "purity", "primary_tissue")]
colnames(mod_tab_small)[1] = "module"
mod_tab_small$type[mod_tab$purity < 0.95] = "Cross-tissue"
mod_tab_small$type[mod_tab$purity >= 0.95] = "Tissue-specific"


write.csv(mod_tab_small, "co-expression/tables/module_tab_min.csv", row.names=FALSE)

# table(mod_tab_small$type)
