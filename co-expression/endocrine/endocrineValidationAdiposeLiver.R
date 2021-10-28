rm(list=ls())
# Input from endocrineValidationEigenegeneAll.R with selection for all endocrines

endo_val = read.csv("co-expression/tables/endo_adipose_liver_validation.csv")


idx = endo_val$tissue %in% c("VAF", "SF") & endo_val$target_tissue_primary == "LIV"
sum(idx)

length(unique(endo_val$gene_symbol[idx]))


pmat = data.matrix(
	endo_val[idx, c("HMDP_HF_p", "HMDP_chow_p", "HMDP_apoe_p")]
)
rownames(pmat) = endo_val$id[idx]

evaluated = apply(pmat, 1, function(pvals) {
	any(!is.na(pvals))
})
sum(evaluated)


# Multiple hypothesis correction for selected endocrine factors (adipose to liver)
pmat_adjust = apply(pmat, 2, p.adjust)

p_cutoff = 0.05
fdr_cutoff = 0.2

sum(apply(pmat < p_cutoff, 1, any), na.rm=TRUE)
sum(apply(pmat_adjust < fdr_cutoff, 1, any), na.rm=TRUE)

# Print list of endocrine factor IDs. target:gene_symbol
rownames(pmat)[which(apply(pmat < p_cutoff, 1, any))]
rownames(pmat)[which(apply(pmat_adjust < fdr_cutoff, 1, any))]


# Count unique endocrine symbols
# ------------------------------------

length(unique(endo_val$gene_symbol))  # total number of unique endocrine factor gene symbols
length(unique(endo_val$gene_symbol[idx]))  # adipose->liver unique symbols

# Count number of adipose->liver genes with P < 0.05
length(unique(
	sapply(
		strsplit(
			rownames(pmat)[which(apply(pmat < p_cutoff, 1, any))],
			":"
		),
		function(x) x[2]
	)
))


# Test signs of correlation coefficients
# ------------------------------------
plot(endo_val$ts_endo_cor, endo_val$HMDP_HF_cor)
cor.test(endo_val$ts_endo_cor, endo_val$HMDP_HF_cor)

plot(endo_val$ts_endo_cor, endo_val$HMDP_chow_cor)
cor.test(endo_val$ts_endo_cor, endo_val$HMDP_chow_cor)

plot(endo_val$ts_endo_cor, endo_val$HMDP_apoe_cor)
cor.test(endo_val$ts_endo_cor, endo_val$HMDP_apoe_cor)

