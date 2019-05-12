rm(list=ls())

library(data.table)
library(WGCNA)
library(RColorBrewer)
library(gplots)

setwd("~/GoogleDrive/projects/STARNET/cross-tissue")

source("src/permuteTest.R")
source("src/parse.R")

data_dir = "/Users/sk/DataProjects/cross-tissue"  # root of data directory


# Load and parse HMDP gene expression data
hmdp = lapply(list.files(file.path(data_dir, "HMDP/chow_HF")),
	function(file_name) {
	file_path = file.path(data_dir, "HMDP/chow_HF", file_name)
	d = fread(file_path)
	d = data.frame(d)


	rownames(d) = d[, 1]  # mouse strain names
	d = data.matrix(d[, -1])

	genes = read.table(file_path, nrows=1)
	genes = t(genes)[, 1]
	genes = unname(genes)

	colnames(d) = genes

	return(d)
})
names(hmdp) = list.files(file.path(data_dir, "HMDP/chow_HF"))

sapply(hmdp, dim)


# Select all of FDR threshold endocrine candidates
# endocrine_sel = "all"
endocrine_sel = "FDR"

# Load endocrine table
if (endocrine_sel == "FDR") {
	endo = fread("co-expression/tables/CT_endocrines_TS_interactions.csv")
} else if (endocrine_sel == "all") {
	endo = fread("co-expression/tables/CT_endocrines_TS_interactions_all.csv")
} else {
	stop("Wrong selector: ", endocrine_sel)
}
# sum(endo$tissue %in% c("SF", "VAF") & endo$target_tissue_primary == "LIV")

# endo[endo$clust == 78 & endo$target_clust == 98, ]
# endo_all[endo_all$clust == 78 & endo_all$target_clust == 98, ]


# Load STARNET modules
modules = fread("co-expression/tables/modules.csv")




# Load mouse homology and map STARNET gene symbols to mouse gene symbols
human_mouse_homology = fread(file.path(data_dir, "MGI/HOM_MouseHumanSequence.rpt"))

# Separate table into human and mouse
human_homology = human_mouse_homology[
	human_mouse_homology[["Common Organism Name"]] == "human", ]

mouse_homology = human_mouse_homology[
	human_mouse_homology[["Common Organism Name"]] == "mouse, laboratory", ]


modules$mouse_symbol = findMouseHomologue(modules$gene_symbol, human_homology, mouse_homology)

endo$mouse_symbol = findMouseHomologue(endo$gene_symbol, human_homology, mouse_homology)



# Calculate all correlation matrices, for mouse genes with homologues in STARNET
cor_mats = lapply(hmdp, matchCor, unique(modules$mouse_symbol))



# Calculate selected cross-tissue correlation matrices. Adipose-liver. HMDP.
# -----------------------------------------------------

liver_adipose_cmats = list()

# Test sample ID consistency
if (!all(rownames(hmdp$Liver_chow_male) == rownames(hmdp$Adipose_chow_male))) {
	stop("HMDP ID mismatch")
}

if (!all(rownames(hmdp$Liver_HF_male) == rownames(hmdp$Adipose_HF_male))) {
	stop("HMDP ID mismatch")
}



# Adipose-liver chow cross-tissue correlation matrix
x1 = hmdp$Liver_chow_male
colnames(x1) = paste("liver", colnames(x1), sep=":")

x2 = hmdp$Adipose_chow_male
colnames(x2) = paste("adipose", colnames(x2), sep=":")

liver_adipose_cmats$chow_male = matchCor(
	cbind(x1, x2),  # combined expression matrix
	allowed=c(
		paste("adipose", unique(modules$mouse_symbol), sep=":"),
		paste("liver", unique(modules$mouse_symbol), sep=":")
	)
)


# Adipose-liver HF correlation
x1 = hmdp$Liver_HF_male
colnames(x1) = paste("liver", colnames(x1), sep=":")

x2 = hmdp$Adipose_HF_male
colnames(x2) = paste("adipose", colnames(x2), sep=":")

liver_adipose_cmats$HF_male = matchCor(
	cbind(x1, x2),  # combined expression matrix
	allowed=c(
		paste("adipose", unique(modules$mouse_symbol), sep=":"),
		paste("liver", unique(modules$mouse_symbol), sep=":")
	)
)


# Make definitions (gene lists) of mouse homologue genes for adipose-liver modules.
# Symbols are prepended with 'adipose:' or 'liver:'
# -------------------------------------------------------------------
ct_modules = c(78)

ct_modules_mouse_symbols = lapply(ct_modules, function(k) {
	liv_symbols = modules$mouse_symbol[modules$clust == k & modules$tissue == "LIV"]
	liv_symbols = unique(na.omit(liv_symbols))
	liv_symbols = paste("liver", liv_symbols, sep=":")

	adipose_symbols = modules$mouse_symbol[modules$clust == k & (modules$tissue == "VAF" | modules$tissue == "SF")]
	adipose_symbols = unique(na.omit(adipose_symbols))
	adipose_symbols = paste("adipose", adipose_symbols, sep=":")
	return(c(liv_symbols, adipose_symbols))
})
names(ct_modules_mouse_symbols) = ct_modules



# 



# permutation tests
m = 10000

ct_perm_tests = lapply(liver_adipose_cmats, function(cmat) {
	mod_tests = lapply(ct_modules_mouse_symbols, function(genes) {
		corPermTest(cmat, genes, m)
	})
	return(mod_tests)
})


# genes = ct_modules_mouse_symbols[["78"]]
# cmat = liver_adipose_cmats$chow_male

# cmat_sub = cmat[rownames(cmat) %in% genes, colnames(cmat) %in% genes]

# rand_sample = sample(nrow(cmat), size=length(genes))
# cmat_rand = cmat[rand_sample, rand_sample]

# x = abs(cmat_sub[lower.tri(cmat_sub)])
# y = abs(cmat_rand[lower.tri(cmat_rand)])

# quantile(x, seq(0, 1, 0.05))
# quantile(y, seq(0, 1, 0.05))


# rand_sample = sample(nrow(cmat), size=length(genes))
# cmat_rand = cmat[rand_sample, rand_sample]
# x = abs(cmat_rand[lower.tri(cmat_rand)])
# # y = abs(cmat_rand[lower.tri(cmat_rand)])

# wilcox.test(x, y)
# t.test(x, y)


# mean(abs(cmat_sub[lower.tri(cmat_sub)]))
# median(abs(cmat_sub[lower.tri(cmat_sub)]))

# heatmap.2(cmat_sub,
# 	trace="none",
# 	col=colorRampPalette(rev(brewer.pal(9, "RdBu")))(100)
# )


# Validate tissue-specific liver modules. Get mouse gene symbols for liver modules targeted by endocrine candidates
# -------------------------------------------------------------------
liver_modules = endo$target_clust[endo$target_tissue == "LIV"]
liver_modules = unique(liver_modules)

tissue = "LIV"
liver_modules_mouse_symbols = lapply(liver_modules, function(k) {
	modules$mouse_symbol[modules$clust == k & modules$tissue == tissue]
})
names(liver_modules_mouse_symbols) = liver_modules

liver_modules_mouse_symbols = lapply(liver_modules_mouse_symbols, na.omit)



m = 10000

# Test per tissue
perm_tests = lapply(cor_mats, function(cmat) {
	# Test per module
	mod_tests = lapply(liver_modules_mouse_symbols, function(genes) {
		corPermTest(cmat, genes, m)
	})
	return(mod_tests)
})





k = which(names(perm_tests) == "Liver_HF_male")
# k = which(names(perm_tests) == "Adipose_HF_male")
mice = names(perm_tests)[k]

i = 8  # mod 98
# i = 4  # mod 96
mod = names(perm_tests[[k]])[i]
pdf(paste0("co-expression/plots/endocrine/validationHMDP/perm_tests/perm_test_", mice, "_", mod, ".pdf"),
	width=5,
	height=3.5)

# hist(perm_tests[[k]][[i]]$stat_perm,
# 	xlim=c(0, 0.5),
# 	col="black",
# 	xlab=expression("Mean connectivity " * bar(c)),
# 	main=paste0(mice, ", module ", mod),
# 	breaks=50)
# points(perm_tests[[k]][[i]]$stat, 0, col=brewer.pal(9, "Set1")[1], pch=16)

plotPermuteTest(perm_tests[[k]][[i]],
	main=paste0(mice, ", module ", mod))

# plotPermuteTest(ct_perm_tests$HF_male[["78"]])

dev.off()



# Same as above, but including cross-tissue module 78
k = which(names(perm_tests) == "Liver_HF_male")
mice = names(perm_tests)[k]

i = 8  # mod 98
mod = names(perm_tests[[k]])[i]

pdf(paste0("co-expression/plots/endocrine/validationHMDP/perm_tests/perm_test_module78_module98.pdf"),
	width=5,
	height=6.5)

par(mfrow=c(2, 1))
plotPermuteTest(ct_perm_tests$HF_male[["78"]],
	main=paste0("Adipose-liver HF_male, CT module 78"))
plotPermuteTest(perm_tests[[k]][[i]],
	main=paste0(mice, ", module ", mod))
legend("topright", pch=15, col="black", legend="Null")
dev.off()




mean_connect = sapply(perm_tests, function(perm) makePermTable(perm)$mean_connect)
p_vals = sapply(perm_tests, function(perm) makePermTable(perm)$p)
rownames(p_vals) = makePermTable(perm_tests$Liver_HF_male)[, 1]


pdf("co-expression/plots/endocrine/validationHMDP/liver_targets_mean_connect_all_tissues.pdf", width=6, height=5)
cols = brewer.pal(9, "Set1")[c(1, 2, 3, 5, 6)]
sel = c(8, 4, 5, 18)
pch_pts = c(21, 22, 23, 24)

par(mar=c(10, 4, 2, 2))
boxplot(mean_connect, las=2, pch=16, ylab=expression("Mean connectivity " * bar(c)),
	border=rgb(0.4, 0.4, 0.4)
)

for (i in 1:length(sel)) {
	lines(mean_connect[sel[i], ], col=cols[i], lwd=1.5)
	points(mean_connect[sel[i], ], col=cols[i], pch=pch_pts[i], bg="white", lwd=1.5)
	# points
}

legend("topright",
	legend=makePermTable(perm_tests$Liver_HF_male)[, 1][sel], 
	pch=pch_pts,
	# pch=16,
	lwd=1.5,
	pt.bg="white",
	col=cols)
dev.off()


pdf("co-expression/plots/endocrine/validationHMDP/liver_targets_perm_test_heatmap.pdf")
heatmap.2(-log10(p_vals), trace="none",
	col=colorRampPalette(brewer.pal(9, "YlOrRd"))(20),
	key.title="",
	key.xlab=expression("-log"[10] *" p"),
	mar=c(12, 16),
	ylab="Liver module",
	main="HMDP connectivity"
)
dev.off()



# i = 8
par(mfrow=c(length(perm_tests_liver_HF), 1), mar=c(1, 1, 1, 1))

for (i in 1:length(perm_tests_liver_HF)) {
	hist(perm_tests_liver_HF[[i]]$stat_perm,
		xlim=c(0, 0.5),
		col="black",
		main=names(perm_tests_liver_HF)[i],
		breaks=50)

	points(perm_tests_liver_HF[[i]]$stat, 0, col=brewer.pal(9, "Set1")[1], pch=16)
}



# Quantify expression of endocrine factors
# -----------------------------------------------------------
# Parse expression matrix
emat = expr_recast[, 3:ncol(expr_recast)]
emat = data.matrix(emat)

meta_genes = expr_recast[, 1:2]
meta_genes = as.data.frame(meta_genes)
meta_genes = parseTranscriptId(meta_genes)


# Reverse pseudo log2 counts
mean_counts = apply(2^emat + 1, 1, mean, na.rm=TRUE)

for (k in 1:nrow(endo)) {
	tissue = endo$tissue[k]
	gene_symbol = endo$gene_symbol[k]

	idx = which(meta_genes$gene_symbol == gene_symbol & meta_genes$tissue == tissue)

	# meta_genes[idx, ]

	if (length(idx) > 1) {
		warning("Multiple gene symbols found")
	}

	endo$mean_norm_counts[k] = mean(mean_counts[idx])
}


# Write HMDP-validated endocrine factors
val_idx = endo$HMDP_HF_ts_endo_cor_p < 0.05 | endo$HMDP_chow_ts_endo_cor_p < 0.05
# val_idx = p.adjust(endo$HMDP_HF_ts_endo_cor_p) < 0.2 | p.adjust(endo$HMDP_chow_ts_endo_cor_p) < 0.2
# val_idx = qvalue(na.omit(endo$HMDP_HF_ts_endo_cor_p)$qvalue < 0.1 | qvalue(endo$HMDP_chow_ts_endo_cor_p)$qvalue < 0.1
sum(val_idx, na.rm=TRUE)

write.csv(endo[val_idx, ], "co-expression/tables/CT_endocrines_TS_interactions_mouse.csv", row.names=FALSE)

length(unique(endo[val_idx, ]$gene_symbol))
length(unique(endo[val_idx, ]$gene_symbol))

sub_endo = endo[val_idx, ]

length(unique(sub_endo$gene_symbol[sub_endo$supernetwork_edge]))

mean(endo$supernetwork_edge)
mean(sub_endo$supernetwork_edge)

sub_endo$target


# hist(log10(endo[val_idx, ]$mean_norm_counts), breaks=50)
# head(endo, 50)


# # gene_symbol = "Lep"
# gene_symbol = "Il6"

# # k = which(colnames(hmdp$Adipose_HF_male) == gene_symbol)
# hist(hmdp$Adipose_HF_male[, k], breaks=40, col="black")
# hist(hmdp$Adipose_chow_male[, k], breaks=40, col="black")


# x = sort(apply(hmdp$Adipose_chow_male, 2, mean), decreasing=TRUE)
x = sort(apply(hmdp$Adipose_HF_male, 2, mean), decreasing=TRUE)

cols = brewer.pal(9, "Set1")

plot(x, type="l")


gene_symbol = "Il6"
k = which(names(x) == gene_symbol)
points(k, x[k], col=cols[1], pch=16)

gene_symbol = "Lep"
k = which(names(x) == gene_symbol)
points(k, x[k], col=cols[2], pch=16)



# hist(endocrine_eigen_cor_adipose_liv_chow$p, breaks=50)
# hist(endo$HMDP_chow_ts_endo_cor_p, breaks=50)

# plot(density(endocrine_eigen_cor_adipose_liv_chow$p, na.rm=TRUE))

# sum(p.adjust(endo$HMDP_HF_ts_endo_cor_p) < 0.1, na.rm=TRUE)

# sum(endo$HMDP_HF_ts_endo_cor_p < 0.05, na.rm=TRUE)



endo[endo$target_clust == 98, ]
endo[endo$target_clust == 18, ]
endo[endo$target_clust == 41, ]


head(endo[order(endo$HMDP_ts_endo_cor_p), ], 100)

# endo[endo$clust == 189, ]


# plot(abs(endo$ts_endo_cor), abs(endo$HMDP_HF_ts_endo_cor))
# cor.test(abs(endo$ts_endo_cor), abs(endo$HMDP_HF_ts_endo_cor))




# Test endocrine-eigengene associations
gene = "Lep"

gene = "Fstl3"
gene = "Lbp"
gene = "Stc2"
mod = 98

gene = "Adamts4"
mod = 181

gene = "S100a9"
gene = "Camp"
gene = "Mmp8"
mod = 18


gene = "Nenf"
mod = 126

gene = "Gzma"
mod = 96


# HF
x = hmdp$Adipose_HF_male[, match(gene, colnames(hmdp$Adipose_HF_male))]
y = eigen_liver_HF$eigengenes[, colnames(eigen_liver_HF$eigengenes) == mod]

plot(x, y)
cor.test(x, y)


# chow
x = hmdp$Adipose_chow_male[, match(gene, colnames(hmdp$Adipose_chow_male))]
y = eigen_liver_chow$eigengenes[, colnames(eigen_liver_chow$eigengenes) == mod]

plot(x, y)
cor.test(x, y)
