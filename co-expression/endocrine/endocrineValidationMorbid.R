rm(list=ls())

library(data.table)
library(RColorBrewer)
library(WGCNA)
library(compiler)
enableJIT(3)

library(GEOquery)
library(Biobase)

library(dplyr)

library(biomaRt)


setwd("~/Google Drive/projects/STARNET/cross-tissue")

data_dir = "~/DataProjects/cross-tissue"  # root of data directory


source("src/permuteTest.R")
source("src/parse.R")


# Get ensembl gene data
message("Loading ensembl")
ensembl = useEnsembl(biomart="ensembl",
	dataset="hsapiens_gene_ensembl")
ensembl_symbol_map = getBM(
	attributes=c("ensembl_gene_id", "ensembl_transcript_id", "hgnc_symbol", "gene_biotype", "entrezgene"),
	mart=ensembl)



# Load STARNET modules
modules = fread("co-expression/tables/modules.csv")
modules$gene_biotype = parseEnsemblBiotype(modules$ensembl)
modules$tissue_biotype = paste(modules$tissue, modules$gene_biotype, sep="_")
modules$ensembl_base = sapply(strsplit(modules$ensembl, "[.]"), function(x) x[1])
modules$tissue_transcript_id = paste(modules$tissue, modules$ensembl_base, sep="_")


# morbid_mat = fread("~/DataBases/MorbidObesity/GSE24335_series_matrix.txt", skip=45, sep="\t")
# morbid_mat = fread("~/DataBases/MorbidObesity/GSE24335_series_matrix.txt", skip=100, sep="\t")
# morbid_mat = fread("~/DataBases/MorbidObesity/GSE24335_series_matrix.txt", skip=45, fill=TRUE)

# morbid_mat = fread("sed -i '$ d' ~/DataBases/MorbidObesity/GSE24335_series_matrix.txt", skip=45)


morbid = getGEO(filename="~/DataBases/MorbidObesity/GSE24335_family.soft.gz", GSEMatrix=TRUE)

# Make sure that platform is identical
gsmplatforms = lapply(GSMList(morbid),function(x) {Meta(x)$platform_id})
table(unlist(gsmplatforms))


platform = getGEO("GPL4372")


# gsmlist = GSMList(morbid)
# Table(gsmlist[[1]])
# Columns(gsmlist[[1]])
# Meta(gsmlist[[1]])

sample_annot = data.frame(
	sample.ID=names(GSMList(morbid)),
	tissue_name = sapply(GSMList(morbid), function(x) Meta(x)$source_name_ch1),
	title = sapply(GSMList(morbid), function(x) Meta(x)$title)
)

sample_annot$MGH.ID = sapply(strsplit(as.character(sample_annot$title), "_"), function(x) x[length(x)])

table(sample_annot$tissue_name)
sample_annot$tissue = recode(sample_annot$tissue_name,
	"omental adipose"="VAF",
	"subcutaneous adipose"="SF",
	"liver"="LIV"
)

table(sample_annot$tissue)
table(table(sample_annot$MGH.ID))

# Gene metadata
# ------------------------------
probeset = Table(GPLList(morbid)[[1]])$ID

# are annotations matched?
stopifnot(all(Table(platform)$ID == probeset))
meta_row = Table(platform)

meta_row$ensembl = ensembl_symbol_map$ensembl_gene_id[match(meta_row$EntrezGeneID, ensembl_symbol_map$entrezgene)]
meta_row$ensembl[is.na(meta_row$EntrezGeneID)] = NA  # ensure that NA values are not matched

meta_row$hgnc_symbol = ensembl_symbol_map$hgnc_symbol[match(meta_row$EntrezGeneID, ensembl_symbol_map$entrezgene)]
meta_row$hgnc_symbol[is.na(meta_row$EntrezGeneID)] = NA  # ensure that NA values are not matched

meta_row$gene_biotype = ensembl_symbol_map$gene_biotype[match(meta_row$EntrezGeneID, ensembl_symbol_map$entrezgene)]



# Get expression matrix for GEO series
# ---------------------------------------
emat = do.call('cbind',
	lapply(GSMList(morbid), function(x) {
		tab = Table(x)
		mymatch = match(probeset, tab$ID_REF)
		return(tab$VALUE[mymatch])
	})
)


# Subset of gene expression matrix with available ensembl IDs for genes
idx = !is.na(meta_row$ensembl)
emat_ensembl = emat[idx, ]
rownames(emat_ensembl) = meta_row$ensembl[idx]
meta_row_ensembl = meta_row[idx, ]


# Reshape data tissue x patient IDs
patient_ids = unique(sample_annot$MGH.ID)

stopifnot(all(patient_ids %in% sample_annot$MGH.ID))

tissues = unique(sample_annot$tissue)
emat_reshape = lapply(tissues, function(tis) {
	# Include
	idx = which(sample_annot$tissue == tis)

	# Match columns to patient ids
	idx = idx[match(patient_ids, sample_annot$MGH.ID[idx])]

	mat = emat_ensembl[, idx]
	colnames(mat) = patient_ids
	# prepend tissue in row names
	rownames(mat) = paste(tis, rownames(mat), sep="_")

	meta_row_tissue = meta_row_ensembl
	meta_row_tissue$tissue = tis

	return(list(mat=mat, meta_row=meta_row_tissue))
})

morbid_reshape = list()
morbid_reshape$meta_row = Reduce(rbind, lapply(emat_reshape, function(x) x$meta_row))
morbid_reshape$mat = Reduce(rbind, lapply(emat_reshape, function(x) x$mat))

morbid_reshape$meta_row$tissue_gene_biotype = paste(morbid_reshape$meta_row$tissue, morbid_reshape$meta_row$gene_biotype, sep="_")

morbid_reshape$mat = t(morbid_reshape$mat)


# Permutation tests 
m = 1000
perm_test = list()

k = 78
mod_idx = modules$clust == k
genes = modules$tissue_transcript_id[mod_idx]
group = modules$tissue_biotype[mod_idx]


# genes[!genes %in% colnames(morbid_reshape$mat)]
# genes[genes %in% colnames(morbid_reshape$mat)]

perm_test[[k]] = corPermTestExprMat(
	expr_mat=morbid_reshape$mat,
	genes=genes,
	m=m,
	mat_group=morbid_reshape$meta_row$tissue_gene_biotype,
	genes_group=group
)


k = 98
mod_idx = modules$clust == k
genes = modules$tissue_transcript_id[mod_idx]
group = modules$tissue_biotype[mod_idx]


# genes[!genes %in% colnames(morbid_reshape$mat)]
# genes[genes %in% colnames(morbid_reshape$mat)]

perm_test[[k]] = corPermTestExprMat(
	expr_mat=morbid_reshape$mat,
	genes=genes,
	m=m,
	mat_group=morbid_reshape$meta_row$tissue_gene_biotype,
	genes_group=group
)


pdf("co-expression/plots/endocrine/validationMorbid/perm_tests/perm_test_78_98.pdf", width=3, height=5)
par(mfrow=c(2, 1))
plotPermuteTest(perm_test[[78]], main="MGH module 78")
plotPermuteTest(perm_test[[98]], main="MGH module 98")
dev.off()



# Cor validation
# ----------------------------------

# source("src/permuteTest.R")
# source("src/parse.R")

# Load STARNET expression data
load(file.path(data_dir, "STARNET/gene_exp_norm_reshape/expr_recast.RData"),
	verbose=TRUE)

starnet = parseExprTable(expr_recast)
rm(expr_recast)

# Rename matrix rows
rownames(starnet$mat) = paste(starnet$meta_row$tissue, starnet$meta_row$ensembl_base, sep="_")


# Morbid obesity daata
morbid_reshape$mat = t(morbid_reshape$mat)
morbid_reshape$meta_row


# k = 98
# k = 78

# k = 123
# k = 197


cor_tests = lapply(1:max(modules$clust), function(k) {
	message(k)

	mod_idx = modules$clust == k
	genes = modules$tissue_transcript_id[mod_idx]
	tissue = modules$tissue[mod_idx]

	out = list()

	# All
	out$all = corModuleTest(mat1=starnet$mat, mat2=morbid_reshape$mat,
		module_genes=genes)

	if (length(table(tissue)) > 1) {
		# Cross-tissue
		out$ct = corModuleTest(mat1=starnet$mat, mat2=morbid_reshape$mat,
			module_genes=genes,
			method="CT")

		for (tis in names(which(table(tissue) > 20))) {
			message(tis)
			genes = modules$tissue_transcript_id[mod_idx & modules$tissue == tis]
			out[[tis]] = corModuleTest(mat1=starnet$mat, mat2=morbid_reshape$mat, module_genes=genes)
		}
	}
	return(out)
})

data_dir = "~/DataProjects/STARNET"
save(cor_tests, file=file.path(data_dir, "moduleVali/valMorbid.RData")


n_genes = sapply(cor_tests, function(x) {
	result = NA
	try({result = x$all$ngenes}, silent=TRUE)
	return(result)
})

r_val = sapply(cor_tests, function(x) {
	result = NA
	try({result = x$all$test$estimate}, silent=TRUE)
	return(result)
})
names(r_val) = 1:length(r_val)


r_val_ct = sapply(cor_tests, function(x) {
	result = NA
	try({result = x$ct$test$estimate}, silent=TRUE)
	if (is.null(result)) {
		result = NA
	}
	return(result)
})
# r_val_ct = unlist(r_val_ct)
names(r_val_ct) = 1:length(r_val_ct)



R2_val = r_val^2
R2_val = R2_val[n_genes > 20 & r_val > 0]

R2_val = sort(R2_val, decreasing=TRUE)
barplot(R2_val, las=2)


R2_val_ct = r_val_ct^2
R2_val_ct = R2_val_ct[n_genes > 20 & r_val_ct > 0]
R2_val_ct = sort(R2_val_ct, decreasing=TRUE)
barplot(R2_val_ct, las=2)



xlab = "STARNET cor."
ylab = "MGH morbid obesity cor."
pdf("co-expression/moduleValidation/morbid/plots/morbid_98_LIV.pdf")
plotCorTest(cor_tests[[98]]$LIV, xlab=xlab, ylab=ylab, main="Module 98 LIV")
dev.off()

pdf("co-expression/moduleValidation/morbid/plots/morbid_78_CT.pdf")
plotCorTest(cor_tests[[78]]$ct, xlab=xlab, ylab=ylab, main="Module 78 cross-tissue")
dev.off()

pdf("co-expression/moduleValidation/morbid/plots/morbid_78_VAF.pdf")
plotCorTest(cor_tests[[78]]$VAF, xlab=xlab, ylab=ylab, main="Module 78 VAF")
dev.off()

pdf("co-expression/moduleValidation/morbid/plots/morbid_78_LIV.pdf")
plotCorTest(cor_tests[[78]]$LIV, xlab=xlab, ylab=ylab, main="Module 78 LIV")
dev.off()

pdf("co-expression/moduleValidation/morbid/plots/morbid_78_SF.pdf")
plotCorTest(cor_tests[[78]]$SF, xlab=xlab, ylab=ylab, main="Module 78 SF")
dev.off()

# hexbinplot(cor_test$ct$r2~cor_test$ct$r1)
# cor_test_ct$test$estimate^2


# morbid_mat_rand = morbid_reshape$mat
# rownames(morbid_mat_rand) = sample(rownames(morbid_mat_rand))
# cor_test_ct_rand = corModuleTest(mat1=starnet$mat, mat2=morbid_mat_rand,
# 	module_genes=genes,
# 	method="CT")
# cor_test_ct_rand$test

# cortest.normal(cor_test_ct$r1, cor_test_ct$r2)
