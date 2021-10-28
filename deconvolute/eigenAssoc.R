rm(list=ls())

library(WGCNA)  # for corAndPvalue()
library(pheatmap)
library(viridis)
library(RColorBrewer)
library(data.table)


setwd("~/GoogleDrive/projects/STARNET/cross-tissue")

# Load CIBERSORT
# -------------------
file_names = list.files("deconvolute/out_freq")
cibersort = lapply(file_names, function(file_name) {
	fread(file.path("deconvolute/out_freq", file_name))
})

names(cibersort) = sapply(strsplit(file_names, "[.]"), function(x) x[4])


# Exclude macrophage and foam cell data
cibersort = cibersort[-which(names(cibersort) %in% c("MAC", "FOC", "COR"))]
names(cibersort)

# Load co-expression module eigengenes
# -----------------------------
eig = read.csv("co-expression/tables/eigengene_mat.csv")


# Patient x module
eig_mat = data.matrix(eig[, -1])
rownames(eig_mat) = eig[, 1]
# rownames(eig_mat)

colnames(eig_mat) = paste0("m", 1:ncol(eig_mat))

# eig_mat[1:10, 1:10]



tests = lapply(1:length(cibersort), function(i) {
	# Parse CIBERSORT data frame into sample x cell type matrix
	df = data.frame(cibersort[[i]], check.names=FALSE)
	end = ncol(df)
	cell_freq_mat = df[, c(-1, -end, -(end-1), -(end-2), -(end-3))]
	rownames(cell_freq_mat) = df[, 1]

	head(cell_freq_mat)


	# Match by patient ID to eigengene matrix
	patient_ids = sapply(strsplit(rownames(cell_freq_mat), "_"), function(x) x[2])

	idx = match(rownames(eig_mat), patient_ids)

	cell_freq_mat_matched = cell_freq_mat[idx, ]

	# Prepend tissue to cibersort column names
	colnames(cell_freq_mat_matched) = paste0(names(cibersort)[i], ":", colnames(cell_freq_mat_matched))
	head(cell_freq_mat_matched)

	cor_tests = corAndPvalue(eig_mat, cell_freq_mat_matched)

	# Store matched cell frequency
	cor_tests$cell_freq = cell_freq_mat_matched

	return(cor_tests)
})


# Combine correlation matrix and P values

# cbindlist

pmat = do.call("cbind", lapply(tests, function(x) x$p))

padjmat = matrix(NA, nrow=nrow(pmat), ncol=ncol(pmat))
padjmat[] = p.adjust(pmat)

sum(padjmat < 0.05, na.rm=TRUE)

# mean(apply(padjmat < 0.05, 1, any, na.rm=TRUE))


sort(pmat[98, ])[1:20]
sort(pmat[78, ])[1:20]
sort(pmat[189, ])[1:20]



cmat = do.call("cbind", lapply(tests, function(x) x$cor))


freq_mat = do.call("cbind", lapply(tests, function(x) x$cell_freq))




# Variance explained
R2mat = cmat^2
R2mat[is.na(R2mat)] = 0

rownames(R2mat) = gsub("m", "", rownames(R2mat))

pdf("deconvolute/plots/eigengene_max_var.pdf", height=4, width=7)
max_var_explained = apply(R2mat, 1, max, na.rm=TRUE)
hist(max_var_explained,
	xlab="Max var explained (R2) by CIBERSORT cell type fraction",
	ylab="Co-expression module frequency",
	main="",
	breaks=40, col=brewer.pal(9, "Set1")[2])
dev.off()

# cell_filter = apply(R2mat > 0.5, 2, any)
# cell_filter = unique(apply(R2mat, 2, which.max))


# pdf("deconvolute/plots/eigengene_CIBERSORT_R2.pdf", width=12, height=12)
png("deconvolute/plots/eigengene_CIBERSORT_R2.png", width=16, height=8, units="in",
	res=300  # DPI
)
cell_filter = apply(R2mat > 0.3, 2, any)

pheatmap(
	t(R2mat[, cell_filter]),
	col=viridis(100),
	border_color="white",
	fontsize_row=5,
	fontsize_col=5
)
dev.off()



sort(R2mat[109, ], decreasing=TRUE)[1:10]
sort(R2mat[78, ], decreasing=TRUE)[1:10]
sort(R2mat[28, ], decreasing=TRUE)[1:10]
sort(R2mat[106, ], decreasing=TRUE)[1:10]
sort(R2mat[171, ], decreasing=TRUE)[1:10]
sort(R2mat[150, ], decreasing=TRUE)[1:10]
sort(R2mat[78, ], decreasing=TRUE)[1:10]
sort(R2mat[121, ], decreasing=TRUE)[1:10]
sort(R2mat[167, ], decreasing=TRUE)[1:10]
sort(R2mat[86, ], decreasing=TRUE)[1:10]
sort(R2mat[81, ], decreasing=TRUE)[1:10]
sort(R2mat[135, ], decreasing=TRUE)[1:10]
sort(R2mat[158, ], decreasing=TRUE)[1:10]
sort(R2mat[13, ], decreasing=TRUE)[1:10]
sort(R2mat[98, ], decreasing=TRUE)[1:10]



plotEigCiber = function(j, feature) {
	red = brewer.pal(9, "Set1")[1]

	i = which(colnames(freq_mat) == feature)

	plot(eig_mat[, j], freq_mat[, i],
		xlab=paste0("Eigengene ", j),
		ylab=colnames(freq_mat)[i])

	test = cor.test(eig_mat[, j], freq_mat[, i])

	legend("topleft",
		bty="n",  # no boundary box
		legend=c(
			paste0("r=", format(test$estimate, digits=3)),
			paste0("p=", format(test$p.value, digits=3))
		),
		text.col=red
	)

	fit = lm(freq_mat[, i]~eig_mat[, j])

	abline(fit, col=red)
}


pdf("deconvolute/plots/eigengene_CIBERSORT_examples.pdf", width=16, height=3)
par(mfrow=c(1, 6))
plotEigCiber(150, "AOR:blood:macrophage")
plotEigCiber(78, "VAF:omentum:preadipocyte")
plotEigCiber(98, "LIV:retina:retinal pigment epithelial cell")
plotEigCiber(121, "VAF:blood:B cell")
plotEigCiber(158, "MAM:blood:neutrophil")
plotEigCiber(13, "MAM:bone marrow:stem cell")
dev.off()
