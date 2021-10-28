rm(list=ls())

library(readxl)
library(RColorBrewer)

setwd("/Users/sk/Google Drive/projects/STARNET-controls")

# Load phenotype data on both cases and controls
pheno = fread(file.path(
	"/Volumes/SANDY/phenotype_data",
	"STARNET_main_phenotype_table.cases.Feb_29_2016.tbl"
))

# Load additional clinical information on.
# Recieved from Raili Ermel
control_info = read_excel("control_info/c1 to c250.xlsx")
# control_info = as.data.frame(control_info)

# STARNET population
cases_info = data.frame(
	id=pheno$starnet.ID,
	bmi=pheno$BMI,
	age=pheno$Age,
	gender=pheno$Sex
)
# Exclude STARNET IDs 1-100. Assumes complete and ordered
cases_info = cases_info[101:nrow(cases_info), ]

# Controls
ctrls_info = data.frame(
	id=control_info[["starnet-ID"]],
	age=control_info$Age,
	bmi=control_info[["BMI(kg/m2)"]],
	gender=tolower(control_info$Sex)
)


# Load tissue availability data
# ---------------------------------------------------------------
avail = list()
for (i in 1:5) {
	avail[[i]] = read_excel("tissue_avail/Tissue sample list-Simon.xlsx",
		sheet=i,
		skip=1)
	avail[[i]] = as.data.frame(avail[[i]])
}
names(avail) = c("LIV", "SF", "AOR", "SKLM", "VAF")

# Parse tissue availability
# Paste ctrl encoding
avail = lapply(avail, function(x) {
	# x[, 1] = paste0("c", x[, 1])
	x[[1]] = paste0("c", x[[1]])
	return(x)
})	

# Load matched patient IDs
rand_patient_ids = read.table("rand_ids/rand_sample_ctrls.txt")
patient_ids = read.table("rand_ids/sample_ctrls.txt")  # non-randomized

ctrls = sapply(patient_ids[,1], function(s) substr(s, 1, 1) == "c")  # bool
patient_ids[!ctrls, 1]

# Primary overlap counts
overlap_counts = data.frame(tissue=names(avail))
ctrls_found = list()  # samples 
samples_found = list()
for (i in 1:length(avail)) {
	ctrls_found[[i]] = intersect(avail[[i]][,1], patient_ids[ctrls, 1])
	overlap_counts$ctrl_found[i] = length(ctrls_found[[i]])

	samples_found[[i]] = intersect(avail[[i]][, 3], patient_ids[!ctrls, 1])
	overlap_counts$samples_found[[i]] = length(samples_found[[i]])
}
names(ctrls_found) = names(avail)
names(samples_found) = names(avail)

par(mfrow=c(2, 1))
barplot(overlap_counts$ctrl_found,
	names.arg=overlap_counts$tissue,
	las=2,
	main="Samples")
barplot(overlap_counts$samples_found,
	names.arg=overlap_counts$tissue,
	las=2,
	main="Ctrls")

cmat = rbind(overlap_counts$ctrl_found, overlap_counts$samples_found)
colnames(cmat) = overlap_counts$tissue


# Barplots of counts
pdf("tissue_avail/plots/tissue_avail_bars.pdf", width=4, height=6)
par(mfrow=c(2, 1))
# cols = c(
# 	rgb(0.4, 0.4, 0.4),  # ctrl
# 	brewer.pal(9, "Pastel1")[1]  # samples
# )

barplot(cmat,
	# col=cols,
	beside=TRUE,  # grouped barplot 
	main="Found samples (out of 250 and 248)",
	ylab="Matched samples"
)

# Patient tissue completion of selection

tissue_complete_ctrls = rev(table(table(unlist(ctrls_found))))
tissue_complete_samples = rev(table(table(unlist(samples_found))))

cmat_tissue = rbind(tissue_complete_ctrls, tissue_complete_samples)

barplot(cmat_tissue, 
	beside=TRUE,
	# col=cols,
	ylab="Matched individuals",
	main="Tissue completion",
	xlab="Tissues"
)

legend("topright", legend=c("Ctrl", "STARNET"),
	pch=22, 
	# cex=0.8,
	pt.bg=gray.colors(2))
dev.off()



# Ctrl IDs with at least on tissue available
all_ctrls = levels(factor(unlist(ctrls_found)))
ctrls_avail = sapply(ctrls_found, function(patient_ids) {
	return(all_ctrls %in% patient_ids)
})
rownames(ctrls_avail) = all_ctrls

# Sample IDs with at least one tissue available
all_samples = levels(factor(unlist(samples_found)))
samples_avail = sapply(samples_found, function(patient_ids) {
	return(all_samples %in% patient_ids)
})
rownames(samples_avail) = all_samples


# Find, count, and list patients with incomplete tissue coverage.
# ----------------------------------------------------------------
min_tissue_coverage = 3
# Count
sum(apply(samples_avail, 1, sum) >= min_tissue_coverage)
sum(apply(ctrls_avail, 1, sum) >= min_tissue_coverage)

# Exclusion lists for 1st tier samples to be sequenced
exclude_samples = rownames(samples_avail)[
	apply(samples_avail, 1, sum) < min_tissue_coverage]

exclude_ctrls = rownames(ctrls_avail)[
	apply(ctrls_avail, 1, sum) < min_tissue_coverage]

# Make missing sample IDs
ctrls_missing = lapply(ctrls_found, function(ids) {
	all_ctrls[!all_ctrls %in% ids]
})

samples_missing = lapply(samples_found, function(ids) {
	all_samples[!all_samples %in% ids]
})

missing = mapply(c, samples_missing, ctrls_missing)
# Prepend tissue
for (i in 1:length(missing)) {
	missing[[i]] = paste0(names(missing)[i], "_", missing[[i]])
}

# unlist(missing)
write.table(unlist(missing),
	file="missing/suggested_samples_ctrls.txt",
	row.names=FALSE, col.names=FALSE, quote=FALSE)


# Coverage plots
cols = c(gray.colors(2)[2], brewer.pal(9, "Pastel1")[2])
pdf("tissue_avail/plots/avail_ctrls_per_indi.pdf", height=1.5)
# pdf("tissue_avail/plots/avail_ctrls_per_indi.pdf", width=12, height=20)
heatmap.2(t(ctrls_avail * 1),  # bool->numeric
	col=c("white", cols[1]),
	trace="none",
	cexCol=0.25,
	cexRow=0.4,
	main="Ctrl",
	xlab="Individual",
	dendrogram="none",
	cellnote=t(ctrls_avail * 1),
	notecex=0.2,
	notecol="black"
)
dev.off()

pdf("tissue_avail/plots/avail_samples_per_indi.pdf", height=1.5)
heatmap.2(t(samples_avail * 1),
	# col=c("white", "grey"),
	col=c("white", cols[2]),
	trace="none",
	cexCol=0.25,
	cexRow=0.4,
	main="STARNET",
	xlab="Individual",
	dendrogram="none",
	cellnote=t(samples_avail * 1),
	notecex=0.2,
	notecol="black"
)
dev.off()



# Phenotype distributions of included 
final_samples = unique(unlist(samples_found))
final_samples = setdiff(final_samples, exclude_samples)

final_ctrls = unique(unlist(ctrls_found))
final_ctrls = setdiff(final_ctrls, exclude_ctrls)


final_ctrls_info = ctrls_info[ctrls_info$id %in% final_ctrls, ]
final_cases_info = cases_info[cases_info$id %in% final_samples, ]


# plot(density(final_ctrls_info$age, na.rm=TRUE))
# lines(density(final_cases_info$age, na.rm=TRUE), col="red")


boxplotPair = function(dlist, cols, cols_subtle, ...) {
	t_test = t.test(dlist[[1]], dlist[[2]])

	boxplot(dlist,
		outline=FALSE,
		ylim=range(dlist, na.rm=TRUE),
		# ylab="Age (years)",
		col=cols_subtle,
		main=paste0("p=", format(t_test$p.value, digits=3)),
		las=2,
		...
	)

	for (i in 1:length(dlist)) {
		points(
			jitter(rep(i, length(dlist[[i]])), amount=0.2),
			pch=16,
			col=cols[i],
			cex=0.5,
			dlist[[i]])
	}
}

cols = c(gray.colors(1), brewer.pal(9, "Set1")[2])
cols_subtle = c(gray.colors(2)[2], brewer.pal(9, "Pastel1")[2])

pdf("distributions/age_bmi_boxplot_post_exclusion.pdf", width=3.0, height=8)
par(mfrow=c(2, 1), 	bty="n")
dlist = list(
	Ctrl=final_ctrls_info$age,
	Cases=final_cases_info$age
)
boxplotPair(dlist, cols, cols_subtle, ylab="Age (years)")


dlist = list(
	Ctrl=final_ctrls_info$bmi,
	Cases=final_cases_info$bmi
)

boxplotPair(dlist, cols, cols_subtle, ylab="BMI")
dev.off()

pdf("distributions/gender_post_exclusion.pdf", width=3.5, height=3)
par(mfrow=c(1, 2))
barplot(table(final_ctrls_info$gender), las=2, ylab="Individuals")
barplot(table(final_cases_info$gender), las=2, col=cols[2])
dev.off()


# Randomize samples across and between tissues
# -------------------------------------------------------------------

# Combine controls and samples
found = mapply(c, samples_found, ctrls_found)

# Exclude controls and samples below tissue coverage threshold
found = lapply(found, function(ids) {
	ids[!ids %in% c(exclude_ctrls, exclude_samples)]
})

# Prepend tissue indicator to each patient
for (i in 1:length(found)) {
	found[[i]] = paste0(names(avail)[i], "_", found[[i]])
}

# Randomize controls and samples. Either with tissues kept seperate or mixed.
found_rand_tissue = unlist(lapply(found, sample))
found_rand_mixed = sample(unlist(found))
write.table(found_rand_tissue,
	file="rand_samples/rand_within_tissue.txt",
	row.names=FALSE,
	col.names=FALSE,
	quote=FALSE)

write.table(found_rand_mixed,
	file="rand_samples/rand_mixed_tissue.txt",
	row.names=FALSE,
	col.names=FALSE,
	quote=FALSE)




