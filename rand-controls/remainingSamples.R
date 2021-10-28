# Tissue-by-tissue statistics for the selected samples
rm(list=ls())

library(data.table)
library(RColorBrewer)
library(readxl)

# Parse data frame of sample IDs in the format: 
parseSamples = function(samples) {
	colnames(samples)[1] = "sample_id"
	samples$sample_id = as.character(samples$sample_id)

	samples$tissue = sapply(
		strsplit(samples$sample_id, "_"),
		function(x) x[1])
	samples$ctrl = sapply(
		strsplit(samples$sample_id, "_"),
		# compare first character
		function(x) substr(x[2], 1, 1) == "c")

	samples$patient_id = sapply(
		strsplit(samples$sample_id, "_"),
		function(x) x[2])
	return(samples)
}

# Compare to sets of numbers using two boxplots and t test.
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


setwd("/Users/sk/Google Drive/projects/STARNET-controls")


# Load selected samples
samples = read.table("rand_samples/rand_within_tissue.txt",
	header=FALSE)

samples = parseSamples(samples)


# Load missing samples due to low volume or RNA concentration. Recieved from Payal
missing_samples = read.table("tissue_avail/low_volume_missing/missing.txt",
	header=FALSE)
missing_samples = as.character(missing_samples[, 1])


# Test if all missing sample IDs are found
if (!all(missing_samples %in% samples$sample_id)) {
	stop("Missing samples not found in selected sample.")
}

# Filter out missing samples
samples = samples[!samples$sample_id %in% missing_samples, ]


# Load phenotype data on both cases
pheno = fread(file.path(
	"/Volumes/SANDY/phenotype_data",
	"STARNET_main_phenotype_table.cases.Feb_29_2016.tbl"
))

cases_info = data.frame(
	id=pheno$starnet.ID,
	bmi=pheno$BMI,
	age=pheno$Age,
	gender=pheno$Sex
)


# Load phenotype data for controls
control_info = read_excel("control_info/c1 to c250.xlsx")

ctrls_info = data.frame(
	id=control_info[["starnet-ID"]],
	age=control_info$Age,
	bmi=control_info[["BMI(kg/m2)"]],
	gender=tolower(control_info$Sex)
)

# Combined annotation table for both controls and samples
all_info = rbind(cases_info, ctrls_info)

info_matched = all_info[match(samples$patient_id, all_info$id), ]


# Barplots of sample counts 
pdf("tissue_avail/plots/tissue_avail_bars_total.pdf", width=8, height=3.5)
par(mfrow=c(1, 2))

cmat = rbind(
	table(samples$tissue[samples$ctrl]),
	table(samples$tissue[!samples$ctrl]))

barplot(cmat,
	beside=TRUE,  # grouped barplot 
	main=paste0("Total samples (n=", nrow(samples), ")"),
	ylab="Patients",
	xlab="Tissue"
)


patient_tissue_n_ctrls = table(table(samples$patient_id[samples$ctrl]))
patient_tissue_n_samples = table(table(samples$patient_id[!samples$ctrl]))

# Check if count catagories match
if (!(all(names(patient_tissue_n_ctrls) == names(patient_tissue_n_samples)))) {
	stop("Count position mismatch")
}

cmat = rbind(rev(patient_tissue_n_ctrls), rev(patient_tissue_n_samples))

barplot(cmat,
	beside=TRUE,  # grouped barplot 
	main=paste0("Total samples (n=", nrow(samples), ")"),
	ylab="Patients",
	xlab="Tissues"
)

legend("topright", legend=c("Ctrls", "Cases"),
	pch=22, 
	# cex=0.8,
	pt.bg=gray.colors(2))

dev.off()


# Tissue order
tissues = levels(factor(samples$tissue))

# Balance test for each tissue the for BMI, age, and gender
pdf("distributions/age_bmi_final.pdf", height=5)
par(bty="n", mfcol=c(2, length(tissues)))

cols = c(gray.colors(1), brewer.pal(9, "Set1")[2])
cols_subtle = c(gray.colors(2)[2], brewer.pal(9, "Pastel1")[2])

for (tissue in tissues) {

	idx_samples = samples$tissue == tissue & !samples$ctrl
	idx_ctrls = samples$tissue == tissue & samples$ctrl

	boxplotPair(
		list(ctrls=info_matched$age[idx_ctrls], cases=info_matched$age[idx_samples]),
		cols, cols_subtle, ylab=paste("Age", tissue)
	)

	boxplotPair(
		list(ctrls=info_matched$bmi[idx_ctrls], cases=info_matched$bmi[idx_samples]),
		cols, cols_subtle, ylab=paste("BMI", tissue)
	)
}
dev.off()

cols = c(gray.colors(1), brewer.pal(9, "Set1")[2])
cols_subtle = c(gray.colors(2)[2], brewer.pal(9, "Pastel1")[2])


pdf("distributions/gender_by_tissue_final.pdf", width=2.5, height=8)
par(mfrow=c(length(tissues), 2),
	mar=c(4, 4, 2, 2))
for (tissue in tissues) {
	idx_samples = samples$tissue == tissue & !samples$ctrl
	idx_ctrls = samples$tissue == tissue & samples$ctrl

	barplot(table(info_matched$gender[idx_ctrls]), las=2, ylab="Individuals", main=tissue)
	barplot(table(info_matched$gender[idx_samples]), las=2, col=cols[2], main=tissue)
}
dev.off()

# gender_ratio = sapply(tissues, function(tissue) {
# 	idx_samples = samples$tissue == tissue & !samples$ctrl
# 	idx_ctrls = samples$tissue == tissue & samples$ctrl

# 	ctrl_tab = table(info_matched$gender[idx_ctrls])
# 	sample_tab = table(info_matched$gender[idx_samples])

# 	c(ctrl_tab[1] / ctrl_tab[2], sample_tab[1] / sample_tab[2])
# })

# plot(gender_ratio[1, ], gender_ratio[2, ],
# 	xlim=c(0.75, 1.25), ylim=c(0.75, 1.25))

# abline(0, 1)
