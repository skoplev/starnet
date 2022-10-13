rm(list=ls())

library(data.table)

setwd("~/GoogleDrive/projects/STARNET/cross-tissue")

# Load dbGap to bam table
sample_bam = fread("case-control/dbGaP_submit/data/sample_bam_files.txt")


# Columns in count matrices that are not interpreted as sample column IDs
exclude_cols = -1:-6

count_names = c(
	"gene_counts_AOR.txt",
	"gene_counts_LIV.txt",
	"gene_counts_SF.txt",
	"gene_counts_SKLM.txt",
	"gene_counts_VAF.txt"
)

# Load gene count matrices
gene_counts = lapply(
	count_names,
	function(file) {
		fread(
			paste0("case-control/data/feature_counts/", file)
		)
})


# Change column IDs to assigned dbGap IDs rather than bam file names
for (i in 1:length(gene_counts)) {

	# Get bam file name
	bam_files = basename(colnames(gene_counts[[i]])[exclude_cols])

	# Map bam file to dbGap ID
	idx = match(bam_files, sample_bam$bam_file)
	stopifnot(any(!is.na(idx)))

	# Replace column names in place
	colnames(gene_counts[[i]])[exclude_cols] = sample_bam$SAMPLE_ID[idx]
}
	

# Write modified gene count matrices
for (i in 1:length(gene_counts)) {
	write.table(gene_counts[[i]],
		file=paste0("case-control/data/feature_counts_dbGapID/", count_names[i]),
		row.names=FALSE,
		quote=FALSE,
		sep="\t")
}