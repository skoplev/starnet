rm(list=ls())

library(data.table)
library(stringr)

setwd("~/GoogleDrive/projects/STARNET/cross-tissue")

# Load STARNET phenotype data
pheno = fread("~/GoogleDrive/projects/STARNET/phenotype/data/current/STARNET_main_phenotype_table.2017_12_03.tsv")


# Load case-control sample IDs
# -------------------------------------------

# Get bam file names, indirectly from gene count matrices
bam_files = lapply(
	c(
		"gene_counts_AOR.txt",
		"gene_counts_LIV.txt",
		"gene_counts_SF.txt",
		"gene_counts_SKLM.txt",
		"gene_counts_VAF.txt"
	),
	function(file) {
		mat = fread(
			paste0("case-control/data/feature_counts/", file)
		)
		return(colnames(mat)[-1:-6])  # strip non-file name columns

})

# Combine from each tissue
bam_files = unlist(bam_files)

# strip file path
bam_files = basename(bam_files)

# Create samples table from bam file catalogue
# -----------------------------------------
# Table mutable throughout rest of script
samples = data.frame(
	bam_file=bam_files,
	stringsAsFactors=FALSE)

samples$starnet_patient_id = sapply(strsplit(samples$bam_file, "_"), function(x) x[3])
samples$tissue = sapply(strsplit(samples$bam_file, "_"), function(x) x[2])


# Load ID map and add new patient IDs assigned with deidentified ID 'B'[0-9]
# -------------------
id_map = fread("case-control/dbGaP_submit/from_oscar_original_submit/dbgap_from_of/deidentified_subjects.txt")
colnames(id_map) = c("starnet_patient_id", "dbgap_patient_id")


# Identify new patients added and assign new random ID
new_patient_ids = samples$starnet_patient_id[!samples$starnet_patient_id %in% id_map$starnet_patient_id]
new_patient_ids = unique(new_patient_ids)


set.seed(42)  # for reproducibility (and consistency) of deidentified patient IDs
new_id_map = data.frame(
	starnet_patient_id=new_patient_ids, 
	dbgap_patient_id=paste0("B", sample(1:1000, size=length(new_patient_ids))),
	stringsAsFactors=FALSE
)

# Append to previous table
new_id_map = rbind(id_map, new_id_map)

write.table(new_id_map, "case-control/dbGaP_submit/data/deidentified_subjects.txt",
	sep="\t",
	row.names=FALSE,
	quote=FALSE
	)
rm(id_map)

# Shorthand named vectors for mapping patient IDs
map = list()
map$starnet_to_dbgap = new_id_map$dbgap_patient_id
names(map$starnet_to_dbgap) = new_id_map$starnet_patient_id


# Annotate samples table with assigned patient IDs
samples$dbgap_patient_id = map$starnet_to_dbgap[
	samples$starnet_patient_id
]


# Write new subject consent table, including previous submission patient IDs.
# Order is maintained from id_map table
# -------------------------------------

# Load previous subject consent table subject samples mapping
# subject_consent = fread("case-control/dbGaP_submit/from_oscar_original_submit/dbgap_from_of/subject_consent_data.txt")

subject_consent = data.frame(
	SUBJECT_ID=new_id_map$dbgap_patient_id,
	CONSENT=2,  # Health/Medical/Biomedical, all original entries are assigned '2'
	SEX=pheno$Sex[
		match(new_id_map$starnet_patient_id, pheno$starnet.ID)
	],
	stringsAsFactors=FALSE
)

# Void gender codes
subject_consent$SEX = replace(subject_consent$SEX,
	subject_consent$SEX == "",
	"NULL")
subject_consent$SEX = replace(subject_consent$SEX,
	subject_consent$SEX == "patient nr male7female",
	"NULL")

# Capitalize gender in line with dbGaP requirement
subject_consent$SEX = stringr::str_to_title(subject_consent$SEX)


write.table(subject_consent, "case-control/dbGaP_submit/data/subject_consent_data.txt",
	sep="\t",
	row.names=FALSE,
	quote=FALSE)


# Subject phenotype table.
# Table was not inncluded in STARNET v1 submission, hence all patient IDs included in v2 dbGaP submission
# ----------------------------------------
stopifnot(new_id_map$dbgap_patient_id == subject_consent$SUBJECT_ID)  # can use new_id_map for STARNET patient IDs

subject_phenotype = subject_consent[, c("SUBJECT_ID", "SEX")]  # without consent

subject_phenotype$AGE = pheno$Age[match(new_id_map$starnet_patient_id, pheno$starnet.ID)]

subject_phenotype$CASE_CTRL = pheno$CAD.status[match(new_id_map$starnet_patient_id, pheno$starnet.ID)]

# Capitalize first letter
subject_phenotype$CASE_CTRL = stringr::str_to_title(subject_phenotype$CASE_CTRL)

write.table(subject_phenotype, "case-control/dbGaP_submit/data/subject_phenotype_data.txt",
	sep="\t",
	row.names=FALSE,
	quote=FALSE)



# Assign unique sample IDs to new samples
# ----------------------

# Load previous subject sample mapping
subject_sample_map = fread("case-control/dbGaP_submit/from_oscar_original_submit/dbgap_from_of/subject_mapping_data.txt")

new_subject_sample_map = data.frame(
	SUBJECT_ID=samples$dbgap_patient_id,
	SAMPLE_ID=paste0(samples$tissue, "_", samples$dbgap_patient_id),  # proposed, will be made unique if duplicate samples between v1 and v2
	SAMPLE_USE="Seq_RNA_Expression",
	stringsAsFactors=FALSE
)

# Append to previous map
new_subject_sample_map = rbind(subject_sample_map, new_subject_sample_map)

# Assign unique sample names to resequenced samples
new_subject_sample_map$SAMPLE_ID = make.unique(new_subject_sample_map$SAMPLE_ID)

write.table(new_subject_sample_map, "case-control/dbGaP_submit/data/subject_mapping_data.txt",
	sep="\t",
	row.names=FALSE,
	quote=FALSE)


# Annotate samples table with assigned unique sample IDs
new_sample_ids = new_subject_sample_map$SAMPLE_ID[-1:-nrow(subject_sample_map)]

samples$SAMPLE_ID = new_sample_ids


# Make sample attribute table for new samples
# --------------------------------------
# sample attribute table, including new 'case-control' variable indication.
# Design='Case-control'/'Cases'
# 

tissue_codes = c(
	AOR="Aortic wall",
	Blood="Blood",
	LIV="Liver",
	MAM="Internal mammary artery",
	SF="Subcutaneous fat",
	SKLM="Skeletal muscle",
	VAF="Visceral fat"
)


sample_attribute = data.frame(
	SAMPLE_ID=samples$SAMPLE_ID,
	BODY_SITE=tissue_codes[samples$tissue],
	BODY_SITE_ABBREV=samples$tissue,
	ANALYTE_TYPE="RNA",
	DESIGN="Case-control",
	stringsAsFactors=FALSE)


write.table(sample_attribute, "case-control/dbGaP_submit/data/sample_attributes_data.txt",
	sep="\t",
	row.names=FALSE,
	quote=FALSE)


# Write sample bam file info, not directly required for dbGap submission.
# -------------------------
write.table(samples[, c("SAMPLE_ID", "bam_file")], "case-control/dbGaP_submit/data/sample_bam_files.txt",
	sep="\t",
	row.names=FALSE,
	quote=FALSE)



# Manual checks for mapping errors
pheno[pheno$starnet.ID == "c59", ]
pheno[pheno$starnet.ID == "313", ]

