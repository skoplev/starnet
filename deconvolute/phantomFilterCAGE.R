# Filter phantom CAGE gene expression data.
rm(list=ls())

library(data.table)
library(readxls)
library(compiler)
enableJIT(3)

library(org.Hs.eg.db)
library(annotate)

data_dir = "/Users/sk/DataProjects/cross-tissue"
setwd("/Users/sk/Google Drive/projects/cross-tissue")


# Returns list with 
# Reads FANTOM annotated gene expression data
# Returns list of mrna matrix and annot table
parseFANTOM = function(file_path) {

	out = list()
	# Load human CAGE data
	out$cage = fread(
		paste("grep -v ^##",
			file_path
		)
	)

	out$annot = out$cage[, 1:7, with=FALSE]
	out$cage = out$cage[, 8:ncol(out$cage), with=FALSE]


	# Get entrez id,
	# first number in format entrez:<#>,entrez<#>
	out$annot$entrezgene_id_single = sapply(
		strsplit(out$annot$entrezgene_id, ","), function(strings) {
			strsplit(strings[1], ":")[[1]][2]
		}
	)

	# Remove peaks with multiple gene associations
	multiple_gene_associations = which(
			sapply(strsplit(out$annot$entrezgene_id, ","),
			length
		) > 1
	)

	out$annot$entrezgene_id_single[multiple_gene_associations] = NA
	# empty strings
	out$annot$entrezgene_id_single[
		out$annot$entrezgene_id_single == ""
		] = NA

	return(out)
}

# Collapses fantom expression matrix by $annot$entrezgene_id_single
# Returns matrix of summed values with rows as unique entrez ids
collapseFANTOM = function(fantom) {
	entrez_id = unique(fantom$annot$entrezgene_id_single)
	entrez_id = entrez_id[!is.na(entrez_id)]

	single_ids = fantom$annot$entrezgene_id_single

	# Index row numbers
	message("Indexing entrez id rows")
	row_index = sapply(entrez_id, function(id) {
		return(which(id == single_ids))
	}, USE.NAMES=TRUE)

	cage_mat = as.matrix(fantom$cage)  # matrix instead of 

	message("Aggregating expression matrix by genes")
	mat = sapply(entrez_id, function(id) {

		# subcage = fantom$cage[single_ids == id, drop=FALSE]
		# subcage = fantom$cage[row_index[[id]], drop=FALSE]
		subcage_mat = cage_mat[row_index[[id]], , drop=FALSE]

		if (nrow(subcage_mat) == 1) {
			# No aggregation
			return(subcage_mat)
		} else {
			# Multiple rows, return column-wise sum
			return(
				colSums(subcage_mat, na.rm=TRUE)
			)
		}
		return(NA)
	})

	mat = t(mat)
	colnames(mat) = colnames(cage_mat)

	return(mat)
}


# Returns matched sample annot table
# fantom is a list that contains $cage columns storing the sample ids in its column names, as
# the last "."-separated item.
matchSampleAnnot = function(fantom, sample_annot) {
	# Get sample id of
	sample_id = sapply(
		strsplit(colnames(human$cage), "[.]"),  # split by "."
		function(vec) vec[length(vec)]  # return last element
	)

	# Match columns of CAGE expression matrix to sample annotation
	sample_annot_matched = sample_annot[match(sample_id, sample_annot[["Source Name"]]),]

	return(sample_annot_matched)
}

filterCols = function(fantom, cols) {
	out = list()
	out$sample_annot = fantom$sample_annot[cols,]
	out$cage = fantom$cage[, cols, with=FALSE]
	out$annot = fantom$annot
	return(out)
}



# Loads FANTOM sample annotation from xlsx file,
# renames data entries
loadSampleAnnotFANTOM = function(file_path) {
	# Load sample annot
	sample_annot = read_excel(file_path, sheet=1)

	# from->to rename
	relabel_cell_types = c(
		"b cell"="B cell",
		"t cell"="T cell",
		"t cell, CD8+"="T cell, CD8+",
		"t cell, CD4+"="T cell, CD4+",
		"t cell, nk, immature"="T cell, nk, immature",
		"t cell, gamma-delta"="T cell, gamma-delta"
	)

	col_name = "Characteristics [Cell type]"  # column to rename in sample_annot
	for (i in 1:length(relabel_cell_types)) {
		from = names(relabel_cell_types[i])
		to = relabel_cell_types[i]

		sample_annot[[col_name]] =
			replace(sample_annot[[col_name]],
				sample_annot[[col_name]] == from,
				to)

	}
	return(sample_annot)
}

# Load FANTOM5 CAGE gene expression data
human = parseFANTOM(
	file.path(data_dir, "FANTOM5/hg19.cage_peak_phase1and2combined_tpm_ann.osc.txt")
)

# Load human sample annotation data
sample_annot = loadSampleAnnotFANTOM(
	file.path(data_dir, "FANTOM5/HumanSamples2.0.sdrf.xlsx"))
# table(sample_annot[["Characteristics [Cell type]"]])

# Match annotation table to sample order
human$sample_annot = matchSampleAnnot(human, sample_annot)

# Filter expression matrix and sample annotation
cols = which(
	human$sample_annot[["Characteristics [Category]"]] == "primary cells"
)

human_filter = filterCols(human, cols)


# Make tissue:cell type ID
tissue_cell_type_id = paste(
	human_filter$sample_annot[["Characteristics[Tissue]"]],
	human_filter$sample_annot[["Characteristics [Cell type]"]], sep=":"
)


# Collapse genes from alternative slice sites by summing values
human_mat = collapseFANTOM(human_filter)


# Annotate and filter gene expression matrix
# ------------------------------------------------------------------------
# Annotate matrix columns with tissue:cell type id
colnames(human_mat) = tissue_cell_type_id

# Reorder matrix based on tissue:cell type id.
# WARNING: order no longer agrees with human_filter$sample_annot
human_mat = human_mat[,order(colnames(human_mat))]

# Remove entries with non-sensical IDs
remove_entries = c(
	"NA:NA",
	"ANATOMICAL SYSTEM:CELL MIXTURE - tissue sample",
	"ANATOMICAL SYSTEM:NA",
	"blood:CELL MIXTURE - tissue sample",
	"Buffy coat:NA"
)

human_mat = human_mat[,!colnames(human_mat) %in% remove_entries]

# Rename genes using gene symbols,
rownames(human_mat) = getSYMBOL(rownames(human_mat), data="org.Hs.eg")
# mapIds(org.Hs.eg.db, keys=rownames(human_mat), column="ENSEMBL", keytype="ENTREZID")

message("Removing CAGE peaks with no HUGO names: ", sum(is.na(rownames(human_mat))))
human_mat = human_mat[!is.na(rownames(human_mat)), ]

# Filter out genes with low expression across samples
epsilon = 0.0001  # small number
low_expression = apply(human_mat, 1, sum) < epsilon

message("low expression genes removed: ", sum(low_expression))

# Filter genes that are not detected
human_mat = human_mat[!low_expression,]


# Write table for use with
write.table(human_mat, file.path(data_dir, "FANTOM5/cell_type_basis/primary_cells.tsv"),
	sep="\t",
	quote=FALSE, col.names=NA)

