# Functions for parsing gene expression matrices.

# Test that the gene names are the same across matrices
testExprRownames = function(mats) {
	for (i in 2:length(mats)) {
		if (!all(rownames(mats[[1]]) == rownames(mats[[i]]))) {
			stop("Gene symbol mismatch")
		}
	}
}


# Loads normalized gene expression matrices from folder.
# Filters genes that varies less than provided standard deviation.
# Returns single data frame with data
loadNormData = function(data_dir, min_sd=0, exclude_files) {
	require(data.table)
	require(reshape2)
	require(plyr)

	# Load expression data
	expr_files = list.files(
		data_dir,
		pattern="*.mat")

	expr_files = expr_files[!expr_files %in% exclude_files]

	# Load all expression data from each tissue
	expr_mats = sapply(expr_files, function(file_name) {
		message("Loading ", file_name)

		data_table = fread(
			file.path(data_dir, file_name),
			sep="\t"
		)

		mat = data_table[, 2:ncol(data_table), with=FALSE]
		mat = as.matrix(mat)
		rownames(mat) = data_table[[1]]  # first column

		if (ncol(mat) > nrow(mat)) {
			warning("Gene expression matrix has more columns, than rows.")
		}

		return(mat)
	})

	# testExprRownames(expr_mats)

	# Prefilter genes
	if (min_sd > 0) {
		message("Filtering data")
		expr_mats = sapply(expr_mats, function(mat) {
			# Gene-wise standard deviation
			std_devs = apply(mat, 1, sd)

			return(mat[std_devs > min_sd, ])
		})
	}

	# Reshape filtered expression matrices
	message("Melting data")
	expr_melted = lapply(expr_mats, function(mat) {
		melted = melt(t(mat))  # data frame of melted expression matrix, efficient melt

		melted = rename(melted, c(
			"Var1"="sample_id",
			"Var2"="transcript_id"))

		# Get sample info
		tissue = sapply(
			strsplit(colnames(mat), "_"),
			function(x) x[1]  # 1st elem
		)
		patient_id = sapply(
			strsplit(colnames(mat), "_"),
			function(x) x[2]  # 2nd elem
		)

		# Annotate sample
		# Get sample number of melted data,
		sample_idx = match(melted[,1], colnames(mat))

		melted = data.frame(
			tissue=tissue[sample_idx],
			patient_id=patient_id[sample_idx],
			melted)

		return(melted)
	})
	names(expr_melted) = names(expr_mats)

	# Combine melted data into single data frame
	expr_melted_all = rbindlist(expr_melted)

	# Recast into table with patients in columns and tissue-specific transcripts in rows
	message("Recasting data")
	expr_recast = dcast(expr_melted_all, tissue + transcript_id ~ patient_id)

	# clean other elements from memory
	rm(expr_melted_all, expr_melted, expr_mats, scaled_mat)
	gc()  # garbage collection

	return(expr_recast)
}

# Parses transcript ids from meta_genes table.
# Changes ensembl and gene_symbol column if input data frame and returns new data frame.
# Supports gene symbols with "_".
parseTranscriptId = function(meta_genes) {
	# Get ensembl IDs as last "_" separated string
	meta_genes$ensembl = sapply(
		strsplit(as.character(meta_genes$transcript_id), "_"),
		function(x) {
			x[length(x)]  # last element
	})

	meta_genes$gene_symbol = sapply(
		strsplit(as.character(meta_genes$transcript_id), "_"),
		function(x) {
			if (length(x) == 2) {
				# Gene symbols does not contain "_"
				return(x[1])
			} else {
				# Gene symbol contains "_"
				gene_symbol = paste(x[-length(x)], collapse="_")  # not last element
				# message(gene_symbol)
				return(gene_symbol)
			}
	})
	return(meta_genes)
}


# Parses module data
parseModuleData = function(mod_env)  {
	# Clusters as integers
	mod_env$clust = as.integer(factor(mod_env$bwnet$colors))
	mod_env$meta_genes = parseTranscriptId(mod_env$meta_genes)

	# Rename eigengene matrix
	eigen_gene_names = substring(colnames(mod_env$bwnet$eigengenes), 3)
	eigen_gene_n = match(eigen_gene_names, levels(factor(mod_env$bwnet$colors)))

	colnames(mod_env$bwnet$eigengenes) = eigen_gene_n

	mod_env$bwnet$eigengenes = mod_env$bwnet$eigengenes[, order(eigen_gene_n)]

	return(mod_env)
}



parseCibersortFiles = function(freq_files, data_dir) {
	require(data.table)
	ciber_freq = lapply(freq_files, function(file_name) {
		file_path = file.path(data_dir, "CIBERSORT/out_freq", file_name)

		tissue = strsplit(file_name, "[.]")[[1]][4]

		# Load CIBERSORT data
		freq = fread(file_path)

		# Parse header separately
		header = read.table(file_path, nrows=1, sep="\t")

		header = unlist(lapply(header, as.character))
		header = c("sample", header)

		# header
		colnames(freq) = header

		# Only numerical fraction entries
		freq_mat = freq[, 2:(ncol(freq) - 3)]
		freq_mat = data.matrix(freq_mat)

		rownames(freq_mat) = freq$sample
		colnames(freq_mat) = paste(tissue, colnames(freq_mat), sep=":")

		return(freq_mat)
	})

	return(ciber_freq)
}

# Returns vector of gene symbols reported for GWAS trait.
# gwas is a data.table from EBI GWAS catalog
getGWAS = function(gwas, trait) {
	exclude_genes = c("intergenic")

	idx = gwas[["DISEASE/TRAIT"]] == trait
	message("Found GWAS: ", sum(idx))

	# Find associated genes
	gwas_genes = trimws(
		unlist(
			strsplit(gwas[["REPORTED GENE(S)"]][idx], ",")
		)
	)
	gwas_genes = unique(gwas_genes)

	gwas_genes = gwas_genes[!gwas_genes %in% exclude_genes]

	return(gwas_genes)
}


getCADGenes = function(data_dir) {
	# Load Deloukas (2013) genes associated with CAD
	delou = read.table(file.path(data_dir, "GWAS/Deloukas/ng.csv"),
		skip=1,
		sep=",",
		header=TRUE)

	# Nikpay GWAS table
	nikpay = read.table(file.path(data_dir, "GWAS/Nikpay/ng/Suppl Table 4-Table 1.csv"),
		skip=3,
		sep=",",
		header=TRUE)
	nikpay$Locus.name

	# Howson GWAS table
	howson = read.table(file.path(data_dir, "GWAS/Howson/ng.csv"),
		skip=3,
		sep=",",
		header=TRUE)
	howson = howson[1:25, ]  # trim table

	# Remove notes
	howson$Nearest.gene = gsub("\\(nsSNP\\)", "", howson$Nearest.gene)


	# Combine CAD nearest loci genes
	# Get CAD gene symbols from proximal loci
	cad_genes = c(
		as.character(delou$Loci_Nearest_Transcript),
		as.character(nikpay$Locus.name),
		as.character(howson$Nearest.gene)
	)

	cad_genes = paste(cad_genes, collapse="/")  # single string
	cad_genes = strsplit(cad_genes, "[/,-]")[[1]]  # separate by / , or -
	cad_genes = trimws(cad_genes)
	cad_genes = unique(cad_genes)
	cad_genes = cad_genes[cad_genes != ""]

	return(cad_genes)
}
