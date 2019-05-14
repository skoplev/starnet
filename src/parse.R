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

# Loads ensembl biotypes by loading Ensembl database
parseEnsemblBiotype = function(ensembl) {
	require(biomaRt)
	ensembl_base = sapply(strsplit(ensembl, "[.]"), function(x) x[1])

	# Get ensembl gene data
	message("Loading ensembl")
	
	ensembl = useEnsembl(biomart="ensembl",
		dataset="hsapiens_gene_ensembl")
	ensembl_symbol_map = getBM(
		attributes=c("ensembl_gene_id", "ensembl_transcript_id", "hgnc_symbol", "gene_biotype"),
		mart=ensembl)

	# Map ensembl IDs to biotype
	gene_biotype = ensembl_symbol_map$gene_biotype[
		match(ensembl_base, ensembl_symbol_map$ensembl_gene_id)
	]

	return(gene_biotype)
}

# Parses transcript ids from meta_genes table.
# Changes ensembl and gene_symbol column if input data frame and returns new data frame.
# Supports gene symbols with "_".
parseTranscriptId = function(meta_genes) {
	require(biomaRt)

	message("Parsing transcript IDs")
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

	meta_genes$tissue_transcript_id = paste(meta_genes$tissue, meta_genes$transcript_id, sep="_")

	meta_genes$gene_biotype = parseEnsemblBiotype(meta_genes$ensembl)

	return(meta_genes)
}


# Parses combined expression data  for all STARNET tissues
parseExprTable = function(expr_recast) {
	# Parse expression matrix
	mat = data.matrix(expr_recast[, 3:ncol(expr_recast)])

	meta_row = expr_recast[, 1:2]
	meta_row = parseTranscriptId(meta_row)

	rownames(mat) = meta_row$tissue_transcript_id

	# Tissue-biotype groups
	meta_row$tissue_biotype = paste(
		meta_row$tissue,
		meta_row$gene_biotype,
		sep="_"
	)

	meta_row$ensembl_base = sapply(strsplit(meta_row$ensembl, "[.]"), function(x) x[1])


	return(list(mat=mat, meta_row=meta_row))
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

countModuleTissueStat = function(modules) {

	out = list()

	out$tissue_counts = sapply(1:max(modules$clust), function(cl) {
		table(modules$meta_genes$tissue[modules$clust == cl])
	})

	out$purity = apply(out$tissue_counts, 2, max) / apply(out$tissue_counts, 2, sum)

	out$n_tissues = apply(out$tissue_counts, 2, function(col) {
		sum(col > 0)
	})

	out$size = apply(out$tissue_counts, 2, sum)

	return(out)
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


loadMorbidObesity = function() {
	require(biomaRt)
	require(dplyr)

	# Get ensembl gene data
	message("Loading ensembl")
	ensembl = useEnsembl(biomart="ensembl",
		dataset="hsapiens_gene_ensembl")
	ensembl_symbol_map = getBM(
		attributes=c("ensembl_gene_id", "ensembl_transcript_id", "hgnc_symbol", "gene_biotype", "entrezgene"),
		mart=ensembl)

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

	return(list(emat=emat, meta_row=meta_row, sample_annot=sample_annot))
}