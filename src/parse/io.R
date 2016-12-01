# Various input/output functionality

# Load output tables from CIBERSORT algorithm in directory
# Returns list of data.table objects.
loadCibersortFreq = function(path) {
	freq = list()
	for (file_name in list.files(path, pattern="*.tsv")) {
		message(file_name)
		file_path = file.path(path, file_name)

		# Read data without header
		d = fread(
			file_path,
			header=FALSE,
			skip=1
		)

		# Read and parse first line as header
		headers = strsplit(
			readLines(file_path, n=1),
			"\t"
		)[[1]]

		colnames(d) = c("tissue_patient_id", headers)

		# save in list
		freq[[file_name]] = d
	}
	return(freq)
}

# Matched set of cibersort frequency tables to STARNET patient IDs
# Trims non-frequency fields by column position.
# input: list of data.tables containing output from CIBERSORT.
matchTrimCibersortFreq = function(cibersort_freq, patient_ids) {
	cibersort_freq_matched = lapply(cibersort_freq, function(freq) {
		# name rows
		rownames(freq) = freq$tissue_patient_id

		# Get patient ids of samples in frequency set
		ids = sapply(
			strsplit(freq$tissue_patient_id, "_"),
			function(x) x[2]
		)
		# reorder rows, and exclude non frequency entries
		freq = freq[match(patient_ids, ids), 2:(ncol(freq) - 3), with=FALSE]

		return(freq)
	})

	return(cibersort_freq_matched)
}