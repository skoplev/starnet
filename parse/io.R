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
