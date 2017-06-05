# Retrieve subset eQTL data from STARNET which contains specified SNP ID's.
# Must be run in project root.

# Files to process
files="
STARNET.eQTLs.MatrixEQTL.AOR.cis.tbl
STARNET.eQTLs.MatrixEQTL.Blood.cis.tbl
STARNET.eQTLs.MatrixEQTL.FC.cis.tbl
STARNET.eQTLs.MatrixEQTL.LIV.cis.tbl
STARNET.eQTLs.MatrixEQTL.MAM.cis.tbl
STARNET.eQTLs.MatrixEQTL.MP.cis.tbl
STARNET.eQTLs.MatrixEQTL.SF.cis.tbl
STARNET.eQTLs.MatrixEQTL.SKLM.cis.tbl
STARNET.eQTLs.MatrixEQTL.VAF.cis.tbl
"

# Original data folder
# in_data_folder="/sc/orga/projects/STARNET/oscar/Matrix_eQTL/adjusted.final/"
in_data_folder="/Users/sk/DataProjects/cross-tissue/eQTL/adjusted.final/"

# Target folder relative to project root.
target_data_folder="eQTL/adjusted/CAD/"
mkdir $target_data_folder

# Query file which contains list of SNP IDs to retrieve
query_file="eQTL/SNPs/CAD_deloukas_nikpay.txt"

for file in $files
do
	echo "Processing "$in_data_folder$file

	# Get first line and overwrite file in target folder
	head -n1 $in_data_folder$file > $target_data_folder$file

	# Grep data lines with workds containing
	grep -w -f $query_file $in_data_folder$file >> $target_data_folder$file
done