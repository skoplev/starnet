#!/bin/bash
# Submit job to server

module load rstudio

mkdir logs

for chr in `seq 1 1 22`; do
for tissue in LIV MAM AOR SF VAF SKLM Blood;do
for num_SV in `seq 1 1 20`;do 
name=${tissue}_SV${num_SV}_chr${chr}; 
echo $name;

bsub -J ${name} \
	-P acc_STARNET \
	-q expressalloc \
	-W 2:00 \
	-n 1 \
	-R rusage[mem=100000] \
	-eo /sc/orga/projects/STARNET/vamsi/Johan/eQTL/logfiles/${name}.stderr \
	-oo /sc/orga/projects/STARNET/vamsi/Johan/eQTL/logfiles/${name}.stdout \
	Rscript /sc/orga/projects/STARNET/vamsi/Johan/eQTL_final.R ${chr} ${tissue} ${num_SV}";
done; 
done; 
done
	
	

