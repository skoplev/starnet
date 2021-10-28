options(stringsAsFactors=F)

rm(list=ls())

library(MatrixEQTL)

args <- commandArgs(trailingOnly=TRUE)
# args = c("22","LIV","1")
chr = as.numeric(as.character(args[[1]]))
tissue = args[[2]]
num_SV = args[[3]]

gene_annotation = read.delim("/sc/orga/projects/STARNET/ariella/MatrixEQTLgenes.txt",header=T,quote="\"")
names(gene_annotation) = c("gene","chr","start_coord","end_coord")
gene_annotation$gene<-sapply(strsplit(as.character(gene_annotation$gene),'[.]'),function(x) x[1])

load('/sc/orga/projects/STARNET/koples01/data/cross-tissue/gene_exp_norm_batch/all.RData')
expression_data = expr_mats_batch[[tissue]]
rownames(expression_data)<-sapply(strsplit(as.character(rownames(expression_data)),'_'),function(x) x[2])
rownames(expression_data)<-sapply(strsplit(as.character(rownames(expression_data)),'[.]'),function(x) x[1])

SNP_file_name = paste("/sc/orga/projects/STARNET/ariella/genotype.dose.",chr,sep="")

## Setting up Matrix eQTL for running all linear models:
# Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS
errorCovariance = numeric();
pvOutputThreshold_cis  = 0.05
cisDist = 1e6;

## Load genotype data
snpsME = SlicedData$new();
snpsME$fileDelimiter = "\t";      # the TAB character
snpsME$fileOmitCharacters = "NA"; # denote missing values;
snpsME$fileSkipRows = 1;          # one row of column labels
snpsME$fileSkipColumns = 9;       # one column of row labels
snpsME$fileSliceSize = 2000;      # read file in slices of 2,000 rows
snpsME$LoadFile(SNP_file_name);

## Load gene expression data
geneME = SlicedData$new();
geneME$fileDelimiter = " ";      # the TAB character
geneME$fileOmitCharacters = "NA"; # denote missing values;
geneME$fileSkipRows = 1;          # one row of column labels
geneME$fileSkipColumns = 1;       # one column of row labels
geneME$fileSliceSize = 2000;      # read file in slices of 2,000 rows
# geneME$LoadFile(expression_file_name);
geneME$CreateFromMatrix(expression_data);


## making sure samples are the same in SNP and Gene Expression:
sample_snps= as.data.frame(cbind(colnames(snpsME))) 
sample_snps = cbind(1:ncol(snpsME),sample_snps)
names(sample_snps)=c("snp.idx","sample")
sample_genes= as.data.frame(cbind(unlist(lapply(colnames(geneME),function(x){strsplit(x,"_",fixed=T)[[1]][2]}))))
sample_genes = cbind(1:ncol(geneME),sample_genes)
names(sample_genes)=c("gene.idx","sample")

samples = merge(sample_snps,sample_genes,by=2)

#taking only the samples that are in both data types:
snpsME$ColumnSubsample( c(samples$snp.idx))
geneME$ColumnSubsample( c(samples$gene.idx))

cat("SNPs and Gene Expression match up?",identical(colnames(geneME),paste(tissue,colnames(snpsME),sep="_")),"\n")




## setting up the snp positions and gene positions
snpspos = read.table(paste("/sc/orga/projects/STARNET/ariella/chr",chr,"_snps",sep=""), header = TRUE, stringsAsFactors = FALSE);
snpspos.name = read.table(paste("/sc/orga/projects/STARNET/ariella/chr",chr,"_snps.name",sep=""),header=TRUE,stringsAsFactors=FALSE)

snp_positions = cbind(snpspos.name, snpspos)
snp_positions$chr= as.numeric(snp_positions$chr)
snp_positions$pos= as.numeric(snp_positions$pos)

file.out.cis <- paste0("/sc/orga/projects/STARNET/vamsi/Johan/eQTL/",tissue,"/cisQTLs_SV",num_SV,"_chr",chr,".cis.txt")
# Running Matrix EQTL: 
me <- Matrix_eQTL_main(
  snps = snpsME, #matrix eqtl formmatted sliced data of snp of interest
  gene = geneME, # expression
  cvrt = SlicedData$new(),
  output_file_name = NULL,
  output_file_name.cis = file.out.cis,
  pvOutputThreshold = 0,
  pvOutputThreshold.cis = 1,
  useModel = modelLINEAR, # modelANOVA or modelLINEAR or modelLINEAR_CROSS
  errorCovariance = numeric(),
  verbose = TRUE,
  pvalue.hist = FALSE,
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = TRUE,
  cisDist = 1000000,
  snpspos = snp_positions, # position of the snp of interest 
  genepos = gene_annotation)





