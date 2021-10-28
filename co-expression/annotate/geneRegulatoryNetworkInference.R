# Infers gene regulatory networks within each co-expression module using GENIE3
#
# --------------------------------------------------

rm(list=ls())

library(data.table)
library(igraph)
library(Mergeomics)
library(GENIE3)
library(reshape2)
library(RColorBrewer)


library(compiler)
enableJIT(3)


setwd("~/GoogleDrive/projects/STARNET/cross-tissue")
data_dir = "/Users/sk/DataProjects/cross-tissue"  # root of data directory

source("src/models/bayesNet.R")
source("src/parse.R")


learnRegNetw = function(
	clust,
	meta_genes,
	emat,
	source_symbols,
	# nCores=6,
	max_size=3000)
{
	if (nrow(meta_genes) != nrow(emat)) {
		stop("Expression matrix and transcript meta data mismatch.")
	}

	if (length(clust) != nrow(meta_genes)) {
		stop("Clust vector length differ from meta_genes.")
	}

	netws = lapply(1:max(clust, na.rm=TRUE), function(k) {
		message("Module ", k)
		row_idx = which(clust == k)

		# Check module size
		if (length(row_idx) > max_size) {
			return(NA)
		}

		if (length(row_idx) < 10) {
			return(NA)
		}

		# Print first 10 transcript IDs
		message("Module: ", paste(meta_genes$transcript_id[row_idx][1:2], collapse=", "), "...(n=", length(row_idx), ")")

		meta_genes[row_idx, ]

		# Get matrix
		submat = t(emat[row_idx,])
		colnames(submat) = meta_genes$tissue_transcript_id[row_idx]

		# Exclude samples that have no observations
		col_idx = apply(submat, 1, function(col) {
			return(!all(is.na(col)))
		})

		submat = submat[col_idx, ]

		# Effective sample size
		message("Sample size: ", sum(col_idx))

		# colnames(submat) = mod_ids

		# Standardize and impute missing to mean
		submat = scale(submat)
		submat[is.na(submat)] = 0.0

		if (sum(meta_genes$regulator[row_idx]) < 2) {
			message("Module ", k, " has less than 2 regulators.")
			return(NA)
		}

		netw = GENIE3(
			t(submat),
			regulators=which(meta_genes$regulator[row_idx]),
			# nCores=6,
			nCores=1,
			# nTrees=1000,
			nTrees=100,
			# nTrees=10,
			verbose=TRUE
		)

		return(netw)
	})
}


# Load cross-tissue modules
# ------------------------------
between = new.env()
load(file.path(data_dir, "modules/between_within-cross-tissue.RData"),
	between,
	verbose=TRUE)

between = parseModuleData(between)


# Load normalized gene expression data, used for evaluating Bayesian networks within cross-tissue modules
# --------------------------------------------------------------
emat_file = "STARNET/gene_exp_norm_reshape/expr_recast.RData"

load(file.path(data_dir, emat_file), verbose=TRUE)

emat = expr_recast[, 3:ncol(expr_recast)]
meta_genes = expr_recast[, 1:2]
meta_genes = as.data.frame(meta_genes)

meta_genes = parseTranscriptId(meta_genes)
meta_genes$tissue_ensembl_id = paste0(meta_genes$tissue, "_",
	sapply(strsplit(meta_genes$ensembl, "[.]"), function(x) x[1])
)
	
rm(expr_recast)

# Test if the gene metadata is the same as for the cross-tissue modules
if (!all(meta_genes$transcript_id == between$meta_genes$transcript_id)) {
	stop("Transcript mismatch")
}


# Load TF definition from Lambert et al
# -------------------------------------------------------
tf_symbols = as.character(read.table("transcription-factors/lambert/TF_names_v_1.01.txt")$V1)


# Load STARNET cis-eQTL, Old data from Science paper
# Returns vector of tissue_ensembl IDs
# ------------------------------------------
getEqtlOrig = function() {
	message("Loading eQTL data...")
	cis_eqtl_dir = "/Users/sk/DataProjects/STARNET/eQTL"
	cis_eqtl_files = list.files(cis_eqtl_dir, pattern="*.cis.tbl")

	# Load cis-eQTL
	cis_eqtl = lapply(
		cis_eqtl_files,
		function(file_name) {
			d = fread(file.path(cis_eqtl_dir, file_name))
			d = d[d$adj.p < 0.05, ]  # FDR < 5%
			d = d[order(d$p), ]
			return(d)
	})
	names(cis_eqtl) = sapply(strsplit(cis_eqtl_files, "[.]"), function(x) x[4])  # name by tissue code

	# Add tissue_ensembl_id column
	cis_eqtl = lapply(1:length(cis_eqtl), function(i) {
		tissue = toupper(names(cis_eqtl)[i])
		message(tissue)

		d = cis_eqtl[[i]]

		# d$tissue_ensembl_id = paste0(tissue, "_", d$gene)

		# Strip version from ensembl IDs
		d$tissue_ensembl_id = paste0(tissue, "_",
			sapply(strsplit(d$gene, "[.]"), function(x) x[1])
		)


		return(d)
	})
	names(cis_eqtl) = sapply(strsplit(cis_eqtl_files, "[.]"), function(x) x[4])  # name by tissue code

	# Get list of all significant cis-eQTL genes
	cis_eqtl_all = unlist(sapply(cis_eqtl, function(x) x$tissue_ensembl_id))
	cis_eqtl_all = unique(cis_eqtl_all)

	return(cis_eqtl_all)
}

# Load STARNET cis-eQTL data, estimated by Vamsi
# -----------------------------------------------------
getEqtlNew = function() {

	cis_eqtl_dir = "~/DataProjects/STARNET/vamsi_eQTL/adjusted.final"
	cis_eqtl_files = list.files(cis_eqtl_dir, pattern="*.tbl")

	tissues = sapply(strsplit(cis_eqtl_files, "_"), function(x) x[1])

	# Rename tissue codes
	tissues[tissues == "SKM"] = "SKLM"
	tissues[tissues == "SUF"] = "SF"
	tissues[tissues == "BLO"] = "BLOOD"


	cis_eqtl = lapply(
		cis_eqtl_files,
		function(file_name) {
			d = fread(file.path(cis_eqtl_dir, file_name))
			# d = d[d$padj_fdr < 0.05, ]  # FDR < 5%
			# d = d[d$padj_fdr < 0.01, ]  # FDR < 1%
			# d = d[d$padj_fdr < 0.001, ]  # FDR < .1%
			d = d[d$padj_fdr < 0.0001, ]  # FDR < .01%
			d = d[order(d[["p-value"]]), ]
			return(d)
	})
	names(cis_eqtl) = tissues

	# Add tissue information to table
	for (i in 1:length(cis_eqtl)) {
		cis_eqtl[[i]]$tissue = tissues[i]
	}

	# Exclude macrophage eQTL
	cis_eqtl = cis_eqtl[-which(names(cis_eqtl) == "MAC")]

	# Combine tables
	cis_eqtl = rbindlist(cis_eqtl)

	cis_eqtl$tissue_ensembl_id = paste(cis_eqtl$tissue, cis_eqtl$gene, sep="_")  # tissue ensembl IDs for matching with module assignments

	return(unique(cis_eqtl$tissue_ensembl_id))
}



# cis_eqtl_all = getEqtlOrig()  # eQTL published in Science paper
cis_eqtl_all = getEqtlNew()  # Vamsi's eQTL

length(cis_eqtl_all)


length(unique(sapply(strsplit(cis_eqtl_all, "_"), function(x) x[2])))

# Specify regulator nodes: transcription factor symbols or eQTL detected at FDR < 0.1
# ------------------------------------------------------
meta_genes$regulator = FALSE
meta_genes$regulator[meta_genes$gene_symbol %in% tf_symbols] = TRUE
sum(meta_genes$regulator)

meta_genes$regulator[meta_genes$tissue_ensembl_id %in% cis_eqtl_all] = TRUE
sum(meta_genes$regulator)

mean(meta_genes$regulator)


# Infer gene regulatory networks, run with nTrees=1000
netws = learnRegNetw(
	clust=between$clust,
	meta_genes,
	emat,
	max_size=3000)


# Infer gene regulatory networks, run with nTrees=100
netws_global = learnRegNetw(
	# clust=rep(1, length(between$clust)),  # single cluster containing all genes
	# clust=between$clust %% 2,  # modulo, 2 clusters total, merged at 'random'
	clust=as.integer(meta_genes$tissue),  # per tissue 'clusters'
	meta_genes,
	emat,
	max_size=1e6)


# Old eQTL
# saveRDS(netws, file="co-expression/annotate/grn/netws.rds")
# netws = readRDS(file="co-expression/annotate/grn/netws.rds")

# New eQTL
# saveRDS(netws, file="co-expression/annotate/grn/netws_vamsi_eqtl.rds")
netws = readRDS(file="co-expression/annotate/grn/netws_vamsi_eqtl.rds")

# saveRDS(netws_global, file="co-expression/annotate/grn/netws_global_vamsi_eqtl.rds")
# netws_global = readRDS(file="co-expression/annotate/grn/netws_global_vamsi_eqtl.rds")



# Convert networks to edge lists
# -------------------

netwsToEdgeList = function(netws, weight_thresh=0.05) {
	# edsgelist
	edge_list = lapply(netws, function(adj) {

		if (is.na(adj)) {
			return(NA)
		} else {
			el = melt(adj)
			return(el)
		}
	})

	# Remove entries for large modules without inferred networks
	edge_list = edge_list[!is.na(edge_list)]


	# Combine edge lists from co-expression modules
	# edge_list = Reduce(rbind, edge_list)
	edge_list = rbindlist(edge_list)

	colnames(edge_list) = c("TAIL", "HEAD", "WEIGHT")
	edge_list = as.data.frame(edge_list)
	plot(density(edge_list$WEIGHT))

	summary(edge_list$WEIGHT)
	sum(edge_list$WEIGHT > weight_thresh)

	# Remove low weight edges
	edge_list = edge_list[edge_list$WEIGHT > weight_thresh, ]
}

edge_list = netwsToEdgeList(netws)

edge_list_global = netwsToEdgeList(netws_global)

edge_list_global_thresh01 = netwsToEdgeList(netws_global, weight_thresh=0.01)


edge_list_merge = merge(edge_list, edge_list_global_thresh01,
	all=TRUE,
	suffixes=c(".module", ".global"),
	by=c("TAIL", "HEAD"))


edge_list_merge[is.na(edge_list_merge)] = 0



# parl(mfrow=c(1, 2))

found_global = edge_list_merge$WEIGHT.global > 0.01

values = c(
	"Global by tissue"=sum(edge_list_merge$WEIGHT.module > 0.05 & found_global),
	"Module or CT exclusive"=sum(edge_list_merge$WEIGHT.module > 0.05 & !found_global)
)

pdf("co-expression/annotate/globalNetwComp/plots/number_edges_global.pdf", width=2)
barplot(
	t(t(values)),
	legend.text=TRUE,
	ylab="GENIE3 edges identified by co-expression module",
	# col=brewer.pal(9, "Pastel1")
)
dev.off()

edge_list_merge = edge_list_merge[order(edge_list_merge$WEIGHT.module, decreasing=TRUE), ]
# edge_list_merge = edge_list_merge[order(edge_list_merge$WEIGHT.global, decreasing=TRUE), ]




pdf("co-expression/annotate/globalNetwComp/plots/top_edges_global.pdf")
par(mfrow=c(2, 1))
n = 2000
plot(edge_list_merge$WEIGHT.module[1:n],
	xlab="",
	ylab="Module GENIE3 weight",
	ylim=c(0,1),
	cex=0.5, col=brewer.pal(9, "Set1")[1])
plot(edge_list_merge$WEIGHT.global[1:n],
	ylab="Global GENIE3 weight",
	xlab="Ranked gene-gene interactions",
	cex=0.5, col=brewer.pal(9, "Set1")[2])
dev.off()



# Key driver analysis
# -------------------------------------------

# Make module table for KDA.
kda_mod_tab = data.frame(
	MODULE=between$clust,
	NODE=meta_genes$tissue_transcript_id
)



# Write edge list to file, used for input to Weighted Key Driver analysis.
# netw_dir = "co-expression/annotate/grn"  # old eQTL
netw_dir = "co-expression/annotate/grn_vamsi_eqtl"  # Vamsi's new eQTL

# netw_dir = "co-expression/annotate/grn_vamsi_eqtl"  # Vamsi's new eQTL
write.table(edge_list,
	file.path(netw_dir, "all.tsv"),
	sep="\t",
	row.names=FALSE,
	quote=FALSE
)

write.table(kda_mod_tab,
	file.path(netw_dir, "nodes.tsv"),
	sep="\t",
	row.names=FALSE,
	quote=FALSE
)


# Weighted key driver analysis (wKDA) of gene regulatory networks within each cross-tissue module
# kda_label = paste0(gene, "_key_driver")
# kda_label = "modules"
kda_label = "modules.directed"
job.kda = list()
job.kda$netfile = file.path(netw_dir, "all.tsv")
job.kda$modfile = file.path(netw_dir, "nodes.tsv")
job.kda$label = kda_label
job.kda$folder = netw_dir
job.kda$nperm = 20


# These 3 are additional arguments in the 'directed' analysis. THey have not been used in other analyses.
job.kda$direction = -1  # downstream directed
job.kda$edgefactor = 1.0  # consider edge weights
job.kda$depth = 3  # search depth


job.kda = kda.configure(job.kda)
job.kda = kda.start(job.kda)
job.kda = kda.prepare(job.kda)
job.kda = kda.analyze(job.kda)
job.kda = kda.finish(job.kda)

# load results table
kda_results = read.table(
	file.path(netw_dir, "kda", paste0(kda_label, ".results.txt")),
	header=TRUE
)


kda_results[kda_results$MODULE == 98, ]
kda_results[kda_results$MODULE == 78, ]


