#!/usr/bin/env python

from sklearn.svm import SVC  # C-support vector classiciation, libsvm
from sklearn.feature_selection import RFE

from sklearn.model_selection import StratifiedKFold
from sklearn.feature_selection import RFECV

import matplotlib.pyplot as plt

import pandas as pd
import numpy
import os

basis_file = "/Users/sk/DataProjects/cross-tissue/FANTOM5/cell_type_basis/primary_cells.tsv"
out_folder = "/Users/sk/DataProjects/cross-tissue/FANTOM5/cell_type_basis"

# load data into panda data frame
# appends .# to columns with identical names.
expr = pd.read_csv(basis_file,
	sep="\t",
	index_col=0  # gene symbols as row names
)

# # Column names of expression data, panda forces unique names
# # Strip colnames of .# added
target = [elem.split(".")[0] for elem in expr.columns]

# Rename column names, for output of expression matrix
expr.columns = target
# len(numpy.unique(target))


# Recursive feature elimination, support vector machine multiclass on the cell types.
svc = SVC(kernel="linear", C=1)  # supervised estimator

rfe = RFE(estimator=svc,
	n_features_to_select=1000,  # stop conditions
	step=100,  # step size, features removed at each iteration
	verbose=1
)

# Run model
rfe.fit(expr.transpose(), target)

# Prints ranking, best features are ranked 1
# rfe.ranking_.tolist()
# Boolean array of selected features
# rfe.support_
# Get matrix of selected
# rfe.ranking_


# Write sub matrices
# 1000 rows
expr[rfe.support_].to_csv(
	os.path.join(out_folder, "primary_cells_svm1000.tsv"),
	sep="\t")

# Aggregate cell types
expr_aggregate = expr.groupby(by=expr.columns, axis=1).mean()

expr_aggregate[rfe.support_].to_csv(
	os.path.join(out_folder, "primary_cells_svm1000_aggregate.tsv"),
	sep="\t")


# gene_rank = 11  # second basis matrix, ~2000 rows
# gene_rank = 21  # second basis matrix
# gene_rank = 31  # second basis matrix
# gene_rank = 41  # second basis matrix


for gene_rank in [11, 21, 31, 41]:
	num_genes = numpy.sum(rfe.ranking_ <= gene_rank)  # number of genes for gene_rank cutoff
	# write to
	expr[rfe.ranking_ <= gene_rank].to_csv(
		os.path.join(
			out_folder,
			"primary_cells_svm" + str(num_genes) + ".tsv"),
		sep="\t")
	# Aggregate multiple columns by their mean
	expr_aggregate[rfe.ranking_ <= gene_rank].to_csv(
		os.path.join(
			out_folder,
			"primary_cells_svm_" + str(num_genes) + "_aggregate.tsv"),
		sep="\t")




# Cross validation to identify the optimal number of features to select
# -----------------------------------------
# rfecv = RFECV(estimator=svc,
# 	step=1000,
# 	cv=StratifiedKFold(2),  # k-fold
# 	scoring="accuracy",
# 	verbose=1
# )

# rfecv.fit(expr.transpose(), target)

# rfecv.n_features_

# plt.figure()
# plt.xlabel("Features selected")
# plt.ylabel("CV accuracy")
# plt.plot(
# 	range(1, len(rfecv.grid_scores_) + 1),
# 	rfecv.grid_scores_)
# plt.show()
