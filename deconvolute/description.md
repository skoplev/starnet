# Generation of CIBERSORT basis matrix

1) Gene expression data collection
PHANTOM5 CAGE gene expression data from mouse and human downloaded from
http://fantom.gsc.riken.jp/5/datafiles/latest/extra/CAGE_peaks/

2) phantomFilterCAGE.r
Filters CAGE gene expression data: primary cells.
Collapses to HUGO gene symbols.
Outputs primary_cells.tsv expression matrix

3) phantomFeatureSelect.py
Performs recursive feature elimination, one-to-one classifier
Outputs a series of reduced basis matrices. Aggregates biological repeats by averaging.


# Preprocessing of STARNET gene expression data

1) starnetNormalize.r
Normalization of STARNET gene expression data across tissue.
