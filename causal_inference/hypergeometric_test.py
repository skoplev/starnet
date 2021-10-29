#!/usr/bin/env python

import pandas as pd
import numpy as np
import scipy.stats as ss
import sys
import qvalue  

#select from 'LIV', 'SKLM', 'AOR', 'MAM', 'VAF', 'SF', 'BLOOD'

t='BLOOD'
print(t)
#import tissue expression data and 10% findr results
if t == 'BLOOD':
	t_exp=pd.read_csv('BLOOD_STARNET_exp_genofilt.mat', delimiter='\t')
else:
	t_exp=pd.read_csv(''+t+'_STARNET_exp_genofilt.mat', delimiter='\t')

tf=pd.read_csv(t+'_pijs_gassist_GRN_cis-anchor_alt_10FDR.tsv', delimiter='\t')
#import GRN modules
mod=pd.read_csv('modules.csv', delimiter=',')
#import crossover results
cross=pd.read_csv(t+'_key-drivers_directed_crossover_FDR10.tsv', delimiter='\t')

#get genes to be tested from crossover analysis
genes=tf['A-genes'].drop_duplicates()

#perform hypergeometric test for every gene identfied from crossover
pv=[]
M_list=[]
n_list=[]
N_list=[]
k_list=[]
cl_list=[]
for g in genes:
	#get GRN module for gene
	kd_mod=mod.loc[(mod['tissue'] == t) & (mod['gene_symbol'] == g)]
	
	#get GRN module
	cl=kd_mod['clust'].values[0]
	cl_list.append(cl)

	#get full GRN module
	cl_mod=mod.loc[(mod['tissue'] == t) & (mod['clust'] == cl)]

	#get Findr 10% FDR network for gene 
	tf_g=tf.loc[tf['A-genes'] == g]

	#get GRN/ Findr crossover for gene
	tf_filt=tf_g.loc[tf_g['B-genes'] != g]
	mod_filt=cl_mod.loc[cl_mod['gene_symbol'] != g]
	tf_mod=pd.merge(tf_filt[['B-genes']], mod_filt[['gene_symbol']], left_on='B-genes',
					right_on='gene_symbol', how='inner')

	#total number of genes in tissue T
	M=len(t_exp)
	M_list.append(M)
	#number of genes from tissue T in G
	n=len(cl_mod)-1
	n_list.append(n)
	#number of targets of KD in Findr 10% network for tissue T
	N=len(tf_g.drop_duplicates('B-genes'))
	N_list.append(N)
	#number of genes in the intersection between (genes from tissue T in G) and (targets of KD in Findr 10% network for tissue T)
	k=len(tf_mod.drop_duplicates())
	k_list.append(k)

	#conduct hypergeometric test
	hpd=ss.hypergeom(M, n, N)
	p=hpd.pmf(k)
	pv.append(p)

pvn=np.array(pv)
qv=qvalue.estimate(pvn)
out=pd.DataFrame({'key_driver': genes, 'clust': cl_list, 'tissue_genes': M_list, 'GRN_size-1': n_list,
				'Findr_targets_10FDR': N_list, 'intersection': k_list, 'p-value': pv, 'q-value': qv})
out.to_csv(t+'_key-driver_FDR10_hypergeometric.tsv', sep='\t', index=None)
