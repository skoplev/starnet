#!/usr/bin/env python

import pandas as pd
import numpy as np
import math
import findr

""""Iterates aross list of genes with best cis-eQTL and calculates all 5 tests from Findr 
using best eQTLs as causal anchors. Manual calculation of pij_gassist, pij_gassist_trad 
and P2*P5 test combination and returns 3 DataFrames"""

#state specific tissue e.g. LIV, SKLM, MAM, AOR, VAF, SF, BLOOD and file containing genes and eQTLs of interest
t='BLOOD'
genes=pd.read_csv(t+'_KD_cis-eQTLs_adjpv0.05.tsv', delimiter='\t')

#import expression data
t_mrna=pd.read_csv(t+'_STARNET_exp_genofilt.mat', delimiter='\t')

if t == 'BLOOD':
	nm=pd.read_csv('BLOOD_gene_annotation_si', delimiter='\t')
else:
	nm=pd.read_csv(t+'_gene_annotation_si', delimiter='\t')

t_mrna.rename(columns={'id':'Gene stable ID'}, inplace=True)

trad_out=[]
novel_out=[]
alt_out=[]

#filter to get eQTL genotypes
df_list=[]
genes.rename(columns={'SNP':'marker_id'}, inplace=True)
rs=genes[['marker_id']]
if t == 'BLOOD':
	for cnk in pd.read_csv('BLOOD.STARNET.si.tsv',
			sep='\t', iterator=True, chunksize=200000):
			mer=pd.merge(cnk, rs, on='marker_id', how='inner')
			df_list.append(mer)
else:
	for cnk in pd.read_csv(t+'.STARNET.si.tsv',
			sep='\t', iterator=True, chunksize=200000):
			mer=pd.merge(cnk, rs, on='marker_id', how='inner')
			df_list.append(mer)
trans_genotypes=pd.concat(df_list)

for index, row in genes.iterrows():
	#get eQTLs
	eQTLs=pd.DataFrame({'Gene_name': [row[1]], 'Gene stable ID': [row[2]], 'SNP': [row[3]], 
		'beta': [row[6]], 't-stat': [row[7]], 'p-value': [row[8]], 'adj.p-value': [row[9]], 
		'tissue': [row[0]]})

	print(eQTLs['Gene_name'].values[0])
	#filter to get row genotype
	eQTLs.rename(columns={'SNP':'marker_id'}, inplace=True)
	rs_row=eQTLs[['marker_id']]
	trans_geno=pd.merge(rs_row, trans_genotypes, on='marker_id', how='inner')

	#obtain expression data from gene list
	trans_nm=pd.DataFrame(eQTLs['Gene stable ID'])
	trans_mrna=pd.merge(trans_nm, t_mrna, on='Gene stable ID', how='left')

	#remove duplicate genes so there is expression data for specificed eQTLs
	mrna_index=t_mrna.loc[t_mrna['Gene stable ID'].isin(trans_mrna['Gene stable ID'].tolist())]
	mrnai=mrna_index.reset_index()
	mrna_index_sorted=pd.merge(trans_nm, mrnai, on='Gene stable ID', how='left')
	mrna_index_v=mrna_index_sorted['index'].values

	#pass conditions and error output
	if math.isnan(mrna_index_v[0]) == True:
		with open(t+'_error.log', 'a') as e_out:
			e_out.write(eQTLs['Gene_name'].values[0]+' failed: skipping...'+'\n')
			e_out.close()
		print(eQTLs['Gene_name'].values[0]+' failed: skipping...')

		continue

	else:

		pass

	t_mrna_sorted=pd.concat([t_mrna.iloc[mrna_index_v,:],
		t_mrna.drop(mrna_index_v, axis=0)], axis=0)


	#remove columns and convert to format for findr
	dg=trans_geno.drop('marker_id', axis=1)
	d=trans_mrna.drop('Gene stable ID', axis=1)
	dt=t_mrna_sorted.drop('Gene stable ID', axis=1)

	dgi=dg.values
	di=d.values
	dti=dt.values

	dg_input=dgi.astype(np.uint8)
	d_input=di.astype(np.float32)
	dt_input=dti.astype(np.float32)

	#pass conditions and error output
	if np.shape(dg_input) != np.shape(d_input):
		with open(t+'_error.log', 'a') as e_out:
			e_out.write(eQTLs['Gene_name'].values[0]+' failed: skipping...'+'\n')
			e_out.close()
		print(eQTLs['Gene_name'].values[0]+' failed: skipping...')

		continue

	else:

		pass

	#run novel causal inference test in Findr
	l=findr.lib(path='libfindr.so',loglv=6,rs=0,nth=0)
	ans=l.pijs_gassist(dg_input,d_input,dt_input,na=None,nodiag=True,memlimit=1000000000)
	ans

	#create x and y labels
	xlab=t_mrna_sorted['Gene stable ID']
	ylab=pd.merge(trans_mrna, genes, on='Gene stable ID', how='left').drop_duplicates(subset='Gene_name')

	#generate and orientate dataframes
	p1=pd.DataFrame(ans['p1'])
	p1t=p1.T.reset_index()
	p2=pd.DataFrame(ans['p2'], columns=xlab)
	p2t=p2.T.reset_index()
	p3=pd.DataFrame(ans['p3'], columns=xlab)
	p3t=p3.T.reset_index()
	p4=pd.DataFrame(ans['p4'], columns=xlab)
	p4t=p4.T.reset_index()
	p5=pd.DataFrame(ans['p5'], columns=xlab)
	p5t=p5.T.reset_index()

	#replace Ensembl ID with genes names
	p1t.columns=['Gene stable ID']+list(ylab['Gene_name'])
	p2t.columns=['Gene stable ID']+list(ylab['Gene_name'])
	p3t.columns=['Gene stable ID']+list(ylab['Gene_name'])
	p4t.columns=['Gene stable ID']+list(ylab['Gene_name'])
	p5t.columns=['Gene stable ID']+list(ylab['Gene_name'])

	p2_comp=pd.merge(p2t, nm, on='Gene stable ID', how='inner')
	p2_comp=p2_comp.drop(['Unnamed: 0', 'Transcript stable ID', 'Chromosome/scaffold name',
	 'Gene start (bp)', 'Gene end (bp)', 'Gene stable ID', 'Gene type'], axis=1)
	p3_comp=pd.merge(p3t, nm, on='Gene stable ID', how='inner')
	p3_comp=p3_comp.drop(['Unnamed: 0', 'Transcript stable ID', 'Chromosome/scaffold name',
	 'Gene start (bp)', 'Gene end (bp)', 'Gene stable ID', 'Gene type'], axis=1)
	p4_comp=pd.merge(p4t, nm, on='Gene stable ID', how='inner')
	p4_comp=p4_comp.drop(['Unnamed: 0', 'Transcript stable ID', 'Chromosome/scaffold name',
	 'Gene start (bp)', 'Gene end (bp)', 'Gene stable ID', 'Gene type'], axis=1)
	p5_comp=pd.merge(p5t, nm, on='Gene stable ID', how='inner')
	p5_comp=p5_comp.drop(['Unnamed: 0', 'Transcript stable ID', 'Chromosome/scaffold name',
	'Gene start (bp)', 'Gene end (bp)', 'Gene stable ID', 'Gene type'], axis=1)

	#reorder columns so Gene Name comes first
	cols=p2_comp.columns.tolist()
	cols=cols[-1:] + cols[:-1]
	p2_comp=p2_comp[cols]
	p3_comp=p3_comp[cols]
	p4_comp=p4_comp[cols]
	p5_comp=p5_comp[cols]

	#remove gene columns for test calculation
	p2_m=p2_comp.drop('Gene_name', axis=1)
	p3_m=p3_comp.drop('Gene_name', axis=1)
	p4_m=p4_comp.drop('Gene_name', axis=1)
	p5_m=p5_comp.drop('Gene_name', axis=1)

	#Get results of P2*P5 test
	p2p5=p2_m*p5_m
	p2p5['Gene_name']=p2_comp['Gene_name']
	p2p5=p2p5[cols]
	p2p5=p2p5.set_index('Gene_name')
	p2p5=p2p5.stack().reset_index()
	p2p5.columns=['B-genes', 'A-genes', 'Findr_score']
	p2p5_df=p2p5[['A-genes', 'B-genes', 'Findr_score']]
	alt_out.append(p2p5_df)

	#Get results for pij_gassist_trad
	trad=p2_m*p3_m
	trad['Gene_name']=p2_comp['Gene_name']
	trad=trad[cols]
	trad=trad.set_index('Gene_name')
	trad=trad.stack().reset_index()
	trad.columns=['B-genes', 'A-genes', 'Findr_score']
	trad_df=trad[['A-genes', 'B-genes', 'Findr_score']]
	trad_out.append(trad_df)

	#Get results for pij_gassist
	gas=(0.5*((p2_m*p5_m)+p4_m))
	gas['Gene_name']=p2_comp['Gene_name']
	gas=gas[cols]
	gas=gas.set_index('Gene_name')
	gas=gas.stack().reset_index()
	gas.columns=['B-genes', 'A-genes', 'Findr_score']
	gas_df=gas[['A-genes', 'B-genes', 'Findr_score']]
	novel_out.append(gas_df)

#export data as dataframes
alt_exp=pd.concat(alt_out)
trad_exp=pd.concat(trad_out)
novel_exp=pd.concat(novel_out)
alt_exp.to_csv(t+'_pijs_gassist_GRN_cis-anchor_alt', sep='\t', index=None)
trad_exp.to_csv(t+'_pijs_gassist_GRN_cis-anchor_trad', sep='\t', index=None)
novel_exp.to_csv(t+'_pijs_gassist_GRN_cis-anchor_novel', sep='\t', index=None)
