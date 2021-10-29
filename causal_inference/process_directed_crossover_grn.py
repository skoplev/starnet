#!/usr/bin/env python

import numpy as np
import pandas as pd

#define tissue
#for Blood use 'BLOOD'
tis=['LIV', 'SKLM', 'AOR', 'MAM', 'VAF', 'SF', 'BLOOD']

#import files
net10=pd.read_csv('full_pijs_gassist_GRN_cis-anchor_alt_10FDR.tsv', delimiter='\t')
grn=pd.read_csv('all_processed.tsv', delimiter='\t')

#align tissue nomenclature
#net10['tissue']=net10['tissue'].str.replace('Blood','BLOOD')

#merge directed interaction based on tissue
for x in tis:
	grn_t=grn.loc[(grn['head_tissue'] == x) & (grn['tail_tissue'] == x)]
	net10_t=net10.loc[net10['tissue'] == x]

	grn_net=pd.merge(grn_t, net10_t, left_on=['head_Gene_name', 'tail_Gene_name'], 
		right_on=['A-genes', 'B-genes'], how='inner')
	grn_abrv=grn_net[['tail_tissue', 'head_Gene_name', 'tail_Gene_name', 
			'head_Gene stable ID', 'tail_Gene stable ID', 'WEIGHT', 'Findr_score']]
	grn_abrv.columns=['tissue', 'head_Gene_name' , 'tail_Gene_name', 
			'head_Gene stable ID', 'tail_Gene stable ID', 'GRN_WEIGHT', 'Findr_score']
	grn_abrv.to_csv(x+'_key-drivers_directed_crossover_FDR10.tsv', sep='\t', index=None)
	print(x)

#merge reverse directed interaction based on tissue with tail and head flipped
for x in tis:
	grn_t=grn.loc[(grn['head_tissue'] == x) & (grn['tail_tissue'] == x)]
	net10_t=net10.loc[net10['tissue'] == x]

	grn_net=pd.merge(grn_t, net10_t, left_on=['tail_Gene_name', 'head_Gene_name'], 
		right_on=['A-genes', 'B-genes'], how='inner')
	grn_abrv=grn_net[['tail_tissue', 'head_Gene_name', 'tail_Gene_name', 
			'head_Gene stable ID', 'tail_Gene stable ID', 'WEIGHT', 'Findr_score']]
	grn_abrv.columns=['tissue', 'head_Gene_name' , 'tail_Gene_name', 
			'head_Gene stable ID', 'tail_Gene stable ID', 'GRN_WEIGHT', 'Findr_score']
	grn_abrv.to_csv(x+'_key-drivers_directed_crossover_reversed_FDR10.tsv', sep='\t', index=None)
	print(x)

