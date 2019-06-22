import re
import os
import sys

import numpy as np
import pandas as pd
from sklearn import preprocessing
from sklearn.decomposition import TruncatedSVD
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt  
import time


type1define = {'lineage_this_cell':3, 'lineage_other_cell':4, 'other_activator':8, 'inhibitor':9, 'inhibitor_lineage_this_cell':10, 'inhibitor_lineage_other_cell':11,
		'lineage_many_this':5, 'lineage_many_other':6,'lineage_all':7, 'terminal_specific_this':12, 'terminal_specific_other':15, 'terminal_2_this':13, 'terminal_2_other':16, 'terminal_all':14, 'tf':1, 'non_tf':2, 'lin':0}
typerev = {b:a for a, b in zip(type1define.keys(), type1define.values())}
cols = {'other_activator':'tab:blue', 'lineage_all':'tab:orange', 'lineage_this_cell':'tab:green', 'inhibitor':'tab:red', 'lineage_other_cell':'tab:purple', 'terminal_specific_this':'tab:brown', 'terminal_all':'tab:pink', 'terminal_specific_other':'tab:gray', 'terminal_2_other':'tab:olive', 'terminal_2_this':'tab:cyan'}

#pp='../cellevolver_shared/motifs_sep1_sepmot3_selby0_filtTrue/4cell_mce1fix_mutb/'
#tabname = '4cell_mce1fix_mutb_180528T130_geneByCellNormPosition_sepThreshold40.csv'

def makeTSNE(pp, tabname, perplexities=[30], reps = 3, outname = 'TSNE.pdf', motlist = None):
	df = pd.read_csv(pp+ '/' + tabname, index_col=0)
	if(motlist is None):
		mat = np.array(df[df.columns[14:]]) #remove motifs of size 2
	else:
		mat = np.array(df[[i for i in df.columns if i.split('_')[0] in motlist]])
	svd2 = TruncatedSVD(n_components=100, n_iter=20, random_state=42)
	mat2 = svd2.fit_transform(mat)
	for perplexity in perplexities:
		for rep in range(reps):
			tt = time.time()
			#x_embedded = mat2[:,:2]
			x_embedded = TSNE(n_components=2, perplexity=perplexity).fit_transform(mat2)
			#cols = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan']
			types = list(set(list(df.type2)))
			scale  = 10
			fig, ax = plt.subplots(figsize=(10,7))
			for i, t in enumerate(types):
				ind = [j for j,k in enumerate(list(df.type2)) if k==t]
				ax.scatter(x_embedded[ind,0], x_embedded[ind, 1], s=scale, c=[cols[typerev[t]] for h in ind], label=re.sub('_', ' ',typerev[t]), alpha=0.8, edgecolors='none')
			ax.legend(loc=4)#bbox_to_anchor=(1.5, 1))
			plt.savefig(outname + '_TSNEperx' + str(perplexity) + '_'+str(rep) + '.pdf')
			plt.close('all')
			print("perplexity " + str(perplexity) + " repetition " + str(rep) + " in " + str((time.time() - tt)/60) + " minutes. SVD total variance: " + str(sum(svd2.explained_variance_)
) + "\n")



def createIndexOfSignificantMotifs(tabfile, conditions, cutoff = 2):
	to_remove = ['4cell_mce0random_5', '4cell_mce0fixrandom_4', '4cell_mce0fixrandom_3', '4cell_mce0_mutbrandom_2'] #a mistake happened during motif search
	msum = pd.read_csv(tabfile, index_col=0)
	msum = msum.loc[msum['motif_size'] > 2]
	msum = msum.drop(to_remove, axis=1)
	for c in conditions:
		rand = np.array(msum[[i for i in msum.columns if c+'random_' in i]])
		randmean = np.mean(rand, axis = 1)
		randstd = np.std(rand, axis = 1)
		zscore = (np.array(msum[c]) - randmean + 1)/(randstd+1)
		normzscore = (zscore+1)/(np.sum(np.sqrt(np.power(zscore, 2)))+1)
		for value, name in zip([randmean, randstd, zscore, normzscore], [c+i for i in ['_randmean', '_randstd', '_zscore', '_normzscore']]):
			msum[name] = value
	msum.to_csv('all_Zscores.csv')
	minzscores = np.min(np.array(msum[[i for i in msum.columns if '_zscore' in i]]), axis=1)
	mots = [m for z,m in zip(minzscores, msum['motif_id']) if z >= cutoff]
	antimots = [m for z, m in zip(minzscores, msum['motif_id']) if z <= -1*cutoff]
	return {'motifs':mots, 'antimotifs':antimots, 'cutoff':cutoff}
		
def getfilename(path, filestring):
	fname = [i for i in os.listdir(path) if filestring in i]
	if(len(fname) == 1):
		fname = fname[0]
	else:
		raise Exception('more than 1 summary file in: ' + path + ' - ' + filestring)
	return fname

def main():
	conditions = ['4cell_mce0', '4cell_mce0fix_mutb', '4cell_mce1fix', '4cell_mce2fix', '4cell_mce2inh5', '4cell_mce0fix', '4cell_mce0_mutb', '4cell_mce1fix_mutb', '4cell_mce2fix_mutb', '4cell_mce2inh5_mutb', '4cell_mce0Xss']
	paths= ['../../cellevolver_shared/motifs_sep1_sepmot3_selby0_filtTrue_ignoreSelfReg_fixed1/' + i for i in conditions]
	lastcond = paths[-1]
	summary_name = '_motifSummary.csv'
	motbycellname = '_geneByCellNormPosition_sep.csv'
	last_summary_file = getfilename(lastcond, summary_name)
	if(len(sys.argv)>1):
		cutoff = int(sys.argv[1])
	else:
		cutoff = 2
	motiflist = createIndexOfSignificantMotifs(paths[-1] + '/' + last_summary_file, conditions, cutoff)
	for c, p in zip(conditions, paths):
		fname = getfilename(p, motbycellname)
		makeTSNE(p, fname,[30, 40], outname = c, motlist=motiflist['motifs'])



if __name__ == "__main__":
    main()

