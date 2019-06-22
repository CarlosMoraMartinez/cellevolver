import os
import sys
import pandas as pd
import numpy as np
import copy as cp
import matplotlib.pyplot as plt
import seaborn as sns; sns.set()
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import gridspec
from collections import namedtuple

#produce combined heatmaps 

Filetypes = namedtuple("Filetypes", "expr ann tfs regtable timexpr mutanttfs mutantsites interactions")
class simulationHolder(object):
	ACTIVATOR_TYPE = '1'
	INHIBITOR_TYPE = '-1'
	filenames = Filetypes("finalExpression",  "annotation", "TFset", "RegulationTable", "timeExpression", "mutantTFs", "mutantSites", "TFinteractions")
	def __init__(self, name, path):
		self.name = name
		#print(name, '\n')
		self.expr = self.readData(name, self.filenames.expr, path)
		self.ann = self.readData(name, self.filenames.ann, path)
		self.tfs = self.readData(name, self.filenames.tfs, path)
		self.timexpr = self.readData(name, self.filenames.timexpr, path)
		self.mutanttfs = self.readData(name, self.filenames.mutanttfs ,path)
		self.mutantsites = self.readData(name, self.filenames.mutantsites, path)
		self.interactions = self.readData(name, self.filenames.interactions, path)
		self.genome_size = self.expr.shape[0]
		self.ncells = len([i for i in self.expr.columns if 'opt' in i])
		#self.regtable = self.readData(name, self.filenames.regtable, path)
		self.regtable = self.makeRegTable()
		self.tf_inds = self.getRegulators()
		self.activators = self.getActivators()
		self.inhibitors = self.getInhibitors()	
		self.type = self.expr.type
		self.type2 = self.__setType2()
		self.type3 = self.__setType3()
		self.motifList = None
	def __setType2(self):
		t = cp.copy(self.expr.type)
		iniexpr = np.array(self.getInitialExpr())
		for i in range(len(t)):
			if(t[i] == self.ACTIVATOR_TYPE or t[i] == self.INHIBITOR_TYPE):
				t[i] = t[i] + 'lin' + ''.join([str(c) for c in np.where(iniexpr[i, :] > 0)])
		return t
	def __setType3(self):
		import re
		t = cp.copy(self.type2)
		for i in range(len(t)):
			m = re.search('(\\[.+?\\])', t[i])
			if(m is not None):
				aux = m.group(1)
				aux = '[' + str(len(re.sub("\D", "", aux))) + ']'
			else:
				aux = '[0]'
			t[i] = re.sub("\\[.*\\]", aux, t[i])
		return t			
	def readData(self, pckname, filetype, path):
		f = os.listdir(path)
		datafile = [i for i in f if i.endswith('.csv') and pckname in i and filetype in i]
		assert len(datafile) == 1, "More than one possible file: " + pckname + ' - ' + filetype
		datafile = datafile[0]
		data = pd.read_csv('/'.join([path, datafile]), sep = "\t", header=0, index_col=0)
		return data
	def makeRegTable(self):
		wt = self.expr
		d = self.ann
		reg = pd.DataFrame()
		reg['gene'] = [i for i in range(self.genome_size)]
		reg['type'] = wt.type
		tf_names = self.getRegulators()
		#generate interaction table
		for tf in tf_names:
			aux = d[:][d.tf_ind == tf]
			if(aux.shape[0] > 0):
				reg['Complete' + str(tf)] = [sum(aux.percent[aux.gene == i]) for i in reg.gene]
				if(reg['type'][tf] == '-1'):
					reg['Complete' + str(tf)] = -1*reg['Complete' + str(tf)]
				for c in range(len([i for i in wt.columns if 'exp' in i])):
					newname = str(tf) + 'Cell ' + str(c)
					reg[newname] = reg['Complete' + str(tf)]*wt['exp'+str(c)][tf]
			else:
				reg['Complete' + str(tf)] = [0 for i in reg.gene]
				for c in range(len([i for i in wt.columns if 'exp' in i])):
					newname = str(tf) + 'Cell ' + str(c)
					reg[newname] = reg['Complete' + str(tf)]
		return reg
	def getTargets(self, tf, full = False):
		r = self.ann[self.ann.tf_ind == tf]
		if (not full):
			return sorted(list(set(r['gene'])))
		else:
			return sorted(r)
	def getTargetTypes(self, tf):
		targets = self.getTargets(tf, False)
		return sorted(self.expr.type[targets])
	def getRegulators(self, target = None, cell = None, filt = True, thresholdExpr = 0, types = ['1', '-1']):
		if(target is None):
			reg = [i for i in range(self.genome_size) if(self.expr.type[i] in types)]
		else:
			reg = [i for i in range(self.genome_size) if(self.expr.type[i] in types and i in self.ann.tf_ind[self.ann.gene == target].tolist())]
		if(cell is not None and filt):
			cellexpr = self.expr['exp'+str(cell)]
			reg = [i for i in reg if cellexpr[i] >  thresholdExpr]
		return sorted(reg)			
	def getActivators(self, target = None, cell = None, filt = True, thresholdExpr = 0):
		return self.getRegulators(target, cell, filt, thresholdExpr, [self.ACTIVATOR_TYPE])
	def getInhibitors(self, target = None, cell = None, filt = True, thresholdExpr = 0):
		return self.getRegulators(target, cell, filt, thresholdExpr, [self.INHIBITOR_TYPE])
	def getRegulatorsExpr(self, target = None, cell = None, filt = True, thresholdExpr = 0, types = ['1', '-1']):
		reg = self.getRegulators(target, cell, filt, thresholdExpr, types)
		exp = self.expr.loc[reg]
		exp = exp[[i for i in exp.columns if 'exp' in i]]
		if(cell is not None):
			exp = exp[[i for i in exp.columns if str(cell) in i]]
		exp['tf_ind'] = reg
		return exp
	def getActivatorsExpr(self, target = None, cell = None, filt = True, thresholdExpr = 0):
		return self.getRegulatorsExpr(target, cell, filt, thresholdExpr, [self.ACTIVATOR_TYPE])
	def getInhibitorsExpr(self, target = None, cell = None, filt = True, thresholdExpr = 0):
		return self.getRegulatorsExpr(target, cell, filt, thresholdExpr, [self.INHIBITOR_TYPE])
	def isActivator(self, gene):
		return True if (self.expr.type[gene] == self.ACTIVATOR_TYPE) else False
	def isInhibitor(self, gene):
		return True if (self.expr.type[gene] == self.INHIBITOR_TYPE) else False
	def __getSubExpr(self, cell = None, gene = None, sub='exp'):
		e = self.expr[[i for i in self.expr.columns if sub in i]]
		if(cell is not None):
			e = e[[i for i in e.columns if str(cell) in i]]
		if(gene is not None):
			e = e.loc[gene]
		return e
	def getOptimalExpr(self, cell = None, gene = None):
		return self.__getSubExpr(cell, gene, 'opt')
	def getInitialExpr(self, cell = None, gene = None):
		return self.__getSubExpr(cell, gene, 'init')
	def getFinalExpr(self, cell = None, gene = None):
		return self.__getSubExpr(cell, gene, 'exp')
	def getSetRegulators(self, genelist = None, cell = None, filt = True, thresholdExpr= 0, types = ['1', '-1']):
		if(genelist is not None):
			l = []
			for g in genelist:
				l = l + self.getRegulators(target=g, cell=cell, filt=filt, thresholdExpr=thresholdExpr, types=types)
			return sorted(list(set(l)))
		else:
			return sorted(self.getRegulators(target=None, cell=cell,filt=filt, thresholdExpr=thresholdExpr, types =  ['1', '-1']))
	def getSetActivators(self, genelist = None, cell = None, filt = True, thresholdExpr= 0):
		return self.getSetRegulators(genelist, cell, filt, thresholdExpr, ['1'])
		
	def getSetInhibitors(self, genelist = None, cell = None, filt = True, thresholdExpr= 0):
		return self.getSetRegulators(genelist, cell, filt, thresholdExpr, ['-1'])
	def getIntersectionRegulators(self, genelist = None, cell = None, filt = True, thresholdExpr= 0, types = ['1', '-1']):
		if(genelist is not None):
			l = self.getRegulators(target=genelist[0], cell=cell, filt=filt, thresholdExpr=thresholdExpr, types=types)
			for g in genelist[1:]:
				aux = self.getRegulators(target=g, cell=cell, filt=filt, thresholdExpr=thresholdExpr, types=types)
				l = [i for i in l if (i in aux)]
			return sorted(list(set(l)))
		else:
			return sorted(self.getRegulators(target=None, cell=cell,filt=filt, thresholdExpr=thresholdExpr, types =  ['1', '-1']))
	def getIntersectionActivators(self, genelist = None, cell = None, filt = True, thresholdExpr= 0):
		return self.getIntersectionRegulators(genelist, cell, filt, thresholdExpr, ['1'])		
	def getIntersectionInhibitors(self, genelist = None, cell = None, filt = True, thresholdExpr= 0):
		return self.getIntersectionRegulators(genelist, cell, filt, thresholdExpr, ['-1'])
	def getSortedExpressionByMean(self):
		pass
	def getSortedExpresionEach(self):
		pass
	def getSortedMutEffectsByMean(self, type1 = None, type2=None, type3=['[1]'], regType1=['1', '-1'], regType2=None, regType3=None, sortByCell=False, regList=None, mutation_type = 'tfs', difference_to_wt=True, error=True, type_used_sort = 1, sort_gene_in_each_cell = True):
		#get targets of required types
		targets = []
		if(type1 is not None):
			for i in range(self.genome_size):
				targets = targets + [i for (i, j) in enumerate(self.type) if j in type1]
		if(type2 is not None):
			for i in range(self.genome_size):
				targets = targets + [i for (i, j) in enumerate(self.type2) if j in type2]
		if(type3 is not None):
			for i in range(self.genome_size):
				targets = targets + [i for (i, j) in enumerate(self.type3) if j in type3]
		targets = sorted(set(targets))
		#get regulators of required types
		regulators = []
		if(regType1 is not None):
			for i in range(self.genome_size):
				regulators = regulators + [i for (i, j) in enumerate(self.type) if j in regType1]
		if(regType2 is not None):
			for i in range(self.genome_size):
				regulators = regulators + [i for (i, j) in enumerate(self.type2) if j in regType2]
		if(regType3 is not None):
			for i in range(self.genome_size):
				regulators = regulators + [i for (i, j) in enumerate(self.type3) if j in regType3]
		regulators = sorted(set(regulators))
		if(mutation_type == 'tfs'):
			mat = self.mutanttfs
			vname = 'tf_mutated'
		else:
			mat = self.mutantsites
			vname = 'mutated_sites'
		mat = mat.iloc[[i for i in range(len(mat.index)) if mat.index[i] in targets and mat[vname].iloc[i] in regulators]]
		indTar = list(mat.index)
		indMut = list(mat[vname])
		exp = np.array(mat[[i for i in mat.columns if 'exp' in i]])
		opt = np.array(mat[[i for i in mat.columns if 'opt' in i]])
		#init = np.array(mat[[i for i in mat.columns if 'init' in i]])
		if(error):
			exp = exp - opt
			if(difference_to_wt):
				wt = np.array(self.expr[[i for i in self.expr.columns if 'exp' in i]])[indTar, :] - opt
				exp = exp - wt
		elif(difference_to_wt):
				exp = exp - np.array(self.expr[[i for i in self.expr.columns if 'exp' in i]])[indTar, :]			
		## Sort exp by phenotype. 
		## Sort each cell separately
		if(type_used_sort < 3):
			typesort = np.array(self.type)
		else:
			typesort = np.array(self.type3)
		output = np.zeros([len(regulators), exp.shape[1], len(targets)])
		for this_reg_ind, this_reg in enumerate(regulators):
			cellorder = np.argsort([np.mean(exp[[j for j in range(len(indMut)) if indMut[j] == this_reg], :][:,k]) for k in range(exp.shape[1])])
			gene_ind = 0
			for this_type in np.unique(typesort[targets]):
				if(sortByCell):
					cellorder = np.argsort([np.mean(exp[[j for j in range(len(indMut)) if indMut[j] == this_reg and typesort[indTar[j]] == this_type], :][:,k]) for k in range(exp.shape[1])])
				for cell_ind, cell in enumerate(cellorder):
					this_type_genes = exp[[j for j in range(len(indMut)) if indMut[j] == this_reg and typesort[indTar[j]] == this_type ], cell]
					if(cell_ind == 0):
						gene_inds_order = np.argsort(this_type_genes)
						this_type_genes = np.sort(this_type_genes)
					elif(sort_gene_in_each_cell):
						this_type_genes = np.sort(this_type_genes)
					else:
						this_type_genes = this_type_genes[gene_inds_order]
					output[this_reg_ind, cell_ind, gene_ind:(gene_ind + this_type_genes.shape[0])] = this_type_genes
				gene_ind = gene_ind + this_type_genes.shape[0]	
		output= output.reshape((output.shape[0], output.shape[1]*output.shape[2]))
		ft = np.array(self.type3[regulators])
		return (output, ft)
class simulationSet(object):
	def __init__(self, path, cond):
		if(path.endswith('/')):
			path = path[:-1] 
		self.path = path
		self.cond = cond
		simnames = list(set(['_'.join(f.split('_')[:-1]) for f in os.listdir(path) if f.endswith('.csv')]))
		simnames = list(set([i+'_' for i in simnames if(not i.endswith('_'))]))
		self.simulations = [simulationHolder(s, path) for s in simnames]
		self.stored_data = []
	def joinRegExpression(self, types = ['1', '-1'], shape2d = True):
		if(shape2d):
			ax = 0
		else:
			ax = 2
		m = np.array(self.simulations[0].getRegulatorsExpr(types = types))
		dim0 = m.shape[0]
		dim1 = m.shape[1]
		m = m.reshape((dim0, dim1, 1))
		sim_types = np.array(self.simulations[0].type3)[self.simulations[0].getRegulators(filt = False, types = types)].reshape((dim0))
		for i in range(1, len(self.simulations)):
			m = np.concatenate((m, np.array(self.simulations[i].getRegulatorsExpr(types = types)).reshape((dim0, dim1, 1)) ), axis=ax)
			sim_types = np.concatenate((sim_types, np.array(self.simulations[i].type3)[self.simulations[i].getRegulators(filt = False, types = types)].reshape((dim0))), axis=0)
		if(shape2d):
			m=m.reshape((m.shape[0], m.shape[1]))
		return (sim_types, m)
	def joinMutEffects(self, type1 = None, type2=None, type3=['[1]'], regType1=['1', '-1'], regType2=None, regType3=None, sortByCell=False, regList=None, mutation_type = 'tfs', difference_to_wt=True, error=True, type_used_sort = 1, sort_gene_in_each_cell=True):
		m, tftype = self.simulations[0].getSortedMutEffectsByMean(type1, type2, type3, regType1, regType2, regType3, sortByCell, regList, mutation_type, difference_to_wt, error, type_used_sort, sort_gene_in_each_cell)
		dim0 = m.shape[0]
		dim1 = m.shape[1]
		ax = 0
		for i in range(1, len(self.simulations)):
			aux1, aux4 = self.simulations[i].getSortedMutEffectsByMean(type1, type2, type3, regType1, regType2, regType3, sortByCell, regList, mutation_type, difference_to_wt, error, type_used_sort, sort_gene_in_each_cell)
			m = np.concatenate((m, aux1), axis=ax)
			tftype = np.concatenate((tftype, aux4), axis=0)
		return (m, tftype)
	def makeHeatmaps(self):
		import copy
		tabcols = ['tab:red', 'tab:blue', 'tab:green', 'tab:orange', 'yellow', 'wheat', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan', 'darkseagreen','crimson' ]
		tabcols2 = ['tab:red', 'tab:blue', 'slateblue', 'tab:orange', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan', 'darkseagreen','crimson' ]
		type3_cmap = dict(zip(list(np.unique(self.simulations[0].type3)), tabcols[0:len(list(np.unique(self.simulations[0].type3)))]))
		type_cmap = dict(zip(list(np.unique(self.simulations[0].type)), tabcols2[0:len(list(np.unique(self.simulations[0].type)))]))
		ncells = self.simulations[0].ncells
		if(np.any(self.simulations[0].type == '[0]')):
			cell_colstart = np.where(np.unique(self.simulations[0].type) == '[0]')[0][0]
		else:
			cell_colstart = np.where(np.unique(self.simulations[0].type) == '[01]')[0][0]
		cellmap = dict(zip([i for i in range(ncells)], tabcols2[cell_colstart:cell_colstart+ncells]))
		t0, x0 = self.joinRegExpression(types = ['1'])
		t0oth = [i for i in t0 if i =='1lin[0]']
		t0lin = [i for i in t0 if i !='1lin[0]']
		x0oth = x0[[i for i in range(len(t0)) if t0[i] =='1lin[0]'],:]
		x0lin = x0[[i for i in range(len(t0)) if t0[i] !='1lin[0]'],:]
		t1, x1 = self.joinRegExpression(types = ['-1'])
		t2, x2 = self.joinRegExpression(types = ['1', '-1'])
		rowcol0 = [type3_cmap[i] for i in t0]
		rowcol0oth = [type3_cmap[i] for i in t0oth]
		rowcol0lin = [type3_cmap[i] for i in t0lin]
		rowcol1 = [type3_cmap[i] for i in t1]
		rowcol2 = [type3_cmap[i] for i in t2]
		colcol = []
		for c in range(ncells):
			colcol = colcol + [cellmap[c] for i in range(int(x0[:, :-1].shape[1]/ncells))]
		g0 = sns.clustermap(np.sort(x0[:,:-1], axis=1)[:, ::-1], col_cluster=False,  row_colors=rowcol0,col_colors=colcol, method='ward', center = 0, vmin = 0, vmax = 4.5)
		g0oth = sns.clustermap(np.sort(x0oth[:,:-1], axis=1)[:, ::-1], col_cluster=False,  row_colors=rowcol0oth,col_colors=colcol, method='ward', center = 0, vmin = 0, vmax = 4.5)
		g0lin = sns.clustermap(np.sort(x0lin[:,:-1], axis=1)[:, ::-1], col_cluster=False,  row_colors=rowcol0lin,col_colors=colcol, method='ward', center = 0, vmin = 0, vmax = 4.5)
		g1 = sns.clustermap(np.sort(x1[:,:-1], axis=1)[:, ::-1], col_cluster=False,  row_colors=rowcol1, col_colors=colcol, method='ward', center = 0, vmin = 0, vmax = 4.5)
		g2 = sns.clustermap(np.sort(x2[:,:-1], axis=1)[:, ::-1], col_cluster=False,  row_colors=rowcol2, col_colors=colcol, method='ward', center = 0, vmin = 0, vmax = 4.5)
		g0.ax_heatmap.set_xticklabels([''])
		g0oth.ax_heatmap.set_xticklabels([''])
		g0lin.ax_heatmap.set_xticklabels([''])
		g1.ax_heatmap.set_xticklabels([''])
		g2.ax_heatmap.set_xticklabels([''])
		g0.ax_heatmap.set_yticklabels([''])
		g0oth.ax_heatmap.set_yticklabels([''])
		g0lin.ax_heatmap.set_yticklabels([''])
		g1.ax_heatmap.set_yticklabels([''])
		g2.ax_heatmap.set_yticklabels([''])
		g0.ax_heatmap.set_ylabel('TF', fontsize = 15)
		g0oth.ax_heatmap.set_ylabel('TF', fontsize = 15)
		g0lin.ax_heatmap.set_ylabel('TF', fontsize = 15)
		g1.ax_heatmap.set_ylabel('TF', fontsize = 15)
		g2.ax_heatmap.set_ylabel('TF', fontsize = 15)
		g0.ax_heatmap.set_xlabel('cell', fontsize = 15)
		g0oth.ax_heatmap.set_xlabel('cell', fontsize = 15)
		g0lin.ax_heatmap.set_xlabel('cell', fontsize = 15)
		g1.ax_heatmap.set_xlabel('cell', fontsize = 15)
		g2.ax_heatmap.set_xlabel('cell', fontsize = 15)
		#g0.ax_heatmap.set_title("Activator Expression", fontsize = 15)
		#g1.ax_heatmap.set_title("Inhibitor Expression", fontsize = 15)
		#g2.ax_heatmap.set_title("All TF xpression", fontsize = 15)
		g0.savefig(self.path + '_activatorExpression.png')
		g0oth.savefig(self.path + '_activatorOthExpression.png')
		g0lin.savefig(self.path + '_activatorLinExpression.png')
		g1.savefig(self.path + '_inhibitorExpression.png')
		g2.savefig(self.path + '_allTFExpression.png')
		plt.close("all")
		## TFS /sites
		## TFS /sites
		x3={}
		for muttype in ['tfs', 'sites']:
			x3[muttype]={}
			for i in np.unique(self.simulations[0].type3):
				if(muttype == 'sites' and (i[0] == '-' or i[0]=='1')):
					continue
				these_types=[self.simulations[0].type.tolist()[j] for j in range(len(self.simulations[0].type.tolist())) if self.simulations[0].type3.tolist()[j] == i]
				#xaux, taux = self.joinMutEffects(type1 = None, type2=None, type3=[i], regType1=['1', '-1'], regType2=None, regType3=None, sortByCell=False, regList=None, mutation_type = muttype, difference_to_wt=True, error=True, type_used_sort = 3, sort_gene_in_each_cell=True) ### use this to plot also inhibitor efects
				xaux, taux = self.joinMutEffects(type1 = None, type2=None, type3=[i], regType1=['1'], regType2=None, regType3=None, sortByCell=False, regList=None, mutation_type = muttype, difference_to_wt=True, error=True, type_used_sort = 3, sort_gene_in_each_cell=True)
				rowcol3 = [type3_cmap[i] for i in taux]
				colcol3 = [[],[]]
				for c in range(ncells):
					colcol3[0] = colcol3[0] + [cellmap[c] for j in range(int(xaux.shape[1]/ncells))]
					if('lin' in i):				#if it's a TF, use type 3 to differentiate between lineage etc
						colcol3[1] = colcol3[1] + [type3_cmap[i] for j in range(int(xaux.shape[1]/ncells))]
					else:					#if it's terminal feature, use type1 to differentiate between features of different cells
						colcol3[1] = colcol3[1] + [type_cmap[these_types[j]] for j in range(int(xaux.shape[1]/ncells))]
				gaux = sns.clustermap(xaux, col_cluster=False,  col_colors=colcol3,row_colors=rowcol3, method='ward', center = 0, vmin = -3.5, vmax = 3.5)
				gaux.ax_heatmap.set_xticklabels([''])
				gaux.ax_heatmap.set_yticklabels([''])
				gaux.ax_heatmap.set_xlabel(muttype + 'effect in '+ i +' per cell', fontsize = 15)
				gaux.ax_heatmap.set_ylabel('TF', fontsize = 15)
				#gaux.ax_heatmap.set_title("Activator Expression", fontsize = 15)
				gaux.savefig(self.path + '_' + muttype + 'Phenotype_' + i + '.png')
				#Join with expression 
				colcol3e = []
				colcol3e.append(['dimgrey' for i in range(len(colcol))] + ['k' for i in colcol3[1]]) 
				colcol3e.append(colcol + colcol3[0])
				colcol3e.append(['w' for i in range(len(colcol))] + colcol3[1]) 
				#x3e = np.concatenate((np.sort(x2[:,:-1], axis=1)[:, ::-1], xaux), axis=1) #this to use inhibitors too
				x3e = np.concatenate((np.sort(x0[:,:-1], axis=1)[:, ::-1], xaux), axis=1)
				gauxe = sns.clustermap(x3e, col_cluster=False,  col_colors=colcol3e,row_colors=rowcol3, method='ward', center = 0, vmin = -3.5, vmax = 3.5)
				gauxe.ax_heatmap.set_xticklabels([''])
				gauxe.ax_heatmap.set_yticklabels([''])
				gauxe.ax_heatmap.set_xlabel(muttype + 'effect in '+ i +' per cell', fontsize = 15)
				gauxe.ax_heatmap.set_ylabel('TF', fontsize = 15)
				#gauxe.ax_heatmap.set_title("Activator Expression", fontsize = 15)
				gauxe.savefig(self.path + '_' + muttype + 'Phenotype_' + i + '_andEXPR.png')
				x3[muttype][i] = [xaux, taux, gaux, gauxe, colcol3]
				## same thing, cell1
				colcol4e = []
				for lev in range(len(colcol3e)):
					colcol4e.append([n for k, n in enumerate(colcol3e[lev]) if colcol3e[1][k] == cellmap[0] or colcol3e[2][k] == 'w'])
				x4e = x3e[:, [k for k in range(len(colcol3e[0])) if cellmap[0] == colcol3e[1][k] or colcol3e[2][k]== 'w']]
				gauxe = sns.clustermap(x4e, col_cluster=False,  col_colors=colcol4e,row_colors=rowcol3, method='ward', center = 0, vmin = -3.5, vmax = 3.5)
				gauxe.ax_heatmap.set_xticklabels([''])
				gauxe.ax_heatmap.set_yticklabels([''])
				gauxe.ax_heatmap.set_xlabel(muttype + 'effect in '+ i +' per cell', fontsize = 15)
				gauxe.ax_heatmap.set_ylabel('TF', fontsize = 15)
				gauxe.savefig(self.path + '_' + muttype + 'Phenotype_' + i + '_andEXPR_cell0.png')
				plt.close("all")
			## combinations		
			tojoin = [j for j in np.unique(self.simulations[0].type3) if 'lin' not in j]
			#x5 = np.sort(x2[:,:-1], axis=1)[:, ::-1] # with inhibitors
			x5 = np.sort(x0[:,:-1], axis=1)[:, ::-1]
			t5 = x3[muttype][tojoin[0]][1] #same for all
			rowcol5 = [type3_cmap[i] for i in t5]
			colcol5 = []
			colcol5.append(colcol)
			colcol5.append(['dimgrey' for i in range(len(colcol))])
			colcol5.append(['w' for i in range(len(colcol))])
			for j in range(len(tojoin)):
				x5 = np.concatenate((x5, x3[muttype][tojoin[j]][0]), axis=1)
				colcol5[0] = colcol5[0] + x3[muttype][tojoin[j]][4][0]	
				colcol5[1] = colcol5[1] + x3[muttype][tojoin[j]][4][1]	
				colcol5[2] = colcol5[2] + [type3_cmap[tojoin[j]] for i in range(len(x3[muttype][tojoin[j]][4][0]))]
			g5 = sns.clustermap(x5, col_cluster=False,  col_colors=colcol5,row_colors=rowcol5, method='ward', center = 0, vmin = -3.5, vmax = 3.5)
			g5.ax_heatmap.set_xticklabels([''])
			g5.ax_heatmap.set_yticklabels([''])
			g5.ax_heatmap.set_xlabel(muttype + 'effect in all terminal features per cell', fontsize = 15)
			g5.ax_heatmap.set_ylabel('TF', fontsize = 15)
			g5.savefig(self.path + '_' + muttype + 'Phenotype_TerminalJoinAndExpr.png')
			# combinations only in cell 0
			colcol6 = []
			for lev in range(len(colcol5)):
				colcol6.append([n for k, n in enumerate(colcol5[lev]) if colcol5[0][k] == cellmap[0] or colcol5[2][k] == 'w'])
			x6 = x5[:, [k for k in range(len(colcol5[1])) if cellmap[0] == colcol5[0][k] or colcol5[2][k]== 'w']]
			g6 = sns.clustermap(x6, col_cluster=False,  col_colors=colcol6,row_colors=rowcol5, method='ward', center = 0, vmin = -3.5, vmax = 3.5)
			g6.ax_heatmap.set_xticklabels([''])
			g6.ax_heatmap.set_yticklabels([''])
			g6.ax_heatmap.set_xlabel(muttype + 'effect in all terminal features per cell', fontsize = 15)
			g6.ax_heatmap.set_ylabel('TF', fontsize = 15)
			g6.savefig(self.path + '_' + muttype + 'Phenotype_TerminalJoinAndExpr_cell0.png')
			#finally, try to plot only terminal_this:
			cellstr=[str(c) for c in range(ncells)]
			try: #cutre, solo servira para mce0 y asi
				relevant_types=[t for t in np.unique(self.simulations[0].type3) if cellstr[0] in t]
				x7 = x6[:, [i for i in range(x6.shape[1]) if (i-ncells<5) or (i-ncells<25 and i-ncells>= 20)]]
				colcol7=[c for c, i in zip(colcol6, range(x6.shape[1])) if (i-ncells<5) or (i-ncells<25 and i-ncells>= 20)]
				g7 = sns.clustermap(x7, col_cluster=False,  col_colors=colcol7,row_colors=rowcol5, method='ward', center = 0, vmin = -3.5, vmax = 3.5)
				g7.ax_heatmap.set_xticklabels([''])
				g7.ax_heatmap.set_yticklabels([''])
				g7.ax_heatmap.set_xlabel(muttype + 'effect in terminal features this cell', fontsize = 15)
				g7.ax_heatmap.set_ylabel('TF', fontsize = 15)
				g7.savefig(self.path + '_' + muttype + 'Phenotype_TerminalJoinAndExpr_cell0selg.png')
			except:
				print("tried to hetmap bad list of genes")
		self.stored_data.append(x3)
		plt.close("all")


def main():
	#conditions = ['4cell_mce0', '4cell_mce0fix_mutb', '4cell_mce1fix', '4cell_mce2fix', '4cell_mce2inh5', '4cell_mce0fix', '4cell_mce0_mutb', '4cell_mce1fix_mutb', '4cell_mce2fix_mutb', '4cell_mce2inh5_mutb']
	conditions = ['4cell_mce11fix']
	#paths= ['/home/cmora/cellevolver_shared/simul_tables/' + i for i in conditions]
	#conditions = ['4cell_mce0Xss', '4cell_mce1inh5Xss', '4cell_mce2inh5Xss']
	#paths= ['../cellevolver_shared/cellevolver_chromatin/simulations/all_tables_chromatin/' + i for i in conditions]
	#conditions = ['4cell_mce4fix', '4cell_mce5fix', '4cell_mce8fix', '4cell_mce9fix', '6cell_mce6fix', '6cell_mce7fix']
	paths= ['./sim_tables/' + i for i in conditions]
	for cond, path in zip(conditions, paths):
		print(path, '\n')
		ss=simulationSet(path, cond)	
		print('read...')
		ss.makeHeatmaps()



if __name__ == "__main__":
    main()

