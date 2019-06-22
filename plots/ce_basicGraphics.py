import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns; sns.set()
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import gridspec


def readData(pckname, filetype):
	f = os.listdir()
	datafile = [i for i in f if i.endswith('.csv') and pckname in i and filetype in i]
	assert len(datafile) == 1, "More than one possible file: " + pckname + ' - ' + filetype
	datafile = datafile[0]
	data = pd.read_csv(datafile, sep = "\t", header=0, index_col=0)
	return data

def finalExpressionHeatmap(pckname, filetype = 'finalExpression'):
	open_pdf = PdfPages(pckname + '_' + filetype + '.pdf')
	d = readData(pckname, filetype)
	t2 = []
	for i in range(len(d.axes[0])):
		if (d['type'][i] == '-1'):
			t2.append(4)
		elif (d['type'][i] == '1'):
			t2.append(3)
		else:
			t2.append(0)
	d['color'] = t2
	sns.set_context("paper")
	paleta = sns.color_palette("Paired")
	#cell names
	names=['cell ' + str(i) for i in range(len([i for i in d.columns if 'exp' in i]))]
	#fig, ax = plt.subplots(figsize=(15,10), ncols=4, nrows=1)
	fig = plt.figure(figsize=(15,10)) #figure
	gs = gridspec.GridSpec(1, 4, width_ratios = [0.2, 1, 1, 1.25]) #needed for different sizes
	left   =  0.125  # the left side of the subplots of the figure
	right  =  0.9    # the right side of the subplots of the figure
	bottom =  0.1    # the bottom of the subplots of the figure
	top    =  0.9    # the top of the subplots of the figure
	wspace =  .5     # the amount of width reserved for blank space between subplots
	hspace =  1.1    # the amount of height reserved for white space between subplots
	# This function actually adjusts the sub plots using the above paramters
	plt.subplots_adjust(
		left    =  left, 
		bottom  =  bottom, 
		right   =  right, 
		top     =  top, 
		wspace  =  wspace, 
		hspace  =  hspace
		)
	# The amount of space above titles
	y_title_margin = 1.2
	#plt.suptitle("Simulation 11", y = 1.09, fontsize=20)
	ax = [plt.subplot(gs[i]) for i in range(4)]
	ax[0].set_title("Gene type", fontsize = 10)
	ax[1].set_title("Initial expression", fontsize = 15)
	ax[2].set_title("Optimal expression (t=10)", fontsize = 15)
	ax[3].set_title("Expression (t=10)", fontsize = 15)
	sns.heatmap(d[['color']], linewidths=.1, cbar = False, cmap=paleta,center = 3, vmin = 0, vmax = 4, yticklabels = 1, ax=ax[0])
	sns.heatmap(d[[i for i in d.columns if 'init' in i]], linewidths=.1, cbar = False, cmap="Blues",center = 2, vmin = 0, vmax = 4, yticklabels = 1, ax = ax[1])
	sns.heatmap(d[[i for i in d.columns if 'opt' in i]], linewidths=.1, cbar = False, cmap="Blues",center = 2, vmin = 0, vmax = 4, yticklabels = 1, ax=ax[2])
	sns.heatmap(d[[i for i in d.columns if 'exp' in i]], linewidths=.1, cbar = True, cmap="Blues",center = 2, vmin = 0, vmax = 4, yticklabels = 1, ax=ax[3])
	# Set all labels on the row axis of subplots for square_feet data to "square_feet"
	[ax[i].set_xticklabels(names, fontsize=10) for i in range(1, 4)]
	[ax[i].set_yticklabels([str(i) for i in d.axes[0]],rotation = 0, fontsize=8) for i in range(0, 4)]
	ax[0].set_ylabel('gene', fontsize = 15)
	ax[0].set_xticklabels([''])
	open_pdf.savefig()
	open_pdf.close()


def mutantExpressionHeatmap(pckname, filetype = ['finalExpression', 'mutantTFs', 'mutantSites']):
	open_pdf = PdfPages(pckname + '_mutants.pdf')
	wt = readData(pckname, filetype[0])
	tfmut = readData(pckname, filetype[1])
	sitemut = readData(pckname, filetype[2])
	tf_names = sorted(tfmut['tf_mutated'].unique())
	t2 = []
	for i in range(len(wt.axes[0])):
		if (wt['type'][i] == '-1'):
			t2.append(4)
		elif (wt['type'][i] == '1'):
			t2.append(3)
		else:
			t2.append(0)
	wt['color'] = t2
	sns.set_context("paper")
	paleta = sns.color_palette("Paired")
	#cell names
	names=['cell ' + str(i) for i in range(len([i for i in wt.columns if 'exp' in i]))]
	for tf in tf_names:
	#fig, ax = plt.subplots(figsize=(15,10), ncols=4, nrows=1)
		fig = plt.figure(figsize=(15,10)) #figure
		gs = gridspec.GridSpec(1, 5, width_ratios = [0.1, 1, 1.25, 1, 1.25]) #needed for different sizes
		left   =  0.125  # the left side of the subplots of the figure
		right  =  0.9    # the right side of the subplots of the figure
		bottom =  0.1    # the bottom of the subplots of the figure
		top    =  0.9    # the top of the subplots of the figure
		wspace =  .5     # the amount of width reserved for blank space between subplots
		hspace =  1.1    # the amount of height reserved for white space between subplots
		# This function actually adjusts the sub plots using the above paramters
		plt.subplots_adjust(
			left    =  left, 
			bottom  =  bottom, 
			right   =  right, 
			top     =  top, 
			wspace  =  wspace, 
			hspace  =  hspace
			)
		# The amount of space above titles
		y_title_margin = 1.2
		#plt.suptitle("Simulation 11", y = 1.09, fontsize=20)
		ax = [plt.subplot(gs[i]) for i in range(5)]
		ax[0].set_title("Gene type", fontsize = 6)
		ax[1].set_title("TF " + str(tf) + " complete ko - expression", fontsize = 10)
		ax[2].set_title("TF " + str(tf) + " sites in non-TF genes ko - expresion", fontsize = 10)
		ax[3].set_title("TF " + str(tf) + " complete ko - difference", fontsize = 10)
		ax[4].set_title("TF " + str(tf) + " sites in non-TF genes ko - difference", fontsize = 10)
		sns.heatmap(wt[['color']], linewidths=.1, cbar = False, cmap=paleta,center = 3, vmin = 0, vmax = 4, yticklabels = 1, ax=ax[0])
		sns.heatmap(tfmut[[i for i in tfmut.columns if 'exp' in i]][tfmut['tf_mutated'] == tf], linewidths=.1, cbar = False, cmap="Blues",center = 2, vmin = 0, vmax = 4, yticklabels = 1, ax = ax[1])
		sns.heatmap(sitemut[[i for i in wt.columns if 'exp' in i]][sitemut['mutated_sites'] == tf], linewidths=.1, cbar = True, cmap="Blues",center = 2, vmin = 0, vmax = 4, yticklabels = 1, ax=ax[2])
		sns.heatmap(tfmut[[i for i in tfmut.columns if 'exp' in i]][tfmut['tf_mutated'] == tf] - wt[[i for i in wt.columns if 'exp' in i]], linewidths=.1, cbar = False, cmap="PuOr",center = 0, vmin = -2, vmax = 1, yticklabels = 1, ax = ax[3])
		
		sns.heatmap(sitemut[[i for i in sitemut.columns if 'exp' in i]][sitemut['mutated_sites'] == tf] - wt[[i for i in wt.columns if 'exp' in i]], linewidths=.1, cbar = True, cmap="PuOr",center = 0, vmin = -2, vmax = 1, yticklabels = 1, ax = ax[4])
		# Set all labels on the row axis of subplots for square_feet data to "square_feet"
		[ax[i].set_xticklabels(names, fontsize=6) for i in range(1, 5)]
		[ax[i].set_yticklabels([str(i) for i in wt.axes[0]],rotation = 0, fontsize=8) for i in range(0, 5)]
		ax[0].set_ylabel('gene', fontsize = 15)
		ax[0].set_xticklabels([''])
		open_pdf.savefig()
	open_pdf.close()

def timeExpression(pckname, colormap = None, filetype = ['finalExpression', 'timeExpression']):
	#One color for each gene and one linetype for each gene type (inhibitor, activator, or terminal feature with specific optimal expression pattern)
	wt = readData(pckname, filetype[0])
	t = readData(pckname, filetype[1])
	plt.figure(figsize=(15,10))
	names = [i for i in t.columns if 'cell' in i]
	tf_num = len([i for i in wt['type'] if '1' ==i or '-1' == i])
	####Colors
	from cycler import cycler
	if(colormap is None):
		colormap = plt.cm.gist_ncar	
	colormap1 = colormap(np.linspace(0,1,tf_num))
	colormap2 = colormap(np.linspace(0,1,len(wt['type']) - tf_num))
	####linestyles
	linestyles = ['-', '--', '-.', ':']
	linestyleDict = {}
	types = wt['type'].unique()
	for i in range(len(types)):
		linestyleDict[types[i]] = linestyles[i%len(linestyles)]
	linestyleVec = [linestyleDict[i] for i in wt['type']]
	#plot
	g=1
	for cell in range(len(names)): 
		ax = plt.subplot(len(names),2,g)
		ax.set_prop_cycle(cycler('color', colormap1) + cycler('linestyle', linestyleVec[:tf_num]))
		if(cell == 0):
			plt.title('Transcription Factors expression', fontsize=15)
		h1 = [plt.plot(t['time'][i], t[names[cell]][i]) for i in range(tf_num)]
		plt.ylabel('cell ' + str(cell), fontsize=15)
		if (cell == len(names)-1):
			plt.xlabel('time', fontsize=15)
		if(cell == 0):
			ax.legend(labels = list(range(0, tf_num)), loc='lower center', bbox_to_anchor=(-0.18, -1.33),ncol=1, fancybox=True, shadow=True, prop={'size': 10})
		g+=1
		ax = plt.subplot(len(names),2,g)
		ax.set_prop_cycle(cycler('color', colormap2) + cycler('linestyle', linestyleVec[tf_num:]))
		if(cell == 0):
			plt.title('Terminal Features expression', fontsize=15)
		h2 = [plt.plot(t['time'][i], t[names[cell]][i]) for i in range(tf_num, len(wt['type']))]
		if (cell == len(names)-1):
			plt.xlabel('time', fontsize=15)
		if(cell == 0):		
			ax.legend(labels = list(range(tf_num, len(wt['type']))), loc='lower center', bbox_to_anchor=(1.1, -2.8),ncol=1, fancybox=True, shadow=True, prop={'size': 10})
		g+=1
	open_pdf = PdfPages(pckname + '_timeExpression.pdf')
	open_pdf.savefig()
	open_pdf.close()
	#plt.show()

def motifGenomeDiagram(pckname, filetype = ['finalExpression', 'annotation', 'TFset']):
	from reportlab.lib import colors
	from reportlab.lib.units import cm
	from Bio.Graphics import GenomeDiagram
	from Bio import SeqIO
	from Bio.SeqFeature import SeqFeature, FeatureLocation
	wt = readData(pckname, filetype[0])
	d = readData(pckname, filetype[1])
	d['strand'] = [-1 if i==0 else 1 for i in d.strand]
	tfset = readData(pckname, filetype[2])
	cols = list(colors.getAllNamedColors().keys())
	cols = [col for col in cols if 'pale' not in col and 'light' not in col and 'white' not in col and 'snow' not in col and 'ivory' not in col]
	cols = np.random.choice(cols,  len(tfset.tf_name.unique()))
	gdd = GenomeDiagram.Diagram(pckname)
	for gene in d.gene.unique():
		this_gene = d[:][d.gene == gene]
		gdt_features = gdd.new_track(len(d.gene.unique()) - gene, greytrack=False)
		gds_features = gdt_features.new_set()
		for motif in this_gene.axes[0]:
			if(this_gene.strand[motif] == 1):
				angle = 0
			else:
				angle = 180
			mot_size = len(tfset.consensus[tfset.tf_name == this_gene.tf_ind[motif]][0])
			feature = SeqFeature(FeatureLocation(int(this_gene.positions[motif]), int(this_gene.positions[motif] + mot_size - 1)), strand=this_gene.strand[motif])
			gds_features.add_feature(feature, name=str(this_gene.tf_ind[motif]), label=True, label_position="start",
                                   label_size = 6, label_angle=angle, color = cols[this_gene.tf_ind[motif]], sigil="BIGARROW", arrowshaft_height=this_gene.percent[motif],  arrowhead_length=1)
	gdd.draw(format='linear', pagesize=(15*cm, 10*cm), fragments=1,start=0, end=len(wt.seqs[0]), orientation = 'portrait', tracklines = 0)
	gdd.write(pckname +"_motifs.pdf", "pdf")
def regulationHeatmap(pckname, filetype = ['finalExpression', 'annotation', 'TFset']):
	open_pdf = PdfPages(pckname + 'RegulationHeatmap.pdf')
	wt = readData(pckname, filetype[0])
	d = readData(pckname, filetype[1])
	reg = pd.DataFrame()
	reg['gene'] = sorted(d.gene.unique())
	reg['type'] = wt.type
	tf_names = sorted(d.tf_ind.unique())
	#generate interaction table
	for tf in tf_names:
		aux = d[:][d.tf_ind == tf]
		reg['Complete' + str(tf)] = [sum(aux.percent[aux.gene == i]) for i in reg.gene]
		if(reg['type'][tf] == '-1'):
			reg['Complete' + str(tf)] = -1*reg['Complete' + str(tf)]
		for c in range(len([i for i in wt.columns if 'exp' in i])):
			newname = str(tf) + 'Cell ' + str(c)
			reg[newname] = reg['Complete' + str(tf)]*wt['exp'+str(c)][tf]
	reg.to_csv(pckname +'_RegulationTable.csv', sep='\t', header=True,decimal='.', float_format='%.10f')
	#Colors
	colortypes = {}
	for i in range(len(reg.type.unique())):
		colortypes[reg.type.unique()[i]] = i
	t2 = [colortypes[t] for t in reg.type]
	reg['color'] = t2
	sns.set_context("paper")
	paleta = sns.color_palette("Paired")
	#cell names
	names=['Cell ' + str(i) for i in range(len([i for i in wt.columns if 'exp' in i]))]
	names.insert(0,'Complete')
	#fig, ax = plt.subplots(figsize=(15,10), ncols=4, nrows=1)
	for name in names:
		fig = plt.figure(figsize=(15,10)) #figure
		gs = gridspec.GridSpec(1, 2, width_ratios = [0.1, 1]) #needed for different sizes
		left   =  0.125  # the left side of the subplots of the figure
		right  =  0.9    # the right side of the subplots of the figure
		bottom =  0.1    # the bottom of the subplots of the figure
		top    =  0.9    # the top of the subplots of the figure
		wspace =  .5     # the amount of width reserved for blank space between subplots
		hspace =  1.1    # the amount of height reserved for white space between subplots
		# This function actually adjusts the sub plots using the above paramters
		plt.subplots_adjust(
			left    =  left, 
			bottom  =  bottom, 
			right   =  right, 
			top     =  top, 
			wspace  =  wspace, 
			hspace  =  hspace
			)
		# The amount of space above titles
		y_title_margin = 1.2
		#plt.suptitle("Simulation 11", y = 1.09, fontsize=20)
		ax = [plt.subplot(gs[i]) for i in range(2)]
		ax[0].set_title("Gene type", fontsize = 18)
		ax[1].set_title(name + " regulation", fontsize = 20)
		sns.heatmap(reg[['color']], linewidths=.1, cbar = False, cmap=paleta, yticklabels = 1, ax=ax[0])
		sns.heatmap(reg[[i for i in reg.columns if name in i]], linewidths=.1, cbar = True, cmap="RdBu_r",center = 0, yticklabels = 1, ax = ax[1])
		# Set all labels on the row axis of subplots for square_feet data to "square_feet"
		[ax[i].set_xticklabels(tf_names, fontsize=15, rotation = 0) for i in range(2)]
		[ax[i].set_yticklabels([str(i) for i in wt.axes[0]],rotation = 0, fontsize=10) for i in range(2)]
		ax[0].set_ylabel('gene', fontsize = 15)
		ax[0].set_xticklabels([''])
		open_pdf.savefig()
	open_pdf.close()
def draw_network(pckname, filetype=['finalExpression', 'annotation', 'TFset']):
	import networkx as nx
	import copy
	plt.clf()
	wt = readData(pckname, filetype[0])
	d = readData(pckname, filetype[1])
	tfset = readData(pckname, filetype[2])
	d['direction'] = [-1 if wt.type[i] == '-1' else 1  for i in d.tf_ind]
	wt['color'] = [np.where(i == wt.type.unique())[0][0] for i in wt.type] 
	G=nx.MultiDiGraph()
	for i in wt.axes[0]:
		G.add_node(i, gene_type = wt.color[i])
	G.add_edges_from([(a, b, {'weight':c, 'color':e}) for a, b, c, e in zip(d.tf_ind, d.gene, d.percent, d.direction)])
	edge_colors=['red' if i == -1 else 'blue' for i in d.direction]
	nx.draw(G, with_labels=True, node_color = wt.color*0.5, edge_color = edge_colors, cmap = 'Paired', vmin = 0, vmax = len(wt.type.unique())+1, label = "Complete network")
	open_pdf = PdfPages(pckname + '_basicNetworkX.pdf')
	open_pdf.savefig()
	plt.clf()
	#plt.show()
	cellexp = [i for i in wt.columns if 'exp' in i]
	c =0
	for cell in cellexp:
		G2 = copy.deepcopy(G)
		e = [i for i in wt.axes[0] if wt[cell][i] == 0]
		G2.remove_nodes_from(e)
		edge_colors=['red' if i == -1 else 'blue' for i in d.direction]
		edge_colors = [i for i, j, k in zip(edge_colors, d.tf_ind, d.gene) if wt[cell][j]>0 and wt[cell][k]>0 ]
		nodecolors = [i for i, j in zip(wt.color, wt[cell]) if j>0 ]
		nx.draw(G2, with_labels=True, node_color = nodecolors, edge_color = edge_colors, cmap = 'Paired', vmin = 0, vmax = len(wt.type.unique())+1, label = "Network active in cell " + str(c))
		open_pdf.savefig()
		plt.clf()
		c+=1
	open_pdf.close()
		

def main():
	graphics = sys.argv[1]
	colormap = plt.cm.gist_ncar
	f = os.listdir()
	datafile = [i[:-4] for i in f if i.endswith('.pck')]
	for pckname in datafile:
		if('e' in graphics):
			finalExpressionHeatmap(pckname)
		if('m' in graphics):
			mutantExpressionHeatmap(pckname)
		if('t' in graphics):
			timeExpression(pckname, colormap)
		if('g' in graphics):
			motifGenomeDiagram(pckname)
		if('r' in graphics):
			regulationHeatmap(pckname)
		if('n' in graphics):
			draw_network(pckname)



if __name__ == "__main__":
    main()

