import matplotlib.pyplot as plt
import numpy as np
import cellEvolver as ce
import pandas as pd
import imp
import sys


### To use this script, you have three possibilities:
### 1) Starting from a .pck file (containing a multiCellEvolver object with chromatin features) that was produced during a simulation run in console mode in python. 
###      You can load this .pck file in console mode and call the functions of this script.
### 2) If the simulation was not run from the python console, i.e., cellEvolver.py was called from unix console, you should instead call cellEvolver.py again with arguments as in this example:
###      python cellEvolver.py plotXsss1 -./180530_chromtrials/180601_4cell_mce1inh5Xsss1g1379.pck eval "import simple_plot_chromatin_from_mce as sp;sp.plotChromatinOpening(mce)"

### Example:
### python cellEvolver.py newname -chromatin_last_pck/4cell_mce0Xss/180601_4cell_mce0Xsss16g2217.pck eval 'import simple_plot_chromatin_from_mce as ckp;ckp.plotChromatinOpening(mce)'



def plotChromatinOpening(mce, time_to_print='max'):
	#imp.reload(ce)
	#fname='try_chromatin1g1191.pck'
	#fname  = './180530_chromtrials/180530_4cell_mce0Xb_try1g1083.pck'
	#mce=ce.multicellTournament.readMCE(fname)
	name = mce.instance_name
	e = np.array([o.error for o in mce.population])
	o = mce.population[np.where(e == np.min(e))[0][0]]
	o.chromatinData(mce.gxcalc)
	xx = o.expression
	writeChromatinTables(xx, name, t=time_to_print)
	ymin = -0.1
	ymax = mce.gxcalc.max_opening + 0.1
	for seqnum in range(mce.template.genome_size):
		ann = o.ann[seqnum]
		pos = ann.sites.getPos()
		contador = 1
		for i in range(mce.gxcalc.tmax):
			for cell in range(len(xx)):
				x = xx[cell].chromatinState
				concentrations = xx[cell].concentrations/np.max(xx[cell].concentrations)
				plt.subplot(mce.gxcalc.tmax, mce.template.ncells, contador)
				if(i == 0):
					plt.title('cell ' + str(cell))
				plt.plot(range(x[i][seqnum].shape[0]),x[i][seqnum])
				if(i !=  mce.gxcalc.tmax-1):
					plt.xticks([])
				if(cell==0):
					plt.ylabel('t = ' + str(i), fontsize = 6)
				else:
					plt.ylabel('')
					plt.yticks([])
				plt.ylim((ymin, ymax)) 
				plt.xlim((0, x[i][seqnum].shape[0]))
				
				for start, end, tf in zip(pos/x[i][seqnum].shape[0], (pos + 3)/x[i][seqnum].shape[0], ann.sites.getTFs()):
					if(concentrations[tf, i] < 0.01):
						cc = 'gray'
						alpha = 0.4
					elif(mce.template.tfs.direction[tf] == mce.template.tfs.INHIBITOR_TYPE):		
						cc = '#d62728'
						alpha = 0.2 + 0.8*concentrations[tf, i]
					else:
						cc=  '#2ca02c' ##'#1f77b4'
						alpha = 0.2 + 0.8*concentrations[tf, i]
					y = (tf/len(mce.template.tfs.inds))/(ymax+0.2)
					plt.axhline(xmin=start,xmax=end, y=y, linewidth=2.8, color=cc, alpha=alpha)
				contador+=1
		plt.savefig(name + 'g' + str(seqnum) + '.pdf')
		plt.close('all')
	#now plot all genes, last instant
	#contador = 1 use this to modify code and plot all genes in a single plot
	for seqnum in range(mce.template.genome_size):
		ann = o.ann[seqnum]
		pos = ann.sites.getPos()
		i = mce.gxcalc.tmax -1
		for cell in range(len(xx)):
			x = xx[cell].chromatinState
			concentrations = xx[cell].concentrations/np.max(xx[cell].concentrations)
			#plt.subplot(mce.template.genome_size, mce.template.ncells, contador)
			plt.subplot(mce.template.ncells, 1, cell+1)
			if(i == 0):
				plt.title('cell ' + str(cell))
			plt.plot(range(x[i][seqnum].shape[0]),x[i][seqnum])
			if(cell !=  len(xx)-1):
				plt.xticks([])
			plt.ylabel('cell ' + str(cell), fontsize = 16)
			plt.ylim((ymin, ymax)) 
			plt.xlim((0, x[i][seqnum].shape[0]))
			
			for start, end, tf in zip(pos/x[i][seqnum].shape[0], (pos + 3)/x[i][seqnum].shape[0], ann.sites.getTFs()):
				if(concentrations[tf, i] < 0.01):
					cc = 'gray'
					alpha = 0.4
				elif(mce.template.tfs.direction[tf] == mce.template.tfs.INHIBITOR_TYPE):		
					cc = '#d62728'
					alpha = 0.2 + 0.8*concentrations[tf, i]
				else:
					cc=  '#2ca02c' ##'#1f77b4'
					alpha = 0.2 + 0.8*concentrations[tf, i]
				y = (tf/len(mce.template.tfs.inds))/(ymax+0.2)
				plt.axhline(xmin=start,xmax=end, y=y, linewidth=2.8, color=cc, alpha=alpha)
			#contador+=1
		plt.savefig(name + '_finalt_' + 'gene' + str(seqnum) + '.pdf')
		plt.close('all')


def writeChromatinTables(chromdata,name, t='max'):
	if(t == 'all'):
		mat= np.zeros([chromdata[0].chromatinState[0][0].shape[0]*len(chromdata[0].chromatinState)*len(chromdata[0].chromatinState[0]), len(chromdata)])
		cond = [(g, t, nt) for g in range(len(chromdata[0].chromatinState[0])) for t in range(len(chromdata[0].chromatinState)) for nt in range(chromdata[0].chromatinState[0][0].shape[0])]
		all_time = range(len(chromdata[0].chromatinState))
	elif(t == 'max'):
		mat= np.zeros([chromdata[0].chromatinState[0][0].shape[0]*1*len(chromdata[0].chromatinState[0]), len(chromdata)])
		cond = [(g, len(chromdata[0].chromatinState)-1, nt) for g in range(len(chromdata[0].chromatinState[0])) for nt in range(chromdata[0].chromatinState[0][0].shape[0])]
		all_time = [len(chromdata[0].chromatinState) - 1]
	cond = np.array(cond)
	for g in range(len(chromdata[0].chromatinState[0])):
		for ti in all_time:
			for cell in range(len(chromdata)):
				mat[np.array([a and b for a, b in zip(cond[:, 0]== g, cond[:, 1] == ti)]) , cell] = chromdata[cell].chromatinState[ti][g]
	df1 = pd.DataFrame(cond, columns = ['gene', 'dev_time', 'nt_position'])
	df2 = pd.DataFrame(mat, columns = ['cell' + str(i) for i in range(mat.shape[1])])
	df3 = pd.concat([df1, df2], axis=1)
	df3.to_csv(name + '_chromatinState_time' + t + '.csv')


	


