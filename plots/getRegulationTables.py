import os
import sys
import pandas as pd
import numpy as np
#import matplotlib.pyplot as plt
#from matplotlib.backends.backend_pdf import PdfPages
#from matplotlib import gridspec

def readData(pckname, filetype):
	f = os.listdir()
	datafile = [i for i in f if i.endswith('.csv') and pckname in i and filetype in i]
	assert len(datafile) == 1, "More than one possible file: " + pckname + ' - ' + filetype
	datafile = datafile[0]
	data = pd.read_csv(datafile, sep = "\t", header=0, index_col=0)
	return data

def regulationTable(pckname, filetype = ['finalExpression', 'annotation', 'TFset']):
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


def main():
        os.chdir("C:/Users/Carlos/Desktop/cellEvolver_all/1802_cellevolver/all_finaldata")
        f = os.listdir()
        for f2 in f:
                os.chdir(f2)
                f3 = os.listdir()
                f3 = [i for i in f3 if i.endswith(".csv")]                
                datafile = [i[0:i.rfind('_')+1] for i in f3]
                for pckname in datafile:
                        regulationTable(pckname)
                os.chdir("../")



if __name__ == "__main__":
    main()

