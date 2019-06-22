import numpy as np
import scipy.special 
from scipy.special import comb
from collections import namedtuple
import itertools
import copy as cp
import sys
import time
import re
import pandas as pd

motifInstance = namedtuple("motifInstance", "motif positionHash intensities")
subnetwork = namedtuple("subnetwork", "mat reducedmat originalmat rownames colnames rowsum colsum tipo simid cell")
subsubnetwork = namedtuple("subsubnetwork", "mat colnames rowsum colsum tipo simid cell sub_colnames")
type1define = {'lineage_this_cell':3, 'lineage_other_cell':4, 'other_activator':8, 'inhibitor':9, 'inhibitor_lineage_this_cell':10, 'inhibitor_lineage_other_cell':11,
		'lineage_many_this':5, 'lineage_many_other':6,'lineage_all':7, 'terminal_specific_this':12, 'terminal_specific_other':15, 'terminal_2_this':13, 'terminal_2_other':16, 'terminal_all':14, 'tf':1, 'non_tf':2, 'lin':0}
restrict = {0:[type1define['lin'], type1define['non_tf']], 
		    1:[type1define['lineage_this_cell'], type1define['terminal_specific_this'], type1define['terminal_all'], type1define['terminal_2_this']],
		    2:[type1define['lineage_this_cell'], type1define['terminal_specific_this'], type1define['terminal_all'], type1define['terminal_2_this']],
		    3:[type1define['tf'], type1define['non_tf']]}

#tuples that contain the functional requirements o be considered for a motif: (type_number, type, inputs_required, outputs_requireds)
functional_definition = [(0, type1define['lin'], 0, 1), (0, type1define['tf'], 1, 1), (0, type1define['non_tf'], 1, 0)]


class Motif(object):
	restrict = restrict
	def __init__(self, subnet = None, type_to_separate=1, simset='', ID=''):
		self.original_instance = subnet.mat
		self.type_to_separate = type_to_separate
		self.types = subnet.tipo[type_to_separate][subnet.colnames]
		self.colsum = subnet.colsum
		self.matched_instances = [subnet]
		self.matches = {simset:1}
		self.comb_index = None
		self.isomorphs= None
		self.permanentId = ID
		self.createIsomorphicInstances()
	#To do some day: check if isomorphic graphs are also symmetric (it is not a problem since all instances will match with the first isomorphic graph of the symmetry group)
	def createIsomorphicInstances(self):
		steps = [(i, j) for i, j in zip(self.types, self.colsum)]
		steps_set = set(steps)
		steps_count = np.array([steps.count(i) for i in steps_set])
		isomorphs = self.original_instance.reshape((self.original_instance.shape[0], self.original_instance.shape[1],1))
		combnames = np.arange(isomorphs.shape[1]).reshape((1, isomorphs.shape[1]))
		if(np.all(steps_count == 1)):
			self.indexcomb = combnames
			self.isomorphs = isomorphs
			self.comb_index	= combnames
		else:
			permuts = [scipy.special.perm(i, i) for i in steps_count]
			z = np.prod(permuts)
			combinations = []
			#isomorphs = np.zeros((original_instance.shape[0], original_instance.shape[1], z))
			isomorphs = isomorphs.repeat(z, axis = 2)
			combnames = combnames.repeat(z, axis = 0)
			for i, j in enumerate(steps_set):
				if(steps_count[i] > 1):
					ind = np.where([k==j[0] and l==j[1] for k, l in zip(self.types, self.colsum)])[0].tolist()
					combinations.append(list(itertools.permutations(ind)))
			aux = combinations[0]
			for cc in range(1, len(combinations)):
				aux = list(itertools.product(aux, combinations[cc]))
				aux = [i+j for i, j in aux]
			indexcomb = np.array(aux)
			assert indexcomb.shape[0] == isomorphs.shape[2], "Isomorphs: index combinations " + str(indexcomb.shape[0]) + ", isomorph matrix " + str(isomorphs.shape[2])
			for cc in range(1, indexcomb.shape[0]):
				isomorphs[:,indexcomb[0, :],cc] = isomorphs[:,indexcomb[cc, :],cc]
				isomorphs[indexcomb[0, :], :,cc] = isomorphs[indexcomb[cc, :], :,cc]
				combnames[cc, indexcomb[0, :]] = combnames[cc, indexcomb[cc, :]]
			self.isomorphs = isomorphs	
			self.comb_index	= combnames
	def compare(self, mat, simset='', memorizeSubnet = True):
		if(not np.all(self.colsum == mat.colsum) or not np.all(self.types == mat.tipo[self.type_to_separate][mat.colnames])):
			return False
		else:
			found =  False
			i = 0
			while(not found and i < self.isomorphs.shape[2]):
				found = np.all(mat.mat == self.isomorphs[:,:,i])	
				i+=1
			if(not found):
				return False
			else:
				self.addMatch(mat, simset, memorizeSubnet = True)
				return True
	def addMatch(self, mat, simset='',memorizeSubnet = True):
		if(memorizeSubnet):
			self.matched_instances.append(mat)
		nmatches = self.matches.get(simset, 0) + 1
		self.matches[simset] = nmatches
		
	def setPermanentId(self, pid=""):
		self.permanentId = pid
	def toString(self):
		return 'mat:' + re.sub('[\s+]','',str(self.isomorphs[:,:,0])) + '|types:'+str(self.types) +'|' + str(self.type_to_separate)
	def getNormalizedPositionNamesNaive(self):
		return [str(self.permanentId) +'_t' + str(self.type_to_separate)+ '.' + str(int(self.types[i])) + '_p' + str(i) for i in range(self.types.shape[0])]
	def getNormalizedPositionNames(self):
		naive_id = self.getNormalizedPositionNamesNaive()
		if(self.isomorphs.shape[2]==1):
			return naive_id
		else:
			list_of_equal_sets = []
			#equals = np.zeros((self.isomorphs.shape[2],self.isomorphs.shape[2]))
			for i in range(self.isomorphs.shape[2]):
				for j in range(i):
					if(np.all(self.isomorphs[:,:,i] == self.isomorphs[:,:,j])):
						#equals[i, j] = 1
						dif_in_order = np.where(self.comb_index[i, :] != self.comb_index[j,:])[0]	#this was a bug
						dif_in_order_names = self.comb_index[i, np.where(self.comb_index[i, :] != self.comb_index[j,:])[0] ]
						#equals2 = np.zeros((dif_in_order.shape[0], dif_in_order.shape[0]))
						for k in range(len(dif_in_order)):
							for l in range(k):
								ind_k_in_i = np.where(self.comb_index[i,:] == dif_in_order_names[k])[0]
								ind_l_in_j = np.where(self.comb_index[j,:] == dif_in_order_names[l])[0]
								if(np.all(self.isomorphs[:,ind_k_in_i,i] == self.isomorphs[:,ind_l_in_j,j]) and np.all(self.isomorphs[ind_k_in_i,:,i] == self.isomorphs[ind_l_in_j,:,j])):
									found_set = False
									for s in list_of_equal_sets:
										if(dif_in_order_names[k] in s):
											found_set = True
											s.add(dif_in_order_names[l])
											break
										elif(dif_in_order_names[l] in s):
											found_set = True
											s.add(dif_in_order_names[k])
											break
									if(not found_set):
										list_of_equal_sets.append({dif_in_order_names[k], dif_in_order_names[l]})
			for i in range(len(list_of_equal_sets)):
				s = list_of_equal_sets[i]
				setname = str(self.permanentId) +'_t' + str(self.type_to_separate)+ '.' + str(int(self.types[next(iter(s))])) + '_set' + str(i)
				for j in s:
					naive_id[j] = setname
		return naive_id

	##returns a hash: {simulation_id:{gene_id:[(cell, normalized_position_in_motif),...]}}
	def getNormalizedPositionsOfInstances(self, removeRandom = True):
		normpos = {}
		posnames = self.getNormalizedPositionNames()
		for match in self.matched_instances:
			if (removeRandom and 'random' in match.simid ):
				continue
			simdict = normpos.get(match.simid, {})
			isorder = self.__getIsomorphOrder(match)
			for g in range(match.colnames.shape[0]):
				simdict.setdefault(match.colnames[g], []).append((match.cell, posnames[isorder[g]]))
			normpos[match.simid] = simdict
		return normpos
			
	def __getIsomorphOrder(self, subnet):
		res = None
		for i in range(self.isomorphs.shape[2]):
			if(np.all(subnet.mat == self.isomorphs[:,:,i])):
				res = self.comb_index[i,:]
				break
		if(res is None):
			raise ValueError("One of the instances did not match the isomorphs!")
		return res
	def resetInstances(self):
		#self.matched_instances = [self.matched_instances[0]]
		self.matched_instances.clear()

class motifContainer(object):
	threshold_bin = 0
	MAX_MOT_SIZE = 10
	def __init__(self, simulSet, motifs = None, type_to_separate = 1, type_to_separate_motifs=3,compress_terminal_features = True, addNewMotifs=True, preprocessedNets=None, motifLastId = None, ignoreSelfReg = False):
		networkInstances = simulSet.simulations
		self.path = simulSet.path
		self.set_name = self.path.split("/")[-1]
		if(self.set_name == ''):
			self.set_name = self.path.split("/")[-2]
		self.type_to_separate = type_to_separate	#only used to join terminal features by group
		self.type_to_separate_motifs = type_to_separate_motifs	#used to constrain nodes that can be occupied by a type in motifs
		self.ignoreSelfReg = ignoreSelfReg
		if(preprocessedNets is None):
			self.networks = [self.preprocessNetwork(n, compress_terminal_features) for n in networkInstances]
		else:
			self.networks = preprocessedNets
		if(motifs is None):
			self.motifs = {}
			for i in range(min(networkInstances[0].genome_size, self.MAX_MOT_SIZE)):
				self.motifs[i] = []
		else:
			self.motifs = motifs
		self.addNewMotifs = addNewMotifs
		self.summaryTables = None
		self.functional_definition = functional_definition
	def preprocessNetwork(self, network, compress_terminal_features = True):
		cellnets = []
		#file = open('testfile.txt','w') 
		#file.write("NETWORK ")
		for cell in range(network.ncells):

			#file.write('\n**********************************************************\ncell ' + str(cell) + '\n**********************************************************\n')
			tipo = self.typedefine(cell, network)
			this_cell = np.array(network.regtable[[i for i in network.regtable.columns if (('Cell ' + str(cell)) in i)]])
			gene = np.array(network.regtable.gene)
			#Add rows for genes that may be lacking(due to TF removal in previous steps)
			lacking_genes = np.array([i for i in range(np.max(gene)) if i not in gene])
			gene = np.concatenate((gene, lacking_genes), axis = 0)
			this_cell = np.concatenate((this_cell, np.zeros((lacking_genes.shape[0], this_cell.shape[1]))), axis = 0)
			sort_gene = np.argsort(gene)
			gene = gene[sort_gene].astype(int)
			this_cell = this_cell[sort_gene, :]
			#file.write('tipo: ' + str(tipo) + '\nthis_cell:\n ' + str(this_cell) + '\ngene: ' + str(gene) + '\ncompressing:')
			#average terminal features
			if(compress_terminal_features):
				ind = [i for i, g in enumerate(gene) if  type1define['non_tf'] != tipo[0][g]] #reduced matrix only with TFs in rows
				redmat = this_cell[ind, :]
				redgene = gene[ind]
				i_ter = ind[-1] + 1
				for terminal_type in list(set(tipo[self.type_to_separate][[tt for tt in range(i_ter, len(tipo[0])) if tipo[0][tt] == type1define['non_tf']]])):
					#print(terminal_type)
					#file.write('\t\t terminal_type: ' + str(terminal_type))
					aux = np.mean(this_cell[[i for i in range(this_cell.shape[0]) if  terminal_type == tipo[self.type_to_separate][gene[i]]] ,:], axis = 0).reshape((1, redmat.shape[1]))
					redmat = np.concatenate((redmat, aux), axis = 0)
					redgene = np.concatenate((redgene, gene[np.where(tipo[self.type_to_separate][gene] == terminal_type)][0].reshape(1)), axis = 0)
				#file.write('ind: ' +str( ind) + '\nredmat: \n' + str(redmat) + '\nredgene: ' + str(redgene)+ '\ncompressing:')
			else:
				## NOT IMPLEMENTED (just fill the redgene matrix)
				raise Exception("Sorry! Not implemented with compress_terminal_features  = False")			
			#binarize
			mat = np.where(np.abs(redmat)>self.threshold_bin, 1, 0)
			#file.write('binarized mat: \n' + str( mat))
			#Add columns for terminal features
			rowgene = redgene
			colgene = np.concatenate((np.array([i for i in range(mat.shape[1])]), rowgene[tipo[0][rowgene] == type1define['non_tf']]))  #in columns, index == number of gene
			mat = np.concatenate((mat, np.zeros((mat.shape[0], np.where(rowgene[tipo[0][rowgene] == type1define['non_tf']])[0].shape[0]))), axis=1)
			redmat = np.concatenate((redmat, np.zeros((redmat.shape[0], np.where(rowgene[tipo[0][rowgene] == type1define['non_tf']])[0].shape[0]))), axis=1)
			#file.write('Added terminal feature rows:\n\trowgene: ' + str(rowgene) + '\n\tcolgene: ' + str(colgene) +'\n\t matrix:\n ' +str(mat) + '\n\tredmat:\n ' + str(redmat))
			#remove genes that are doing nothing (0s in rows and columns)
			coltoretain = np.sort(np.where(np.sum(mat, axis = 0) > 0)[0])
			rowtoretain = np.sort(np.where(np.sum(mat, axis = 1) > 0)[0])
			toretain = np.sort(np.union1d(coltoretain, rowtoretain))
			redmat = redmat[toretain, :][:, toretain]
			mat = mat[toretain, :][:, toretain]
			rowgene = redgene[toretain]
			colgene = colgene[toretain]	
			#file.write('Zeros removed:\n\tcoltoretain: ' + str(coltoretain)+ '\n\trowtoretain: '+str( rowtoretain) +'\n\ttoretain: ' + str(toretain) + '\n\tredmat:\n '+ str(redmat) + '\n\tmat: ' + str(mat) + '\n\trowgene: ' + str(rowgene) +'\n\tcolgene: ' +  str(colgene))
			assert mat.shape[0] == mat.shape[1], 'matrix with different dimension sizes'
			assert np.all(rowgene == colgene), 'different order in matrix rows and columns'
			# Sort by type
			sortorder = np.argsort(tipo[self.type_to_separate_motifs][colgene])
			colgene = colgene[sortorder]
			rowgene = rowgene[sortorder]
			mat = mat[sortorder,:]
			mat = mat[:, sortorder]
			redmat = redmat[sortorder,:]
			redmat = redmat[:, sortorder]
			# Remove self regulation if necessary:
			if(self.ignoreSelfReg):
				for i in range(mat.shape[0]):
					mat[i,i] = 0
					#redmat[i,i] = 0
			#finally, order by number of output connections to break symmetry
			colsum = np.sum(mat, axis=0)
			#sort cols and rows by col
			for tsort in np.sort(np.unique(tipo[self.type_to_separate_motifs][colgene])):
				ind = np.where(tipo[self.type_to_separate_motifs][colgene] == tsort)[0]
				if(ind.shape[0] > 1):
					indsorted = ind[np.argsort(colsum[ind])]
					mat[:, ind] = mat[:, indsorted]
					mat[ind, :] = mat[indsorted, :]
					redmat[:, ind] = redmat[:, indsorted]
					redmat[ind, :] = redmat[indsorted, :]
					colsum[ind] = colsum[indsorted]
					colgene[ind] = colgene[indsorted]
			rowsum = np.sum(mat, axis=1)
			rowgene = colgene
			#file.write('ORDERED:\n\tcolsum: ' + str(colsum) + '\n\trowsum: ' + str(rowsum) +'\n\tredmat: ' + str(redmat) + '\n\tmat: ' + str(mat) + '\n\trowgene: ' + str(rowgene) +'\n\tcolgene: ' + str(colgene))
			#add this cell network
			cellnets.append(subnetwork(mat, redmat, this_cell, rowgene, colgene, rowsum, colsum, tipo, network.name, cell))
		#file.close()
		return cellnets
			
	def typedefine(self, cell, network):
		t = np.array(network.type3)
		t0 = np.zeros(len(network.type3))
		t1 = np.zeros(len(network.type3))
		t2 = np.zeros(len(network.type3))
		t3 = np.zeros(len(network.type3))
		t4 = np.zeros(len(network.type3)) #will not be used
		t5 = np.zeros(len(network.type3))
		init = np.array(network.getInitialExpr(cell))
		opt = np.array(network.getOptimalExpr())
		dict_of_combos = {}
		first_new_value = max(type1define.values()) + 1
		for i in range(t.shape[0]):
			if('lin' in t[i] and '-1' not in t[i]):
				t3[i] = type1define['tf']
				if('[1]' in t[i] and init[i] > 0):
					t0[i] = type1define['lin']
					t1[i] = type1define['lineage_this_cell']
					t2[i] = type1define['lineage_this_cell']
				elif('[1]' in t[i] and init[i] == 0):
					t0[i] = type1define['tf']
					t1[i] = type1define['other_activator']
					t2[i] = type1define['lineage_other_cell']
				elif('[' + str(network.ncells) + ']' in t[i]):
					t0[i] = type1define['lin']
					t1[i] = type1define['lineage_all']
					t2[i] = type1define['lineage_all']
				elif('[0]' in t[i]):
					t0[i] = type1define['tf']
					t1[i] = type1define['other_activator']
					t2[i] = type1define['other_activator']
				else:
					if(init[i] > 0):
						t0[i] = type1define['lin']
						t1[i] = type1define['lineage_many_this']
						t2[i] = type1define['lineage_many_this']
					else:
						t0[i] = type1define['tf']
						t1[i] = type1define['other_activator']
						t2[i] = type1define['lineage_many_other']
			elif('lin' in t[i] and '-1' in t[i]):
				t3[i] = type1define['tf']
				t0[i] = type1define['tf']
				if(init[i] > 0):
					t1[i] = type1define['inhibitor_lineage_this_cell']
					t2[i] = type1define['inhibitor_lineage_this_cell']
				elif(init[i] == 0):
					t1[i] = type1define['inhibitor']
					t2[i] = type1define['inhibitor']
			else:
				t0[i] = type1define['non_tf']
				t3[i] = type1define['non_tf']
				l = len(np.where(opt[i, :] > 0)[0])
				if(l == 1 and opt[i, cell] > 0):
					t1[i] = type1define['terminal_specific_this']
					t2[i] = type1define['terminal_specific_this']
				elif(l == 1 and opt[i, cell] == 0):		
					t1[i] = type1define['terminal_specific_other']
					t2[i] = type1define['terminal_specific_other']
				elif(l > 1 and l < network.ncells and opt[i, cell] > 0):
					t1[i] = type1define['terminal_2_this']
					t2[i] = type1define['terminal_2_this']
					dict_of_combos[i] = ''.join([str(j) for j in np.where(opt[i,:] > 0)[0]])
				elif(l > 1 and  l < network.ncells and opt[i, cell] == 0):		
					t1[i] = type1define['terminal_2_other']
					t2[i] = type1define['terminal_2_other']
				elif(l == network.ncells):
					t1[i] = type1define['terminal_all']
					t2[i] = type1define['terminal_all']
		set_of_combos = set(dict_of_combos.values())
		t5 = cp.deepcopy(t1)
		if(len(set_of_combos) > 1):
			for c, cnum in zip(set_of_combos, range(len(set_of_combos))):
				for gnum, cc in zip(dict_of_combos.keys(), dict_of_combos.values()):
					if(c == cc):
						t5[gnum] = first_new_value + cnum
		return (t0, t1, t2, t3, t4, t5)
	def findMotifs(self, starting_n = 0, selectBy=0, functionalFilter = True):
		start_time = time.time()
		sim_cont= 0
		for sim in self.networks:
			cell_cont=0
			for subnet in sim:
				subnet_generator = subNetworkGenerator(subnet, selectBy=selectBy, type_to_separate=self.type_to_separate_motifs, starting_n=starting_n)
				try:
					while(not subnet_generator.complete): 
						sub = subnet_generator.next()
						if(functionalFilter):
							if(not self.functionalSubnetwork(sub)):
								continue
						if(sub is None):
							break
						n = sub.mat.shape[0]
						sub_present = False
						for mot in self.motifs[n]:
							sub_present = mot.compare(sub, simset=self.set_name, memorizeSubnet = self.addNewMotifs)
							if(sub_present):
								break
						if(not sub_present and self.addNewMotifs):
							newMotif = Motif(sub, self.type_to_separate_motifs, simset=self.set_name, ID = 's' + str(n) + 'mot' + str(len(self.motifs[n])))
							self.motifs[n].append(newMotif)
				except:
					print("Error while searching for motifs: sim ", sim_cont, " cell ", cell_cont, " ", self.set_name,"\n")
					print(sys.exc_info())
					return subnet_generator
					
				sim_cont += 1
			cell_cont+=1
		return time.time() - start_time
	def functionalSubnetwork(self, sub, definition = 'def0'):
		isfunctional = True
		g = 0
		while(isfunctional and g<sub.colnames.shape[0]):
			for f in self.functional_definition:
				if(sub.tipo[f[0]][sub.colnames[g]] != f[1]):
					continue
				else:
					if(f[2] > sub.rowsum[g] - sub.mat[g,g] or f[3] > sub.colsum[g] - sub.mat[g,g]):
						isfunctional=False
						break
			g+=1
		return isfunctional
	def __iter__(self):
		for i in self.motifs.keys():
			for j in range(len(self.motifs[i])):
				yield self.motifs[i][j]
	def printMotifSummary(self):
		print("Number of motifs: ", [len(self.motifs[i])for i in self.motifs.keys()])
		for i in self.motifs.keys():
			print("********* ", str(i), "  *********")
			for j in range(len(self.motifs[i])):
				print("___ ", str(j), " _____")
				if(len(self.motifs[i][j].matched_instances) == 0):
					print(self.motifs[i][j].original_instance, "\ntypes: ", self.motifs[i][j].types, '\n')
				else:
					print(self.motifs[i][j].original_instance, "\ntypes: ", self.motifs[i][j].types, "\nnames: ", self.motifs[i][j].matched_instances[0].colnames, "\n", self.motifs[i][j].matched_instances[0].simid, " ", self.motifs[i][j].matched_instances[0].cell, "\n" )
	def writeMotifSummary(self, fname = None, outpath = ''):
		if(fname is None):
			fname = outpath + self.set_name + '_motifSummary.txt'
		else:
			fname = outpath + fname
		file = open(fname, 'w')
		file.write("Type definitions:\n" + str(type1define))
		file.write("\n\nNumber of motifs: " + str([len(self.motifs[i])for i in self.motifs.keys()]))
		for i in self.motifs.keys():
			file.write("\n\n********* " + str(i) + "  *********\n")
			for j in range(len(self.motifs[i])):
				file.write("___ " + str(j) + " _____\n")
				file.write(str(self.motifs[i][j].original_instance) + "\ntypes: " + str(self.motifs[i][j].types) + '\n')
				if(len(self.motifs[i][j].matched_instances) > 0):
					file.write("\nnames: " + str(self.motifs[i][j].matched_instances[0].colnames) +
						"\n" + 'simid: ' + str(self.motifs[i][j].matched_instances[0].simid) + 
						" cell: " + str(self.motifs[i][j].matched_instances[0].cell) + 
						'\nn_matches: ' + str(self.motifs[i][j].matches) +"\n")
		file.close()
	def makeTables(self, addToName='', min_matches = 40, write=True):
		start_time=time.time()
		try:
			self.summaryTables = {}
			self.__makeMotifTable()
			self.__makeTFTable(min_matches)
			self.__JoinGenesByCell()
			self.__makeCellMotifTable()
			self.__makeType2ByPosition()
			if(write):
				self.__writeOutputTables(addToName, min_matches)
		except ValueError as e:
			print('Value error - Error making tables:\n', e.args)
		except:
			print('Some kind of error making tables')
		return time.time() - start_time
	def __setPermanentIds(self):
		for motsize in self.motifs.keys():
			mot_id = 0
			for mot in self.motifs[motsize]:
				mot.setPermanentId(str(motsize) + '_mot' + str(mot_id))
				mot_id+=1
	def __makeMotifTable(self):
		i = 0
		df = pd.DataFrame(columns=['motif_id', 'motif_str', 'motif_size', 'functional'])
		for j in self.motifs.keys():
			for k in range(len(self.motifs[j])):
				sims  = np.array(list(self.motifs[j][k].matches.keys()))
				newsims = sims[[i for i in range(sims.shape[0]) if not sims[i] in df.columns.tolist()]]
				oldsims = sims[[i for i in range(sims.shape[0]) if sims[i] in df.columns.tolist()]]
				newrow = [self.motifs[j][k].permanentId, self.motifs[j][k].toString(), j, str(np.any([self.functionalSubnetwork(mi) for mi in self.motifs[j][k].matched_instances]))] + [self.motifs[j][k].matches[s] if s in oldsims else 0 for s in df.columns[4:]] + [self.motifs[j][k].matches[s] for s in newsims]
				for s in newsims:
					df[s] = [0 for ii in df.index]
				df.loc[i] = newrow
				i+=1
		try:
			self.summaryTables['motif_table'] = df
			self.summaryTables['motifIdList'] = df.motif_id.tolist()
		except:
			print("self.summaryTables not initialized - data frame not stored, call makeTables method\n")
			return df3
	def __makeTFTable(self, min_matches=40): 
		assert 'motif_table' in self.summaryTables.keys(), 'call __makeMotifTable first'
## rows with TF/cell identification
		tuplas = [[re.sub('_$','',net[cell.cell].simid) + '_' + str(genenum),re.sub('s[0-9]+_*' ,'',net[cell.cell].simid),re.sub('_$','',net[cell.cell].simid),cell.cell, genenum] + [tt[genenum] for tt in net[cell.cell].tipo]  for net in self.networks for cell in net for genenum in cell.colnames]
		df  = pd.DataFrame(tuplas, columns =  ['gene_unique_id', 'simset','simid','cell', 'gene'] + ['type' + str(i) for i in range(len(self.networks[0][0].tipo))])
		df = df.set_index(['simid', 'cell', 'gene'], drop=False)
		df['number'] = [i for i in range(len(df.index))]
		#columns with motif id
		#motpos_ids = [posname for motlen in self.motifs.keys() for mot in self.motifs[motlen] for posname in set(mot.getNormalizedPositionNames())] #making a set of NormalizedPositionNames was necessary after making position sets
		motpos_ids = [posname for mot in self for posname in set(mot.getNormalizedPositionNames())] #making a set of NormalizedPositionNames was necessary after making position sets
		mot_ids = {i:j for i, j in zip(motpos_ids, range(len(motpos_ids)))}
		#matrix to store the output
		mat=np.zeros((len(tuplas), len(mot_ids)))
		for j in self.motifs.keys():
			for mot in self.motifs[j]:
				##returns a hash: {simulation_id:{gene_id:[(cell, normalized_position_in_motif),...]}}
				mot_matches = mot.getNormalizedPositionsOfInstances()
				for sid in mot_matches.keys():
					for gid in mot_matches[sid].keys():
						for match in mot_matches[sid][gid]:
							sub = re.sub('_$','',sid)
							mat[df.loc[sub,match[0],gid].number, mot_ids[match[1]]]+=1
		df2 = pd.DataFrame(mat, columns = list(mot_ids.keys()))
		df3 = pd.concat([df.reset_index(drop=True), df2], axis=1)
		try:
			self.summaryTables['tfs_cells_sep'] = df3
			#Make also a table with filtered matches
			if(min_matches > 0):
				motids = self.summaryTables['motif_table']['motif_id'].loc[self.summaryTables['motif_table'][self.set_name]>min_matches]
				preserve_columns = [i for i in df3.columns if i.split('_')[0] in list(motids)]
				self.summaryTables['tfs_cells_sep_thresholded'] = df3[list(df3.columns[:10]) + preserve_columns]
		except KeyError:
			print("self.summaryTables not initialized - data frame not stored, call makeTables method\n")
			return df3
	def __JoinGenesByCell(self,  min_matches=40):
		try:
			df = self.summaryTables['tfs_cells_sep'] 
		except KeyError:
			print("tfs_cells_sep not calculated - data frame not stored, call makeTables method\n")
		colnames = list(df.columns[10:])
		#since cells are going to be joined, put other values in type4
		newmap = {x:x for x in type1define.values()}
		newmap[type1define['lineage_other_cell']] = type1define['lineage_this_cell']
		newmap[type1define['terminal_2_other']] = type1define['terminal_2_this']
		newmap[type1define['terminal_specific_other']] = type1define['terminal_specific_this']
		df['type4'] = [newmap[i] for i in df.type2]
		groupers= ['gene_unique_id','simset', 'simid', 'gene', 'type3', 'type4']
		df2 = df[groupers + colnames].groupby(groupers).sum()
		df2 = df2.reset_index(drop=False)
		#since df is a reference, remember to remove type4, but store it for later
		self.summaryTables['type4'] = self.summaryTables['tfs_cells_sep']['type4']
		self.summaryTables['tfs_cells_sep'] = self.summaryTables['tfs_cells_sep'].drop('type4', axis=1)
		#df2 = df2.set_index(['gene_unique_id'],drop=False)
		self.summaryTables['tfs_cells_joint'] = df2
		if(min_matches > 0):
			motids = self.summaryTables['motif_table']['motif_id'].loc[self.summaryTables['motif_table'][self.set_name]>min_matches]
			preserve_columns = [i for i in df2.columns if i.split('_')[0] in list(motids)]
			self.summaryTables['tfs_cells_joint_thresholded'] = df2[list(df2.columns[:6]) + preserve_columns]
		
	def __makeCellMotifTable(self):
		tuplas = [[re.sub('s[0-9]+_*' ,'',net[cell.cell].simid),re.sub('_$','',net[cell.cell].simid),cell.cell] for net in self.networks for cell in net]
		df  = pd.DataFrame(tuplas, columns =  ['simset','simid','cell'])
		df = df.set_index(['simid', 'cell'], drop=False)
		df['number'] = [i for i in range(len(df.index))]
		mot_ids = {m[1].permanentId:m[0] for m in enumerate(self)}

		mat=np.zeros((len(tuplas), len(mot_ids)))
		for m in self:
			for subnet in m.matched_instances:
				mat[df.loc[re.sub('_$','',subnet.simid),subnet.cell].number, mot_ids[m.permanentId]]+=1
		df2 = pd.DataFrame(mat, columns = list(mot_ids.keys()))
		df3 = pd.concat([df.reset_index(drop=True), df2], axis=1)
		try:
			self.summaryTables['cellMotifTable'] = df3
		except:
			print("self.summaryTables not initialized - data frame not stored, call makeTables method\n")
			return df3
	def __makeType2ByPosition(self):
		try:
			tab = self.summaryTables['tfs_cells_sep'].set_index('type2', drop=False)
		except:
			print('No tfs_cells_sep Table or SummaryTables not initialized - call makeTables method')
			return
		types = np.unique(tab.type2)
		typenames = {b:a for a, b in type1define.items()}
		colnames = tab.columns[10:]
		original_mat = tab[colnames]
		mat = np.zeros([len(types), len(colnames)])
		for i in range(len(types)):
			mat[i, :] = np.sum(np.array(tab.loc[types[i]][colnames]), axis = 0)
		mat_sum = np.sum(mat, axis = 1).reshape((len(types),1))
		mat_sum2 = np.sum(mat, axis = 0).reshape((1,len(colnames)))
		mat_norm = mat/mat_sum
		mat_norm2 = mat/mat_sum2
		df0 = pd.DataFrame([(i, typenames[i]) for i in types], columns = ['type2', 'typename'])
		df1 = pd.DataFrame(mat, columns = colnames)
		df2 = pd.DataFrame(mat_norm, columns = colnames)
		df3 = pd.DataFrame(mat_norm2, columns = colnames)
		df1 = pd.concat([df0, df1], axis=1)
		df2 = pd.concat([df0, df2], axis=1)
		df3 = pd.concat([df0, df3], axis=1)
		self.summaryTables['motifsByType2_count'] = df1
		self.summaryTables['motifsByType2_prop_groups'] = df2
		self.summaryTables['motifsByType2_prop_motifs'] = df3
		
	def __writeOutputTables(self, name='', min_mots=40):
		name = '_'.join([self.set_name, name])
		self.summaryTables['motif_table'].to_csv(name + '_motifSummary.csv')
		self.summaryTables['tfs_cells_sep'].to_csv(name + '_geneByCellNormPosition_sep.csv')
		if(min_mots > 0):
			self.summaryTables['tfs_cells_sep_thresholded'].to_csv(name + '_geneByCellNormPosition_sepThreshold'+ str(min_mots) +'.csv')
			self.summaryTables['tfs_cells_joint_thresholded'].to_csv(name + '_geneNormPosition_sumThreshold'+ str(min_mots) +'.csv')
		self.summaryTables['tfs_cells_joint'].to_csv(name + '_geneNormPosition_sum.csv')
		self.summaryTables['cellMotifTable'].to_csv(name + '_CellMotifTable_summary.csv')
		self.summaryTables['motifsByType2_count'].to_csv(name + '_motifsByType2_count.csv')
		self.summaryTables['motifsByType2_prop_groups'].to_csv(name + '_motifsByType2_prop_groups.csv')
		self.summaryTables['motifsByType2_prop_motifs'].to_csv(name + '_motifsByType2_prop_motifs.csv')
		#self.writeMotifSummary(fname = name+'motifSummary.txt')
	def resetInstances(self):
		for n in self.motifs.keys():
			for m in self.motifs[n]:
				m.resetInstances()
	def saveSelf(self, name = ''):
		import _pickle as pck
		filename = self.set_name + '_'+name+ '_motifContainer.pck'
		with open(filename, 'wb') as output: 
 			pickler = pck.Pickler(output, -1)
 			pickler.dump(self)
	@staticmethod
	def readMC(filenamein):
		with open(filenamein, 'rb') as f:
			print('opening ' + filenamein)
			import _pickle as pck
			mc = pck.load(f)
		return mc
				
#subnetwork definition:
# subnetwork namedtuple("subnetwork", "mat reducedmat originalmat rownames colnames rowsum colsum tipo simid cell")
#subsubnetwork = namedtuple("subsubnetwork", "mat colnames rowsum colsum tipo simid cell subcolnames") #subcolnames will contain indices with respect to the original_mat (eg, 0:7), while colnames correspond to the original colnames (eg., 0, 1, 5, 10).
class subNetworkGenerator(object):
	restrict = restrict
	def __init__(self, net, selectBy = 0,type_to_separate=1, starting_n = 0):
		self.original_net = net
		self.selectBy = selectBy
		self.type_to_separate = type_to_separate
		self.original_mat = net.mat
		self.original_nodes = net.colnames
		self.original_tipo = net.tipo[selectBy][self.original_nodes]
		self.original_tipo_sort = net.tipo[type_to_separate][self.original_nodes]	
		self.necessary_types = np.array([i for i in self.restrict[selectBy] if (i in self.original_tipo)])
		self.current_n = self.original_nodes.shape[0] if starting_n ==0 else starting_n
		self.complete = False if self.current_n > len(self.necessary_types) else True
		if(starting_n ==0):
			self.current_combinations = np.array([1 for i in self.original_nodes]).reshape((1, self.original_nodes.shape[0]))
			self.current_combinations_shape = 1
			self.current_row = 0
		else:
			self.current_combinations=np.zeros(10)
			self.current_combinations_shape = 0
			self.current_row = 0
			self.generateCombinations()
	def next(self):
		suitable_combination = None
		while(not self.complete and suitable_combination is None):
			if(self.current_combinations is not None):
				suitable_combination = self.current_combinations[self.current_row, :]
			self.current_row += 1
			if(self.current_row >=  self.current_combinations_shape):
				if(self.current_n <= len(self.necessary_types) or self.current_n <= 2):
					self.complete = True
				else:
					self.current_n-=1
					self.current_row = 0
					self.generateCombinations()
		if(suitable_combination is not None):
			suitable_combination = np.where(suitable_combination == 1)[0]
			suitable_combination = suitable_combination[np.argsort(self.original_tipo_sort[suitable_combination])]
			newmat = self.original_mat[suitable_combination, :][:, suitable_combination]
			newcolsum = np.sum(newmat, axis=0)
			for tsort in np.sort(np.unique(self.original_tipo_sort[suitable_combination])):
				ind = np.where(self.original_tipo_sort[suitable_combination] == tsort)[0]
				if(ind.shape[0] > 1):
					indsorted = ind[np.argsort(newcolsum[ind])]
					suitable_combination[ind] = suitable_combination[indsorted]
			newmat = self.original_mat[suitable_combination, :][:, suitable_combination]
			newcolsum = np.sum(newmat, axis=0)
			newnodes = self.original_nodes[suitable_combination]		
			return subsubnetwork(newmat, newnodes, np.sum(newmat, axis=1), newcolsum, self.original_net.tipo, self.original_net.simid, self.original_net.cell, suitable_combination)
		else:
			return None	
	def generateCombinations(self):
		combmat = np.zeros((comb(self.original_nodes.shape[0], self.current_n, exact=True), self.original_nodes.shape[0]))
		retain = [True for i in range(combmat.shape[0])]
		it = itertools.combinations([i for i in range(self.original_nodes.shape[0])], self.current_n)
		i = 0
		for combin in it:
			if(not np.all([True if t in self.original_tipo[list(combin)] else False for t in self.necessary_types])):
				retain[i] = False
			else:
				auxmat = self.original_mat[combin, :][:, combin] 
				scol = np.where(np.sum(auxmat, axis = 0) > 0)[0]
				srow = np.where(np.sum(auxmat, axis = 1) > 0)[0]
				both = np.union1d(scol, srow)
				if(not len(combin) == both.shape[0]):
					retain[i]=False
				else:
					combmat[i, combin] = 1
			i+=1
		combmat = combmat[retain, :]
		if(combmat.shape[0] > 0):
			self.current_combinations = combmat
			self.current_combinations_shape = combmat.shape[0]
		else:
			self.current_combinations = None
			self.current_combinations_shape = 0
			if(self.current_n <= len(self.necessary_types) or self.current_n <= 2):
				self.complete=True

					
class randomizedContainer(motifContainer):
	#Takes an object of the class motifContained and a string that is used as an ID of the set
	def __init__(self, container, addName = "random_0", randomizeNodeOutputs = True, randomizeNodeInputs = False, addNewMotifs = False, cpmotifs = False):
		self.path = container.path
		self.set_name = container.set_name + addName
		self.type_to_separate = container.type_to_separate
		self.type_to_separate_motifs = container.type_to_separate_motifs
		self.networks = cp.deepcopy(container.networks)
		self.ignoreSelfReg = container.ignoreSelfReg
		self.randomizeNetworks(randomizeNodeOutputs, randomizeNodeInputs)
		if(cpmotifs):
			self.motifs = cp.deepcopy(container.motifs)
		else:
			self.motifs = container.motifs
		self.addNewMotifs = addNewMotifs
		self.functional_definition = functional_definition
	def randomizeNetworks(self, outputs = True, inputs = False):
		for i in range(len(self.networks)):
			for j in range(len(self.networks[i])):
				if(outputs):
					for k in range(self.networks[i][j].mat.shape[1]):
						mask = np.array([False if l ==k and self.ignoreSelfReg else True for l in range(self.networks[i][j].mat.shape[1])])
						neworder=np.arange(self.networks[i][j].mat.shape[0])[mask]
						np.random.shuffle(neworder)
						self.networks[i][j].mat[mask,k] = self.networks[i][j].mat[neworder, k]
						self.networks[i][j].reducedmat[mask,k] = self.networks[i][j].reducedmat[neworder, k]
				if(inputs):
					for k in range(self.networks[i][j].mat.shape[0]):
						mask = np.array([False if l ==k and self.ignoreSelfReg else True for l in range(self.networks[i][j].mat.shape[0])])
						neworder=np.arange(self.networks[i][j].mat.shape[1])[mask]
						np.random.shuffle(neworder)
						self.networks[i][j].mat[k, mask] = self.networks[i][j].mat[k, neworder]
						self.networks[i][j].reducedmat[k, mask] = self.networks[i][j].reducedmat[k, neworder]
				self.networks[i][j] = subnetwork(self.networks[i][j].mat, self.networks[i][j].reducedmat, self.networks[i][j].originalmat, self.networks[i][j].rownames, self.networks[i][j].colnames, np.sum(self.networks[i][j].mat, axis = 1), np.sum(self.networks[i][j].mat, axis = 0), self.networks[i][j].tipo, self.networks[i][j].simid, self.networks[i][j].cell)


class motifContainerHash(motifContainer):
	def __init__(self, simulSet, motifs = None, type_to_separate = 1,type_to_separate_motifs=3, compress_terminal_features = True, addNewMotifs=True, preprocessedNets=None, motifLastId = None, ignoreSelfReg = False):
		networkInstances = simulSet.simulations
		super().__init__(simulSet, motifs, type_to_separate, type_to_separate_motifs, compress_terminal_features, addNewMotifs, preprocessedNets, ignoreSelfReg = ignoreSelfReg)
		self.isplain = False 
		if(motifs is None):
			self.motifs = {}
			self.motif_lastID = {}
			for i in range(min(networkInstances[0].genome_size, self.MAX_MOT_SIZE)):
				self.motifs[i] = {}
				self.motif_lastID[i] = 0
		else:
			assert(motifLastId is not None), "need to provide a dictionary with last IDs for motifs"
			self.motifs = motifs
			self.motif_lastID = motifLastId
	def findMotifs(self, starting_n = 0, selectBy=0, functionalFilter=True):
		assert not self.isplain, "motif hash is plain, use super method"
		start_time = time.time()
		sim_cont= 0
		for sim in self.networks:
			cell_cont=0
			for subnet in sim:
				subnet_generator = subNetworkGenerator(subnet, selectBy=selectBy, type_to_separate=self.type_to_separate_motifs, starting_n=starting_n)
				try:
					while(not subnet_generator.complete): 
						sub = subnet_generator.next()
						if(functionalFilter):
							if(not self.functionalSubnetwork(sub)):
								continue
						if(sub is None):
							break
						n = sub.mat.shape[0]
						signature = ','.join(str(sub.colsum))
						sub_present = False
						motsig_list = self.motifs[n].get(signature, [])
						for mot in motsig_list:
							sub_present = mot.compare(sub, simset=self.set_name, memorizeSubnet = self.addNewMotifs)
							if(sub_present):
								break
						if(not sub_present and self.addNewMotifs):
							newMotif = Motif(sub, self.type_to_separate_motifs, simset=self.set_name, ID = 's' + str(n) + 'mot' + str(self.motif_lastID[n]))
							self.motif_lastID[n]+=1
							motsig_list.append(newMotif)
							self.motifs[n][signature] = motsig_list
				except:
					print("Error while searching for motifs: sim ", sim_cont, " cell ", cell_cont, " ", self.set_name,"\n")
					print(sys.exc_info())
					return subnet_generator
					
				sim_cont += 1
			cell_cont+=1
		return time.time() - start_time
	def writeMotifSummary(self, fname = None, outpath = ''):
		if(self.isplain):
			super().writeMotifSummary(fname, outpath)
		else:
			mots_reserva = self.motifs
			self.motifs = self.getPlainMotifList()
			self.isplain = True
			super().writeMotifSummary(fname, outpath)
			self.motifs = mots_reserva
			self.isplain = False
			file.close()
	def getPlainMotifList(self):
		if(self.isplain):
			return self.motifs
		plainmotifs = {}
		for i in self.motifs.keys():
			plainmotifs[i]=[]
			for k in self.motifs[i].keys():
				plainmotifs[i].extend(self.motifs[i][k])
		return plainmotifs
	def becomePlain(self):
		self.motifs = self.getPlainMotifList()
		self.isplain = True
	def __iter__(self):
		for i in self.motifs.keys():
			if(not self.isplain):
				d = []
				for subkey in self.motifs[i].keys():
					d.extend(self.motifs[i][subkey])
			else:
				d = self.motifs[i]
			for mot in d:
				yield mot
	def makeTables(self, addToName = '', min_matches=40):
		if(self.isplain):
			super().makeTables()
		else:
			mots_reserva = self.motifs
			self.motifs = self.getPlainMotifList()
			self.isplain = True
			t = super().makeTables(addToName, min_matches)
			self.motifs = mots_reserva
			self.isplain = False
		return t
	def __setPermanentIds(self):
		assert self.isplain, "not plain, can't perform action"
		super().__setPermanentIds()
	def resetInstances(self):
		if(self.isplain):
			super().resetInstances()
		else:
			for n in self.motifs.keys():
				for sig in self.motifs[n].keys():
					for m in self.motifs[n][sig]:
						m.resetInstances()

class randomizedContainerHash(motifContainerHash):
	def __init__(self, container, addName = "random_0", randomizeNodeOutputs = True, randomizeNodeInputs = False, addNewMotifs = False, cpmotifs = False):
		self.path = container.path
		self.set_name = container.set_name + addName
		self.type_to_separate = container.type_to_separate
		self.type_to_separate_motifs = container.type_to_separate_motifs
		self.networks = cp.deepcopy(container.networks)
		self.ignoreSelfReg = container.ignoreSelfReg
		self.randomizeNetworks(randomizeNodeOutputs, randomizeNodeInputs)
		if(cpmotifs):
			self.motifs = cp.deepcopy(container.motifs)
		else:
			self.motifs = container.motifs
		self.motif_lastID = container.motif_lastID
		self.isplain = container.isplain
		self.addNewMotifs = addNewMotifs
		self.original_container_threshold_bin = container.threshold_bin
		self.functional_definition = functional_definition
	randomizeNetworks = randomizedContainer.__dict__['randomizeNetworks']


class containerRandomizer():
	def __init__(self, container, randsize = 10, randomizeNodeOutputs = True, randomizeNodeInputs = False, addNewMotifs = False):
		self.set_name = container.set_name + 'random_all'
		self.randsize = randsize
		self.original_container = container
		self.addNewMotifs = addNewMotifs
		if(type(container) == motifContainer):
			self.randomclass = randomizedContainer
		elif(type(container) == motifContainerHash):
			self.randomclass = randomizedContainerHash
		self.randomizedContainers = [self.randomclass(container, addName = "random_"+str(i), randomizeNodeOutputs = randomizeNodeOutputs, randomizeNodeInputs = randomizeNodeInputs, addNewMotifs = addNewMotifs, cpmotifs = False) for i in range(randsize)]
	def findMotifsPar(self, starting_n = 0, selectBy=0, functionalFilter=True):
		import multiprocessing
		random_results = multiprocessing.Queue()
		random_times = multiprocessing.Queue()
		jobs = [multiprocessing.Process(target = self.__findMotifsSingleContainer, args = (rcont, starting_n, selectBy, functionalFilter, random_times, random_results)) for rcont in self.randomizedContainers]
		for j in jobs:
			j.start()
		self.randomizedContainers = [random_results.get() for j in jobs]
		self.mergeMotifCountsToOriginalContainer()
		return [random_times.get() for j in jobs]
	def __findMotifsSingleContainer(self,cont, starting_n, selectBy, functionalFilter, random_times, random_results):
		t=cont.findMotifs(starting_n, selectBy, functionalFilter)
		random_times.put(t)
		random_results.put(cont)
	def findMotifs(self, starting_n = 0, selectBy=0, functionalFilter=True):
		times = []
		for rcont in self.randomizedContainers:
			t= rcont.findMotifs(starting_n, selectBy, functionalFilter)
			times.append(t)
		self.mergeMotifCountsToOriginalContainer()
		return times
	def mergeMotifCountsToOriginalContainer(self):
		omot_hash = {o.permanentId:o for o in self.original_container}
		for rcontainer in self.randomizedContainers:
			for rc in rcontainer:
				try:
					omot_hash[rc.permanentId].matches = {**omot_hash[rc.permanentId].matches, **rc.matches}
				except KeyError:
					omot_hash[rc.permanentId] = rc	#this makes nothing actually
		
			
class multiSetMotifFinder:
	def __init__(self, path, conditions, type_to_separate = 1, type_to_separate_motifs=3,compress_terminal_features = True, starting_n = 3, selectBy=0, functionalFilter = True,randsize = 10, randomizeNodeOutputs = True, randomizeNodeInputs = False, addNewMotifs = False, cpmotifs = False, hashSearch = True, randomParallel = True, ignoreSelfReg = False):
		self.paths= {i:''.join((path, i)) for i in conditions}
		self.type_to_separate = type_to_separate	#This affects sorting of nodes in subnetwork; they are subdivided according to type[type_to_separate]. If compress_terminal_features, terminal features are summed
		self.type_to_separate_motifs= type_to_separate_motifs	#this affects comparison between subnetwork and motifs. type[type_to_separate_motifs] must be the same as in Motif for every node in subnetwork
		self.compress_terminal_features = compress_terminal_features # Whether to sum up connexions to terminal features of one type, according to type[type_to_separate]
		self.starting_n = starting_n	#Maximum size of motif
		self.selectBy= selectBy		#When iterating over subnetworks, subnetworks that don't have at least one gene of each kind in restrict[selectBy], are discarded
		self.functionalFilter = functionalFilter #When searching, if True, subnetwork that don't include an input and an output for every TF are discarded (lineage TFs are only required to have outputs)
		self.randsize = randsize	#number of random sets per original set
		self.randomizeNodeOutputs = randomizeNodeOutputs #randomize outputs in random sets or not
		self.randomizeNodeInputs = randomizeNodeInputs	#randomize inputs in random sets or not
		self.addNewMotifs = addNewMotifs		#whether to add new motifs when searching in random networks. If multiprocessing with random, it does not work
		self.cpmotifs = cpmotifs			#when creating random network datasets, whether to use deepcopy to copy motifs or not
		self.randomParallel = randomParallel
		self.ignoreSelfReg = ignoreSelfReg
		self.mc = None
		self.ss = None
		if(hashSearch):
			self.containerClass = motifContainerHash
		else:
			self.containerClass = motifContainer
	def fullMotifAnalysis(self, extraname = '', min_matches=40, save_mc = False, save_adex = False):
		import simulationAnalyzerServer as sana		### remember to check this before sbmitting!
		motifs = None
		last_mot_id = None
		for cond in self.paths.keys():
			cond_time = time.time()
			print(cond, ' started...')
			self.ss= sana.simulationSet(self.paths[cond], cond)
			self.mc = self.containerClass(self.ss, motifs = motifs, type_to_separate = self.type_to_separate, type_to_separate_motifs=self.type_to_separate_motifs,compress_terminal_features = self.compress_terminal_features, addNewMotifs=True, preprocessedNets=None, motifLastId=last_mot_id, ignoreSelfReg = self.ignoreSelfReg)
			print('\t', cond, ' read')
			t = self.mc.findMotifs(self.starting_n, self.selectBy, self.functionalFilter)
			print('\t',cond, ' motifs found in ', str(t))
			rc = containerRandomizer(self.mc, randsize=self.randsize, randomizeNodeOutputs = self.randomizeNodeOutputs, randomizeNodeInputs = self.randomizeNodeInputs, addNewMotifs = self.addNewMotifs)
			print('\t',cond, ' random sets generated')
			if(self.randomParallel and not self.addNewMotifs):
				t=rc.findMotifsPar(self.starting_n, self.selectBy, self.functionalFilter)
			else:
				if(self.addNewMotifs):
					print(cond, " Warning!: Not using multiprocessing because addNewMotifs == True")
				t=rc.findMotifs(self.starting_n, self.selectBy, self.functionalFilter)
			print('\t',cond, ' random sets analyzed in ', str(t))
			try:
				t = self.mc.makeTables(extraname, min_matches)
				print('\t',cond, ' tables made in ', str(t))
			except:
				print("Some error making tables")
			finally:
				if(save_mc):
					self.mc.saveSelf(extraname)
					print('\t',cond, ' saved')
			adex = auxiliaryDataExtractor(self.mc, self.ss)
			t = adex.makeTables(extraname)
			if(save_adex):
					adex.saveSelf(extraname)
			print('\t',cond, ' auxiliary tables made in ', str(t))
			self.mc.resetInstances()
			print('\t', cond, ' motifs resetted\n')
			motifs = self.mc.motifs
			if(self.containerClass is motifContainerHash):
				last_mot_id = self.mc.motif_lastID
			print('condition ', cond, 'finished in: ', time.time() - cond_time, '\n***\n')

class auxiliaryDataExtractor():
	def __init__(self, net_container, simul_set):
		self.netcontainer = net_container
		self.simset = simul_set
		self.stored_results = {}
	def makeTables(self, extraname = ''):
		start_time = time.time()
		self.makeTabSkeleton()
		self.getExpressionBySim()
		#print('a')
		self.getMutantPhenotypeByCell('tfs')
		#print('b')
		self.getMutantPhenotypeByCell('sites')
		#print('c')
		self.getNumberOfRegulatorsByGeneByCell()
		#print('d')
		self.getMeanPhenotypeByGeneByCell('tfs')
		#print('e')
		self.getMeanPhenotypeByGeneByCell('sites')
		#print('f')
		self.getCorrelationBetweenTerminalTypes()
		#print('g')
		self.getOverlapOfRegulationSets()
		#print('h')
		self.writeTables(extraname)
		#print('i')
		return time.time() - start_time
	def makeTabSkeleton(self):
		groupers= ['gene_unique_id','simset', 'simid', 'gene', 'type3'] #'type4' removed
		try:
			self.stored_results['skeleton'] = self.netcontainer.summaryTables['tfs_cells_sep'][groupers]	
			self.stored_results['skeleton']= self.stored_results['skeleton'].assign(type4=list(self.netcontainer.summaryTables['type4']))
		except KeyError:
			print('makeTabSkeleton KeyError - No tfs_cells_sep')
	def getExpressionBySim(self):
		# Select only TF genes
		df = self.stored_results['skeleton'].loc[self.stored_results['skeleton']['type3'] == type1define['tf']]
		df = df.drop_duplicates()	#in skeleton table there is a row per gene per cell; now only want a row per cell
		#mat = np.zeros([self.stored_results['skeleton'].shape[0], self.simset.simulations[0].ncells]) # store results
		mat = np.zeros([df.shape[0], self.simset.simulations[0].ncells])
		for sim in self.simset.simulations:
			df_index = [i for i,j in enumerate(df.simid) if j == re.sub('_$' ,'',sim.name)]
			mat[df_index, :] = sim.expr[[i for i in sim.expr.columns if 'exp' in i]].loc[list(df['gene'].iloc[df_index])]
		mat = pd.DataFrame(mat, columns = [i for i in sim.expr.columns if 'exp' in i])
		mat_sorted = np.array(mat)
		mat_sorted.sort(axis=1)
		mat_sorted = pd.DataFrame(mat_sorted[:,::-1], columns = ['sorted'+i for i in sim.expr.columns if 'exp' in i])
		self.stored_results['expression_uniqueGeneIds'] = pd.concat([df.reset_index(drop=True), mat, mat_sorted], axis = 1)			
	def getMutantPhenotypeByCell(self, mutation_type = 'tfs'):		
		df = self.stored_results['skeleton'].loc[self.stored_results['skeleton']['type3'] == type1define['tf']]
		df = df.drop_duplicates()
		typeset = sorted(set(self.simset.simulations[0].type3))
		mat = np.zeros([df.shape[0], self.simset.simulations[0].genome_size*self.simset.simulations[0].ncells])
		t3 = self.simset.simulations[0].type3 #outer name (more general and invariant between cells)
		t = self.simset.simulations[0].type #inner name
		mat_colnames = ['_'.join(['cell'+str(c), t3[g], t[g], str(g)]) for tt in typeset for c in range(self.simset.simulations[0].ncells) for g in range(self.simset.simulations[0].genome_size)  if t3[g]==tt]
		for sim in self.simset.simulations:
			simphen =np.concatenate([sim.getSortedMutEffectsByMean(type3 = [i], mutation_type = mutation_type)[0] for i in typeset], axis=1)
			df_index = [i for i,j in enumerate(df.simid) if j == re.sub('_$' ,'',sim.name)]
			mat[df_index, :] = simphen[list(df['gene'].iloc[df_index]),:]
		self.stored_results[mutation_type + 'effect_uniqueGeneIds'] = pd.concat([df.reset_index(drop=True), pd.DataFrame(mat, columns = mat_colnames)], axis=1)		
	def getNumberOfRegulatorsByGeneByCell(self):
		diftypes = sorted(set(self.netcontainer.networks[0][0].tipo[2][self.netcontainer.networks[0][0].tipo[3]==type1define['tf']]))
		tuplas = [[[re.sub('_$' ,'',net.simid)]*net.rownames.shape[0],[net.cell]*net.rownames.shape[0], net.rownames] + [np.concatenate([np.sum(net.mat[:, net.tipo[2][net.colnames]==t], axis=1).reshape([net.colnames.shape[0],1]) for t in diftypes], axis=1)] for i in self.netcontainer.networks for net in i]
		df= self.netcontainer.summaryTables['tfs_cells_sep'][['gene_unique_id', 'simset', 'simid', 'cell', 'gene', 'type0', 'type1','type2', 'type3', 'number']]
		## Check that data can be directly concatenated:
		t2 = np.concatenate([t[2] for t in tuplas], axis = 0)
		t3 = np.array([i for t in tuplas for i in t[0]])
		if(np.all(t2 == df.gene) & np.all(t3 ==  df.simid)):
			inv_map = {v: k for k, v in type1define.items()}
			df2 = pd.DataFrame(np.concatenate([tup[3] for tup in tuplas], axis = 0) , columns = [inv_map[t] for t in diftypes])
			df2 = pd.concat([df, df2], axis = 1)
		else:
			#Not equal! make loop step by step
			print("not equal order in getNumberOfRegulatorsByGeneByCell")
			#Not implemented
		df2['all_activators'] = df2.lineage_this_cell + df2.lineage_other_cell + df2.other_activator
		self.stored_results['regulatorsByGeneByCell']  = df2
		self.stored_results['regulatorsByGeneByCell_means']  = df2.groupby('type2').mean()
	def getMeanPhenotypeByGeneByCell(self, mut_type = 'tfs'):
		#df= self.netcontainer.summaryTables['tfs_cells_sep']['gene_unique_id', 'simset', 'simid', 'cell', 'gene', 'type0', 'type1','type2', 'type3', 'number']
		regtypes = sorted(set(self.netcontainer.networks[0][0].tipo[2][self.netcontainer.networks[0][0].tipo[3]==type1define['tf']]))
		tipo_by_cell = [self.netcontainer.networks[0][cell].tipo for cell in range(len(self.netcontainer.networks[0]))]
		inv_map = {v: k for k, v in type1define.items()}
		mean_phenotype = []
		for sim in self.simset.simulations:
			if(mut_type == 'tfs'):
				mut = sim.mutanttfs
				mutname = 'tf_mutated'
			else:
				mut = sim.mutantsites
				mutname = 'mutated_sites'
			wt=sim.expr.iloc[mut.index] #here .iloc is equivalent to .loc
			phenotype = mut[[i for i in mut.columns if 'exp' in i]] - wt[[i for i in wt.columns if 'exp' in i]]
			phenotype[mutname] = mut[mutname]
			#phenotype = phenotype.reset_index(drop=False)
			for cell in range(sim.ncells):
				for g in set(phenotype.index):
					aux_a =  [re.sub('s[0-9]+_*' ,'',sim.name), re.sub('_$' ,'',sim.name), cell, g] + [t[g] for t in tipo_by_cell[cell]]
					aux_b = []
					this_gene= phenotype.loc[g][['exp' + str(cell), mutname]].iloc[list(phenotype.loc[g]['exp' + str(cell)] != 0)]
					for r in regtypes:
						rr = this_gene.iloc[tipo_by_cell[cell][2][list(this_gene[mutname])] == r]
						if(rr.shape[0] > 0):
							aux_b.append(np.mean(np.array(rr[[i for i in rr.columns if 'exp' in i]])))
						else:
							aux_b.append(np.nan)
					all_activators = [True  if xx in [type1define['other_activator'], type1define['lineage_this_cell']] else False for xx in tipo_by_cell[cell][1][list(this_gene[mutname])]]
					rr = this_gene.iloc[all_activators]
					if(rr.shape[0] > 0):
						aux_b.append(np.mean(np.array(rr[[i for i in rr.columns if 'exp' in i]])))
					else:
						aux_b.append(np.nan)
					mean_phenotype.append(aux_a + aux_b)
		colnames = ['simset', 'simid', 'cell', 'gene'] + ['type' + str(i) for i in range(len(tipo_by_cell[0]))]	+ [inv_map[i] for i in regtypes] + ['all_activators']
		df= pd.DataFrame(mean_phenotype, columns = colnames)
		self.stored_results['meanPhenotypeByGene_' + mut_type]  = df
		self.stored_results['meanPhenotypeByGene_' + mut_type + '_means']  = df.groupby('type2').mean()
				
	def getCorrelationBetweenTerminalTypes(self):
		import itertools
		tipo_by_cell = [self.netcontainer.networks[0][cell].tipo for cell in range(len(self.netcontainer.networks[0]))]
		regtypes = sorted(set(self.netcontainer.networks[0][0].tipo[2][self.netcontainer.networks[0][0].tipo[3]==type1define['tf']]))
		tertypes = sorted(set(self.netcontainer.networks[0][0].tipo[2][self.netcontainer.networks[0][0].tipo[3]==type1define['non_tf']]))
		regnames = np.where(tipo_by_cell[0][3] == type1define['tf'])[0]
		ternames = np.where(tipo_by_cell[0][3] == type1define['non_tf'])[0]
		inv_map = {v: k for k, v in type1define.items()}
		typepairs = list(itertools.combinations_with_replacement(tertypes, 2))
		colnames = ['simset', 'simid', 'cell'] + ['Act_' + inv_map[a] + ':' + inv_map[b] if a != b else 'Act_' + inv_map[a] + '_selfcor' for a, b in typepairs] + ['Inh_' + inv_map[a] + ':' + inv_map[b] if a != b else 'Inh_' + inv_map[a] + '_selfcor' for a, b in typepairs]
		tuplas = []
		for sim in self.simset.simulations:
			for cell in range(sim.ncells):
				this_sim = [re.sub('s[0-9]+_*' ,'',sim.name),re.sub('_$' ,'',sim.name), cell]
				this_sim_inh = []
				mat=np.zeros([ternames.shape[0], regnames.shape[0]])
				mat_inh=np.zeros([ternames.shape[0], regnames.shape[0]])
				#Get regulators for each gene
				for ter_ind, ter in enumerate(ternames):
					mat[ter_ind, sim.getActivators(target=[ter], cell=cell)] = 1	
					mat_inh[ter_ind, sim.getInhibitors(target=[ter], cell=cell)] = 1
				mat_cor = np.zeros([ternames.shape[0], ternames.shape[0]])
				mat_cor_inh = np.zeros([ternames.shape[0], ternames.shape[0]])
				#calculate correlations between each pair of genes	
				for ter_ind in range(mat.shape[0]):
					for ter_ind2 in range(ter_ind):
						y = np.where(mat[ter_ind,:] + mat[ter_ind2,:] > 0)[0].shape[0]
						y_inh = np.where(mat_inh[ter_ind,:] + mat_inh[ter_ind2,:] > 0)[0].shape[0]
						if (y == 0):
							mat_cor[ter_ind, ter_ind2] = np.nan
						else:					
							mat_cor[ter_ind, ter_ind2] = np.where(mat[ter_ind,:] + mat[ter_ind2,:] == 2)[0].shape[0]/y
						if (y_inh == 0):
							mat_cor_inh[ter_ind, ter_ind2] = np.nan
						else:
							mat_cor_inh[ter_ind, ter_ind2] = np.where(mat_inh[ter_ind,:] + mat_inh[ter_ind2,:] == 2)[0].shape[0]/y_inh
				mat_cor = mat_cor + np.transpose(mat_cor)
				mat_cor_inh = mat_cor_inh + np.transpose(mat_cor_inh)
				for a, b in typepairs:
					if(a != b):
						mean_cor = np.mean(mat_cor[np.where(tipo_by_cell[cell][2][ternames] == a)[0],:][:, np.where(tipo_by_cell[cell][2][ternames] == b)[0]])
						mean_cor_inh = np.mean(mat_cor_inh[np.where(tipo_by_cell[cell][2][ternames] == a)[0],:][:, np.where(tipo_by_cell[cell][2][ternames] == b)[0]])
					else:
						aux_mat = mat_cor[np.where(tipo_by_cell[cell][2][ternames] == a)[0],:][:, np.where(tipo_by_cell[cell][2][ternames] == b)[0]]
						aux_mat_inh = mat_cor_inh[np.where(tipo_by_cell[cell][2][ternames] == a)[0],:][:, np.where(tipo_by_cell[cell][2][ternames] == b)[0]]
						mask = np.ones(aux_mat.shape, dtype=bool)
						np.fill_diagonal(mask, 0)
						mean_cor = np.mean(aux_mat[mask])
						mean_cor_inh = np.mean(aux_mat_inh[mask])
					this_sim.append(mean_cor)
					this_sim_inh.append(mean_cor_inh)
				tuplas.append(this_sim + this_sim_inh)
		df = pd.DataFrame(tuplas, columns = colnames)
		self.stored_results['regCorrelationsBetweenTypes']  = df
	def getOverlapOfRegulationSets(self):
		import itertools
		tipo_by_cell = [self.netcontainer.networks[0][cell].tipo for cell in range(len(self.netcontainer.networks[0]))]
		regtypes = sorted(set(self.netcontainer.networks[0][0].tipo[2][self.netcontainer.networks[0][0].tipo[3]==type1define['tf']]))
		tertypes = sorted(set(self.netcontainer.networks[0][0].tipo[2][self.netcontainer.networks[0][0].tipo[3]==type1define['non_tf']]))
		regnames = np.where(tipo_by_cell[0][3] == type1define['tf'])[0]
		ternames = np.where(tipo_by_cell[0][3] == type1define['non_tf'])[0]
		inv_map = {v: k for k, v in type1define.items()}
		typepairs = list(itertools.permutations(tertypes, 2))
		colnames = ['simset', 'simid', 'cell'] + ['Act_' + inv_map[a] + ':' + inv_map[b] for a, b in typepairs] + ['Inh_' + inv_map[a] + ':' + inv_map[b] for a, b in typepairs]
		tuplas = []
		for sim in self.simset.simulations:
			for cell in range(sim.ncells):
				this_sim = [re.sub('s[0-9]+_*' ,'',sim.name),re.sub('_$' ,'',sim.name), cell]
				this_sim_inh = []
				for ta, tb in typepairs:
					genes_a = [g for g in ternames if tipo_by_cell[cell][2][g] == ta]
					genes_b = [g for g in ternames if tipo_by_cell[cell][2][g] == tb]
					aa=[]
					ai=[]
					ba=[]
					bi=[]
					for g in genes_a:
						aa = aa +sim.getActivators(target=[g], cell=cell)
						ai = ai + sim.getInhibitors(target=[g], cell=cell)
					aa= set(aa)
					ai = set(ai)
					for g in genes_b:
						ba = ba +sim.getActivators(target=[g], cell=cell)
						bi = bi + sim.getInhibitors(target=[g], cell=cell)
					ba= set(ba)
					bi = set(bi)
					if(len(aa)>0):
						overlap_act = len(aa.intersection(ba))/len(aa)
					else:
						overlap_act = np.nan
					if(len(ai)>0):
						overlap_inh = len(ai.intersection(bi))/len(ai)
					else:
						overlap_inh = np.nan
					this_sim.append(overlap_act)
					this_sim_inh.append(overlap_inh)
				tuplas.append(this_sim + this_sim_inh)
		df = pd.DataFrame(tuplas, columns = colnames)
		self.stored_results['regIntersectionBetweenTypes']  = df		
	def writeTables(self, extraname = ''):
		#simset = re.sub('s[0-9]+_*' ,'',self.simset.simulations[0].name) #should be the same as current line
		if(extraname != ''):
			extraname = '_' + extraname
		simset =  self.netcontainer.set_name  + extraname
		for name, table in self.stored_results.items():
			if('skeleton' in name):
				continue
			table.to_csv(simset + '_' + name + '.csv')
		#self.saveSelf(extraname)
	def saveSelf(self, name=''):
		import _pickle as pck
		filename = self.netcontainer.set_name + name +'_extraDataObject.pck'
		with open(filename, 'wb') as output: 
 			pickler = pck.Pickler(output, -1)
 			pickler.dump(self)
	@staticmethod
	def load(filenamein):
		with open(filenamein, 'rb') as f:
			print('opening ' + filenamein)
			import _pickle as pck
			dex = pck.load(f)
		return dex						
def main():
	#default_path = '../cellevolver_shared/'
	#all_conditions = ['4cell_mce0', '4cell_mce0fix_mutb', '4cell_mce1fix', '4cell_mce2fix', '4cell_mce2inh5', '4cell_mce0fix', '4cell_mce0_mutb', '4cell_mce1fix_mutb', '4cell_mce2fix_mutb', '4cell_mce2inh5_mutb', '4cell_mce0Xss']
	default_path= './sim_tables/'
	all_conditions=['4cell_mce4fix', '4cell_mce5fix', '4cell_mce8fix', '4cell_mce9fix', '6cell_mce6fix', '6cell_mce7fix']

	path = default_path
	conditions = all_conditions
	mm = multiSetMotifFinder(path, conditions,  type_to_separate = 1, type_to_separate_motifs=3, selectBy=0, starting_n = 3, functionalFilter = True, hashSearch=True, randsize=10, ignoreSelfReg = True)
	mm.fullMotifAnalysis(extraname = sys.argv[1])


if __name__ == "__main__":
    main()
