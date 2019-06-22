import geneExpression as gx
import numpy as np
import sys
from scipy.stats import gmean
import abc
import copy
import _pickle as pck 
#import cPickle as pck

'''Read files with .cemod extension and produce a dictionary with simulation parameters'''

a=1+1
class modelReader(object):
	def readModelFile(self, model_def_file = None):
		assert(model_def_file.endswith('.cemod'))
		import re
		d = {}
		with open(model_def_file) as f:
			for line in f:
				line = line.rstrip('\n')
				if(line[0] != '#' and line != ''):
					(key, val) = line.split(':')
					val = self.__processLineVal(key, val)
					d[key] = val
		keys = d.keys()
		for k in reversed(list(d.keys())):
			d[k]= self.__elementReplace(d[k], d)
		for k in list(d.keys()):
			if(k[0] == '$' or k[0] == '_'):
				del d[k]
		return d
	def __elementReplace(self, el, source):
		if(type(el) is list):
			for k in range(len(el)):
				if(type(el[k]) is not str):
					continue
				elif(el[k][0] == '_' or el[k][0] == '$'):
					el[k] = source[el[k]]
		elif(type(el) is dict):
			for k in el.keys():
				if(type(el[k]) is not str):
					continue
				if(el[k][0] == '_' or el[k][0] == '$'):
					el[k] = source[el[k]]
		elif(type(el) is str):
			if(el[0] == '_' or el[0] == '$'):
				el = source[el]
		else:
			pass
		return el
	def __processLineVal(self, key, val):
		if(key[0] != '$' and key[0] != '_'):
			val = self.__setStringValue(val)
		elif(key[0] == '$'):
			val = [self.__setStringValue(i) for i in val.split()]
		elif(key[0] == '_'):
			newdict = {}
			for s in val.split(','):
				(newkey, newval) = s.split()
				newval = self.__setStringValue(newval)
				newdict[newkey] = newval
			val = newdict
		return val
	def __setStringValue(self, val):
		if(val.isdigit()):
			val = int(val)	
		elif(self.__is_numeric(val)):
			val = float(val)
		elif(type(val) is str):
			if(val == 'True'):
				val = True
			elif(val == 'False'):
				val = False
			elif(val == 'None'):
				val = None
		return val
	@staticmethod
	def __is_numeric(s):
	    try:
	        float(s)
	        return True
	    except (ValueError, TypeError):
	        return False

'''This class is only to evolve 1 DNA sequence; it works fine but it is useless for the rest of the model'''
class singleSeqEvolver():
	def __init__(self, N, length, background):
		self.N = N
		self.length = length
		self.a_priori = background
		self.seqs = self.generateRandomSeqs(N,length, background)
	def generateRandomSeqs(self, N,length, bck):
		alphabet = list(bck.keys())
		probs = list(bck.values())
		letters = np.random.choice(alphabet, N*length, replace=True, p = probs).reshape(N, length)
		seqs  = [''.join(r) for r in letters]
		return np.array(seqs)
	def newGeneration(self, probs, mut_rate = 0.1, exp = 2):
		assert(len(probs) == len(self.seqs))
		if(exp>1):
			probs = probs**exp
		probs = probs/(np.sum(probs))
		newseqs = np.random.choice(self.seqs, len(self.seqs), replace = True, p = probs)
		if(mut_rate > 0):
			newseqs = self.mutateSeqs(newseqs, mut_rate)
		self.seqs = newseqs
	def mutateSeqs(self, seqs, mut_rate = 0.25):
		charlist = [list(s) for s in seqs]

		mutate = np.random.choice(range(0, self.N*self.length), round(mut_rate*self.N*self.length), replace=False)

		rows = np.floor_divide(mutate, self.length)
		cols = np.mod(mutate, self.length)

		alphabet = list(self.a_priori.keys())
		probs = list(self.a_priori.values())

		newchars = np.random.choice(alphabet, len(mutate), replace=True, p = probs)
		for i in range(len(mutate)):
			charlist[rows[i]][cols[i]] = newchars[i]
		return np.array([''.join(r) for r in charlist])

''' Uses singleSeqEvolver to evolve a DNA sequence'''
def evolveSingleSeq():
	N = 6
	length = 600
	bck = {'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25}
	tfs = gx.TFSet(range(0, 10), alphasbtm = 3, interactions = [gx.interactionTuple(0, 1, 3), gx.interactionTuple(2, 3, 3)])
	ex = gx.ExpressionCalculator()
	ss = singleSeqEvolver(N, length, bck)
	ann=tfs.annotate(ss.seqs)
	final_value = 2

	e1 = np.zeros(len(ann))
	genmean = 0
	allExpr = []
	gen = 0
	while(genmean < final_value):
		for i in range(len(ann)):
			e1[i] = ex.activate(ann[i], interactions = True)
		ss.newGeneration(probs=e1, mut_rate = 0.005, exp=4)
		allExpr.append(e1)
		ann=tfs.annotate(ss.seqs)
		genmean = np.mean(e1)
		print('gen ', gen, ' mean: ', genmean, ' sd: ', np.std(e1),'\n')
		gen+=1





class CellBasic():
	def __init__(self, cellname, subtype, initial_concentrations):
		self.name = cellname
		self.subtype = subtype
		self.initial_concentrations = initial_concentrations

'''This is a template class which implementations should specify a way to evaluate an organism fitness'''
class PhenotypeEvaluator():
	def __init__(self, args):
		self.args = args
		print('Instance of abstract class created')		
	@abc.abstractmethod
	def generateCells():
		return None
	@abc.abstractmethod
	def evaluateOrganism():
		return None 

'''In this class I define a way to evaluate the fitness of individuals based on the squared error between a predetermined level of expression in each cell and the actual level of expression.
When an instance is created, a vector of cells is generated, each with its predetermined optimal expression levels. This class depends strongly on the OrganismTemplate class'''
class SimplePhenotypeEvaluator(PhenotypeEvaluator):
	INITIAL_CONCENTRATIONS = 1
	OPTIMAL_CONCENTRATIONS = 3
	def __init__(self, model_file = None, template = None, args = None):
		self.organism_template = template
		self.pan_genes, self.specific_genes, self.ntypes, self.cells_per_type, cellstrings = args
		self.cells = self.generateCells(cellstrings)
	def generateCells(self, cell_strings = None):
		cells = []
		if(cell_strings is None):
			assert(self.organism_template.tf_num + self.pan_genes + self.specific_genes*self.ntypes <= self.organism_template.genome_size and len(self.organism_template.tf_lineage) >= self.organism_template.ncells)
			pan_genes_index = range(self.organism_template.tf_num, self.organism_template.tf_num + self.pan_genes)
			cell_id = 0
			specific_so_far = 0
			for t in range(len(self.cells_per_type)):
				for i in range(self.cells_per_type[t]):
					initial = np.zeros(self.organism_template.genome_size)
					initial_vals = [self.INITIAL_CONCENTRATIONS if i == cell_id else 0 for i in self.organism_template.tf_lineage]
					initial[self.organism_template.tf_lineage] = initial_vals
					optimal = np.zeros(self.organism_template.genome_size)
					optimal[pan_genes_index] = self.OPTIMAL_CONCENTRATIONS
					if(self.specific_genes > 0):
						specific_genes = range(self.organism_template.tf_num + self.pan_genes + specific_so_far, self.organism_template.tf_num + self.pan_genes + specific_so_far + self.specific_genes)
						optimal[specific_genes] = self.OPTIMAL_CONCENTRATIONS
						specific_so_far += self.specific_genes #cells of the same kind will have the same specific genes and differet lineage_tf
					cells.append(self.SimpleCell(cell_id, t, initial ,optimal))
					cell_id +=1
				
		else:
			for c in range(sum(self.cells_per_type)):
				subt = cell_strings[2][c]
				initial = np.array(cell_strings[0][c])
				optimal = np.array(cell_strings[1][c])
				cells.append(self.SimpleCell(cellname = c, subtype = subt, initial_concentrations = initial,  optimal_concentrations = optimal))
		return cells
	def getInitialConcentrations(self, i):
		return self.cells[i].initial_concentrations
	def evaluateOrganism(self, organism):
		#assert(self.organism_template == organism.template)
		error = np.zeros(self.organism_template.ncells)
		for cell in range(organism.expression.shape[1]):
			error[cell] = self.squared_error(self.cells[cell].optimal_concentrations[self.organism_template.tf_num:], organism.expression[self.organism_template.tf_num:, cell])
		#return gmean(error + 0.0001)
		return np.mean(error)
		#return np.max(error)
	@staticmethod
	def squared_error(a, b):
		a = a.reshape(a.shape[0])
		b = b.reshape(b.shape[0])
		err = np.sum(np.square(a - b))/a.shape[0]
		return err
	class  SimpleCell(CellBasic):
		def __init__(self, cellname, subtype, initial_concentrations, optimal_concentrations):
			super().__init__(cellname, subtype, initial_concentrations)
			self.optimal_concentrations = optimal_concentrations
		 
'''OrganismTemplate implements these methods'''
class AbstractOrganismTemplate():
	@abc.abstractmethod	
	def evaluateOrganism():
		return None
	@abc.abstractmethod
	def generateRandomSeqs(N, length, bck):
		return None
	@abc.abstractmethod
	def __mutateSeqs(mut_rate):
		return None
	@abc.abstractmethod
	def __recombine(parent1, parent2, recomb_rate):
		return None
	@abc.abstractmethod
	def reproduce():
		return None
	@abc.abstractmethod
	def generatePopulation(self, indclass, Ninds):
		return None


'''OrganismTemplate class holds all the information that is common to all the organisms (i.e., species-level information) and that does not vary during a simulation: genome size, number of Transcription Factors and an instance of geneExpression.TFSet which holds PWMs, Kmax, alphasBtm, and other relevant information about this species transcription factors. This class is used to avoid replicating all the data in each instance of Organism'''
class OrganismTemplate(AbstractOrganismTemplate):
	DEFAULT_RECOMBINATION_PERKB = 0.2
	DEFAULT_MUTATION_RATE = 0.005
	def __init__(self, cell_def_mode = None, model_def_file = None, pheno_class = None):
		options = {
			0:self.readModelFile,
			1:self.typeAModel,
			2:self.typeBModel,
			3:self.typeZeroModel,
			4:self.useModelDirectly
		}
		self.cell_def_mode = cell_def_mode
		self.model_file = model_def_file
		self.pheno_class = pheno_class
		model = options[cell_def_mode](model_def_file)
		self.genome_size = model['gsize']
		self.seq_length = model['seqlen']
		self.tf_num = model['tf_num']
		self.tfs = gx.TFSet(range(0, model['tf_num']),direction = model['tf_dir'], kmax = model['tf_kmax'], alphasbtm = model['tf_alphasbtm'], interactions = model['tf_interactions'], difKmax = model['tf_difKmax'], min_inhibitor_index = model['tf_lineage'])
		self.tf_lineage = range(0, model['tf_lineage'])
		self.ncells = model['num_cells']
		self.background = model['background']
		self.mut_rate = model['mut_rate']
		self.rec_rate = model['rec_rate']
		self.environment = self.pheno_class(None, template = self, args = model['selection_args'])
	def readModelFile(self, model_def_file = None):
		return modelReader().readModelFile(model_def_file)
	def getInitialConcentrations(self, i):
		assert(i < self.ncells)
		return self.environment.getInitialConcentrations(i)
	def typeAModel(self, model_def_file = None):
		model = {'gsize':40,
			 'seqlen':150,
			 'tf_num':15,
			 'tf_dir':None,
			 'tf_kmax':None,
			 'tf_alphasbtm':None,
			 'tf_interactions':None, #[gx.interactionTuple(10, 11, 3), gx.interactionTuple(12, 13, 3)],
			 'tf_difKmax':False,
			 'tf_lineage':4,
			 'num_cells':4,
			 'mut_rate':self.DEFAULT_MUTATION_RATE,
			 'rec_rate':self.DEFAULT_RECOMBINATION_PERKB,
			 'background':{'A':0.25,'C':0.25,'G':0.25,'T':0.25},
			 'selection_args':[5, 5, 4, [1, 1, 1, 1], None]
			}
		return model
	def typeZeroModel(self, model_def_file = None):
		model = {'gsize':40,
			 'seqlen':150,
			 'tf_num':15,
			 'tf_dir':None,
			 'tf_kmax':None,
			 'tf_alphasbtm':None,
			 'tf_interactions':None,#[gx.interactionTuple(0, 1, 3), gx.interactionTuple(2, 3, 3)],
			 'tf_difKmax':False,
			 'tf_lineage':4,
			 'num_cells':1,
			 'mut_rate':self.DEFAULT_MUTATION_RATE,
			 'rec_rate':self.DEFAULT_RECOMBINATION_PERKB,
			 'background':{'A':0.25,'C':0.25,'G':0.25,'T':0.25},
			 'selection_args':[15, 0, 1, [1], None] #pan_genes, specific_genes, celltypes, cells_per_type, cellstrings 
			}
		return model
	def typeBModel(self, model_def_file = None):
		model = {'gsize':400,
			 'seqlen':200,
			 'tf_num':40,
			 'tf_dir':None,
			 'tf_kmax':None,
			 'tf_alphasbtm':None,
			 'tf_interactions':None, #[gx.interactionTuple(0, 1, 3), gx.interactionTuple(2, 3, 3),gx.interactionTuple(4, 5, 3),gx.interactionTuple(6, 7, 3),gx.interactionTuple(7, 8, 3)],
			 'tf_difKmax':True,
			 'tf_lineage':8,
			 'num_cells':8,
			 'background':{'A':0.25,'C':0.25,'G':0.25,'T':0.25},
			 'selection_args':[80, 80, 4, [2, 2, 1, 3], None]
			}
		return model
	def useModelDirectly(self, model):
		return model
	def evaluateOrganism(self, organism):
		self.environment.evaluateOrganism(organism)
	def reproduce(self, indcls, parent1, parent2 = None, recomb_rate = None, mut_rate = None, intra_seq_recomb = True):
		if(recomb_rate is None):
			recomb_rate = self.rec_rate
		if(mut_rate is None):
			mut_rate = self.mut_rate
		if(parent2 is None):
			offspringSeqs = self.__mutateSeqs(parent1.seqs, mut_rate)
		else:
			#assert(parent1.template == parent2.template)#for some reason self was different
			offspringSeqs = self.__mutateSeqs(self.__recombine(parent1.seqs, parent2.seqs, recomb_rate, intra_seq_recomb), mut_rate)		
		return indcls(offspringSeqs, self, mut_rate)
	def generatePopulation(self, indclass, Ninds):
		individuals = []
		for i in range(Ninds):
			individuals.append(indclass(self.generateRandomSeqs(self.genome_size, self.seq_length, self.background), self))
		return individuals
	@staticmethod
	def generateRandomSeqs(N, length, bck):
		alphabet = list(bck.keys())
		probs = list(bck.values())
		letters = np.random.choice(alphabet, N*length, replace=True, p = probs).reshape(N, length)
		seqs  = [''.join(r) for r in letters]
		return np.array(seqs)
	def __mutateSeqs(self, seqs, mut_rate = 0.25):
		charlist = [list(s) for s in seqs]

		mutate = np.random.choice(range(0, self.genome_size * self.seq_length), int(mut_rate * self.genome_size * self.seq_length), replace=False)

		rows = np.floor_divide(mutate, self.seq_length)
		cols = np.mod(mutate, self.seq_length)

		alphabet = list(self.background.keys())
		probs = list(self.background.values())

		newchars = np.random.choice(alphabet, len(mutate), replace=True, p = probs)
		for i in range(len(mutate)):
			charlist[rows[i]][cols[i]] = newchars[i]
		return np.array([''.join(r) for r in charlist])

	def __recombine(self, parent1, parent2, recomb_rate, intraseq = True):
		from_p1 = np.random.choice(range(0, self.genome_size), np.floor_divide(self.genome_size, 2), replace = False)
		newseqs = copy.deepcopy(parent2)
		newseqs[from_p1] = copy.deepcopy(parent1[from_p1])
		if(intraseq):
			num_rec_points = np.random.poisson(lam = 0.001*recomb_rate*self.genome_size*self.seq_length, size = 1)
			rec_points = np.random.choice(range(0, self.genome_size * self.seq_length), num_rec_points)
			previous_seq = None
			for i in rec_points:
				num_seq = np.floor_divide(i, self.seq_length)
				position = cols = np.mod(i, self.seq_length)				
				if(position > 0 and position < self.seq_length - 1 and num_seq != previous_seq):
					newseqs[num_seq] = parent1[num_seq][:position] + parent2[num_seq][position:]
					previous_seq = num_seq
		return newseqs

''' Each instance of this class represents an individual. Holds information that is different between individuals and that changes during the simulations: a matrix of size genome_size x number_of_cells, whith the expression levels of each gene, the DNA sequeneces, TFBS annotations and the error (which is calculated by a PhenotypeEvaluator instance). Also holds a reference to the OrganismTemplate class to which the individual belongs. '''
class Organism():	
	def __init__(self, seqs = None, template = None, mut_rate = None):
		self.template = template
		self.expression = np.zeros((self.template.genome_size, self.template.ncells))
		self.seqs = seqs
		self.ann = None
		self.error = 999999
		self.mut_rate = template.mut_rate if mut_rate is None else mut_rate
	def getAnnotation(self):
		self.ann = self.template.tfs.annotate(self.seqs)
	def equilibriumExpression(self, ode_calculator, storeCon = False):
		self.getAnnotation()
		for cell in range(self.template.ncells):
			self.expression[:,cell] = ode_calculator.run(self.ann, self.template.getInitialConcentrations(cell), storeCon = storeCon)
	def timeExpression(self, ode_calculator, time_steps = 10):
		self.getAnnotation()
		self.expression = np.zeros((self.template.genome_size, time_steps, self.template.ncells))
		for cell in range(self.template.ncells):
			self.expression[:,:,cell] = ode_calculator.run(self.ann, self.template.getInitialConcentrations(cell), storeCon = True)
	def mutantExpression(self, ode_calculator, tf_ind):
		self.getAnnotation()
		for g in range(self.template.tf_num, self.template.genome_size):
			self.ann[g].Qonpartial[self.ann[g].sites.getTFs() == tf_ind] = 0
			self.ann[g].Qoffpartial[self.ann[g].sites.getTFs() == tf_ind] = 0
		self.expression = np.zeros((self.template.genome_size, self.template.ncells))
		for cell in range(self.template.ncells):
			self.expression[:,cell] = ode_calculator.run(self.ann, self.template.getInitialConcentrations(cell), storeCon = False)
	def getError(self):
		self.error = self.template.environment.evaluateOrganism(self)
		self.mut_rate = np.power(self.error, 1)*self.template.mut_rate if self.error < 1 else self.error*self.template.mut_rate
		return self.error
	def chromatinData(self, ode_calculator):
		self.getAnnotation()
		self.expression = []
		for cell in range(self.template.ncells):
			self.expression.append(ode_calculator.run(self.ann, self.template.getInitialConcentrations(cell), storeCon = True, fullSeq=True, seqlen = self.template.seq_length))


''' This class performs the simulations. The method getNewGeneration contains the core of the genetic algorithm. With parrun a simulation is run with multiprocessing; with basicrun no multiprocessing is used. Each generation, the organisms to reproduce are picked randomly with a probability that is proportional to 1/error^n, with n being the competitive_power variable (intensity of competition)'''
class multicellEvolver():
	DEFAULT_POPULATION_SIZE=24
	DEFAULT_GENERATIONS = 1500
	DEFAULT_ERROR = 0.01
	DEFAULT_SIM_PARAMS = [0.2, False, 1, 10, 0.001]
	DEFAULT_COMPETITIVE_POWER = 5
	def __init__(self, organismclass, template_org, phenotype_eval, competition, sim_params = None, pop_size = None, generations = None, error = None, modelFile = None, sexual = False, predefined_model = 1, instance_name = 'mce0'):
		if(modelFile is not None and predefined_model != 0):
			self.model = modelReader().readModelFile(modelFile) 
			self.template = template_org(4, self.model, phenotype_eval)
		elif(predefined_model == 0):
			self.model = None
			self.template = template_org(0, modelFile, phenotype_eval)
		else:
			self.model = None
			self.template = template_org(predefined_model, None, phenotype_eval)
		self.pop_size = pop_size if pop_size is not None else self.DEFAULT_POPULATION_SIZE
		self.generations = generations if generations is not None else self.DEFAULT_GENERATIONS
		self.error = error if error is not None else self.DEFAULT_ERROR
		if (self.model is None):
			self.organismclass = organismclass
		elif('ORGANISM_CLASS' not in self.model.keys()):
			self.organismclass = organismclass
		else:
			self.organismclass = eval(self.model['ORGANISM_CLASS'])
		self.population = self.template.generatePopulation(self.organismclass, self.pop_size)
		self.competitive_power = competition if competition is not None else self.DEFAULT_COMPETITIVE_POWER
		if(self.model is None):
			sim_params = self.DEFAULT_SIM_PARAMS if sim_params is None else sim_params
			self.gxcalc_class = gx.ODERunner
			self.gxcalc = gx.ODERunner(self.template.tfs, betas = sim_params[0] , with_interactions = sim_params[1], h = sim_params[2], tmax = sim_params[3], min_variation = sim_params[4])
		else:
			sim_params = self.model['sim_params']
			self.gxcalc_class = eval('gx.' + self.model['ODE_RUNNER_CLASS'])
			self.gxcalc = self.gxcalc_class(self.template.tfs, *sim_params)
			
		self.sexual = sexual
		self.last_generation = 0
		self.error_means_acum = []
		self.instance_name = instance_name
		self.param_history = [{'generation':0,'mut_rate':self.template.mut_rate,'competitive_power':self.competitive_power}]
	def parrun(self, save = True, saving_freq = 100):
		import multiprocessing
		competitivity = np.zeros(self.pop_size)
		mean_comp = 0
		if(save):
			self.saveSelf()
		for g in range(self.last_generation, self.generations):
			#Heuristics to edit parameters
			self.paramHeuristics(g, mean_comp, np.var(competitivity))
			#True algorithm: calculate expression and error with multiprocessing
			current_pop = multiprocessing.Queue()
			jobs = [multiprocessing.Process(target = self.parstep, args = (i, current_pop)) for i in range(len(self.population))]
			for j in jobs:
				j.start()
			self.population = [current_pop.get() for j in jobs]
			for i in range(len(self.population)):
				competitivity[i] = self.population[i].error				
			mean_comp = np.mean(competitivity)
			self.error_means_acum.append(mean_comp)
			#Save and print
			if(save and (g+1)%saving_freq == 0):
				self.saveSelf()
				print('generation ', g, ': mean = ', mean_comp, ', best = ', min(competitivity),', mut_rate = ',self.template.mut_rate, ', k = ',self.k)
			if(mean_comp <= self.error):
				break
			#Get new generation
			self.population = self.getNewGeneration(competitivity)
			self.last_generation = g
		return (g, mean_comp)		
	def parstep(self, i, current_pop):
		o = self.population[i]
		o.equilibriumExpression(self.gxcalc)
		err = o.getError()
		current_pop.put(o)
		#print('org ', i, ' ', err, '\t')
	def parstepTime(self, i, current_pop):
		o = self.population[i]
		o.timeExpression(self.gxcalc, self.sim_params[3])
		err = o.getError()
		current_pop.put(o)
	def basicrun(self, save = True, saving_freq = 100):
		competitivity = np.zeros(self.pop_size)
		mean_comp = 0
		for g in range(self.generations):
			if(save and (g+1)%saving_freq == 0):
				self.saveSelf()
			self.paramHeuristics(g, mean_comp, np.var(competitivity))
			i = 0
			for o in self.population:
				o.equilibriumExpression(self.gxcalc)
				competitivity[i] = o.getError()
				print('org ', i, ' ', competitivity[i], '\t')
				i+=1
			mean_comp = np.mean(competitivity)
			print('generation ', g, ': ', mean_comp)
			if(mean_comp <= self.error):
				break
			self.population = self.getNewGeneration(competitivity)
		return (g, mean_comp)
	def getNewGeneration(self, error):
		efficacy = 1/np.power(error, self.competitive_power)
		efficacy = efficacy/np.sum(efficacy)
		new_gen = []
		parents = []
		for i in range(self.pop_size):
			if(self.sexual):
				parents = np.random.choice(self.population, 2, replace = True, p = efficacy)
				mut_rate = 0.5*(parents[0].mut_rate + parents[1].mut_rate)
				new_gen.append(self.template.reproduce(self.organismclass, parents[0], parents[1], mut_rate = mut_rate))
			else:
				parent = np.random.choice(self.population, 1, p = efficacy)[0]
				new_gen.append(self.template.reproduce(self.organismclass, parent, None, mut_rate = parent.mut_rate))
		return new_gen
	def saveSelf(self):
		filename = self.instance_name + 'g' + str(self.last_generation) + '.pck'
		with open(filename, 'wb') as output: 
 			pickler = pck.Pickler(output, -1)
 			pickler.dump(self)
	def setComp(self, cp):
		self.competitive_power = cp
		self.param_history.append({'generation':self.last_generation,'mut_rate':self.template.mut_rate,'k':self.k})
	def setMutRate(self, mut_rate):
		self.mut_rate = mut_rate
		self.param_history.append({'generation':self.last_generation,'mut_rate':self.template.mut_rate,'k':self.competitive_power})
	def paramHeuristics(self, g, mean_comp, var_comp):
		if(g == 100):
			self.setMutRate(0.001)
		if(g > 300 and (g+1)%20 == 0):
			self.setMutRate(1.5/(self.template.genome_size*self.template.seq_length*mean_comp))
		if(g > 500 and (g+1)%20 == 0):
			self.setMutRate(1.25/(self.template.genome_size*self.template.seq_length*mean_comp))
		if(var_comp == 0 and mean_comp > 0 and g > 10):
			self.setMutRate(1.5/(self.template.genome_size*self.template.seq_length*mean_comp))
	def printOrg(self, i = 0):
		o = self.population[i]
		o.equilibriumExpression(self.gxcalc)
		er = o.getError()
		print("error of org ", i, ": ", er)
		print("SITES/EXPRESSION:")
		for i in range(self.template.genome_size):
			if(i >= self.template.tf_num):
				if((i - self.template.tf_num)%5 == 0):
					print('__________\n')
			if(o.ann[i] is not None):
				print(i, ': ', o.ann[i].sites.getTFs(), '||||', o.expression[i, :])
			else:
				print(i, ': None', '||||', o.expression[i, :])
		print("INHIBITORS: ", np.where(self.template.tfs.alphasbtm<1))
		print("INTERACTIONS: ", self.template.tfs.interactions)
		print("MEAN BY TYPE:")
		for j in [range(15,20), range(20,25), range(25,30), range(30,35), range(35,40)]:
			print([np.mean(o.expression[j, i]) for i in range(4)])
		print('__________\n')
	@staticmethod
	def readMCE(filenamein, iname="default_name"):
		with open(filenamein, 'rb') as f:
			print('opening ' + filenamein)
			import _pickle as pck
			mce = pck.load(f)
			mce.last_generation += 1
			print('template mut_rate: ', mce.template.mut_rate, ', k:', mce.k, ',tf type: ', mce.template.tfs.direction, '\n')
			if(iname != mce.instance_name):
				mce.instance_name = mce.instance_name + '_' + iname
		return mce

	def produceDataFrames(self, options, simname):
		import multiprocessing
		global pd
		import pandas as pd
		self.fix()	## Kmax of TF were not being multiplied. In order to be able to shut them off, they are now multiplied but all should be equal to one 
		competitivity = np.zeros(self.pop_size)
		for i in range(len(competitivity)):
			competitivity[i] = self.population[i].error
		if(any(competitivity > 10000)):
			print('recalculating all expression... ', simname)
			current_pop = multiprocessing.Queue()
			jobs = [multiprocessing.Process(target = self.parstep, args = (i, current_pop)) for i in range(len(self.population))]
			for j in jobs:
				j.start()
			self.population = [current_pop.get() for j in jobs]
			for i in range(len(self.population)):
				competitivity[i] = self.population[i].error	
		winner_ind = np.where(competitivity == np.min(competitivity))[0][0]			
		winner = self.population[winner_ind]
		error = winner.getError()
		outname = simname.replace('.pck', '_o'+str(winner_ind)+'err'+str(round(error, 4)))
		if(options.find('b') >= 0):
			d = self.finalExpressionToDF(winner, outname)
			d = self.annotationToDF(winner, outname)
			d = self.tfsToDF(outname)
		if(options.find('t') >= 0):
			d = self.timeExpressionToDF(winner, outname)
		if(options.find('m') >= 0):
			d = self.mutantAnalysisToDF(winner, outname)
		if(options.find('s') >= 0):
			d = self.mutantSitesByGroupToDF(winner, outname)
	def finalExpressionToDF(self, o, outname = ''):	
		df = pd.DataFrame(o.expression, columns = ['exp' + str(i) for i in range(o.template.ncells)])
		for i in range(o.template.ncells):
			df['opt' + str(i)] = o.template.environment.cells[i].optimal_concentrations
		for i in range(o.template.ncells):
			df['init' + str(i)] = o.template.environment.cells[i].initial_concentrations
		cols = [col for col in df.columns if 'opt' in col]
		df['type'] = df[cols].apply(lambda row: np.array2string(np.where(row>0)[0], separator=''), axis=1)
		df['type'][0:self.template.tf_num] = [1 if i == self.template.tfs.ACTIVATOR_TYPE else -1 for i in self.template.tfs.direction]
		df['seqs'] = o.seqs
		if(outname != ''):
			df.to_csv(outname +'_finalExpression.csv', sep='\t', header=True,decimal='.', float_format='%.10f')
		return df
	def annotationToDF(self, o, outname = ''):
		df = pd.DataFrame()
		for i in range(o.template.genome_size):
			if(o.ann[i] is not None):
				aux = pd.DataFrame(o.ann[i].sites.get())
				aux['gene'] = o.ann[i].ind
				aux['QonPartial'] = o.ann[i].Qonpartial
				aux['QoffPartial'] = o.ann[i].Qoffpartial
				df = df.append(aux)
		if(outname != ''):
			df.to_csv(outname +'_annotation.csv', sep='\t', header=True,decimal='.', float_format='%.10f')
		return df
	def tfsToDF(self, outname):
		df = pd.DataFrame()
		m = self.template.tfs.getPSSMs()
		for i in range(self.template.tf_num):
			aux = pd.DataFrame()
			for l in m[i].pos.alphabet.letters:
				aux[l] = m[i].pos.get(l)
			aux['tf_name'] = self.template.tfs.getSourceIndex(i)
			aux['alpha_btm'] = self.template.tfs.alphasbtm[self.template.tfs.getSourceIndex(i)]
			aux['k_max'] = self.template.tfs.kmax[self.template.tfs.getSourceIndex(i)]
			aux['type'] = self.template.tfs.direction[self.template.tfs.getSourceIndex(i)]
			aux['consensus'] = str(m[i].pos.consensus)
			df = df.append(aux)
		if(outname != ''):
			df.to_csv(outname +'_TFset.csv', sep='\t', header=True,decimal='.', float_format='%.10f')
		df2 = self.interactionsToDF(outname)
		return df
	def interactionsToDF(self, outname):
		df = pd.DataFrame()
		interactions = self.template.tfs.interactions
		tf1 = []
		tf2 = []
		weight = []
		dir_tf1 = []
		dir_tf2 = []
		if(len(interactions) > 0):
			for i in interactions:
				tf1.append(i.tf1)
				tf2.append(i.tf2)
				weight.append(i.weight)
				dir_tf1.append(self.template.tfs.direction[i.tf1])
				dir_tf2.append(self.template.tfs.direction[i.tf2])
		df['tf1'] = tf1
		df['tf2'] = tf2
		df['weight'] = weight
		df['tf1_type'] = dir_tf1
		df['tf2_type'] = dir_tf2
		if(outname != ''):
			df.to_csv(outname +'_TFinteractions.csv', sep='\t', header=True,decimal='.', float_format='%.10f')
		return df

		
	def timeExpressionToDF(self, o, outname):
		expr = o.expression
		time = int(self.gxcalc.tmax/self.gxcalc.h)
		o.timeExpression(self.gxcalc, time)
		df = pd.DataFrame()
		for i in range(time):
			aux = pd.DataFrame(o.expression[:,i,:], columns = ['cell' + str(n) for n in range(o.template.ncells)])
			aux['time'] = i*self.gxcalc.h
			aux['gene'] = range(0,self.template.genome_size)
			df = df.append(aux)
		if(outname != ''):
			df.to_csv(outname +'_timeExpression.csv', sep='\t', header=True,decimal='.', float_format='%.10f')
		o.expression = expr
		return df
	def mutantCalculate(self, o, i, mutant_exp):
		true_kmax = o.template.tfs.kmax[i]
		o.template.tfs.kmax[i] = 0
		o.equilibriumExpression(self.gxcalc)
		err = o.getError()
		aux = self.finalExpressionToDF(o, '')
		aux['tf_mutated'] = i
		aux['overall_error'] = err
		o.template.tfs.kmax[i] = true_kmax
		mutant_exp.put(aux)
	def mutantAnalysisToDF(self, o, outname):
		print('mutant analysis beginning... ', outname)
		df = pd.DataFrame()
		import multiprocessing
		#for i in range(self.template.tf_num):
		mutant_exp = multiprocessing.Queue()
		jobs = [multiprocessing.Process(target = self.mutantCalculate, args = (o, i, mutant_exp)) for i in range(self.template.tf_num)]
		for j in jobs:
			j.start()
		for j in jobs:
			aux = mutant_exp.get()
			df = df.append(aux)
		if(outname != ''):
			df.to_csv(outname +'_mutantTFs.csv', sep='\t', header=True,decimal='.', float_format='%.10f')
		return df
	def mutSitesCalculate(self, o, i, mutant_exp):
		o.mutantExpression(self.gxcalc, i)
		aux = self.finalExpressionToDF(o, '')
		aux['mutated_sites'] = i
		mutant_exp.put(aux)			
	def mutantSitesByGroupToDF(self, o, outname):
		print('mutating sites... ', outname)
		df = pd.DataFrame()
		import multiprocessing
		mutant_exp = multiprocessing.Queue()
		jobs = [multiprocessing.Process(target = self.mutSitesCalculate, args = (o, i, mutant_exp)) for i in range(self.template.tf_num)]
		for j in jobs:
			j.start()
		for j in jobs:
			aux = mutant_exp.get()
			df = df.append(aux)	
		if(outname != ''):
			df.to_csv(outname +'_mutantSites.csv', sep='\t', header=True,decimal='.', float_format='%.10f')
		return df
	def fix(self):
		for i in range(self.template.tf_num):
			self.template.tfs.kmax[i] = self.template.tfs.DEF_KMAX
		self.template.tfs.kmax = np.array(self.template.tfs.kmax)
		for p in self.population:
			p.template.tfs.kmax = self.template.tfs.kmax
''' This class performs simulations just like multicellEvolver. The only difference is that it uses tournament algorithm to produce the next generation. It seems to work better'''
class multicellTournament(multicellEvolver):
	DEFAULT_K = 2
	def __init__(self, organismclass, template_org, phenotype_eval, competition, sim_params = None, pop_size = None, generations = None, error = None, modelFile = None, sexual = False, predefined_model = 1, instance_name = 'mce0', k = None):
		super().__init__( organismclass, template_org, phenotype_eval, competition, sim_params, pop_size, generations, error, modelFile, sexual, predefined_model, instance_name)
		self.k = self.DEFAULT_K if k is None else k
		self.changedK = False
		self.param_history = [{'generation':0,'mut_rate':self.template.mut_rate,'k':self.k}]
	def getNewGeneration(self, error):
		efficacy = 1/error
		efficacy = efficacy/np.sum(efficacy)
		new_gen = []
		parents = []
		for i in range(self.pop_size):
			if(self.sexual):
				parents = [self.population[self.__tournament(efficacy)], self.population[self.__tournament(efficacy)]]
				mut_rate = 0.5*(parents[0].mut_rate + parents[1].mut_rate)
				new_gen.append(self.template.reproduce(self.organismclass, parents[0], parents[1], mut_rate = mut_rate))
			else:
				parent = self.population[self.__tournament(efficacy)]
				new_gen.append(self.template.reproduce(self.organismclass, parent, None, mut_rate = parent.mut_rate))
		return new_gen
	def __tournament(self, efficacy):
		participants = np.random.choice(range(len(efficacy)), self.k, replace = False)
		m = np.argmax(efficacy[participants])
		return participants[m]
	def setK(self, k):
		self.k = k
		self.param_history.append({'generation':self.last_generation,'mut_rate':self.template.mut_rate,'k':self.k})
	def setMutRate(self, mut_rate):
		self.template.mut_rate = mut_rate
		self.param_history.append({'generation':self.last_generation,'mut_rate':self.template.mut_rate,'k':self.k})
	def paramHeuristics(self, g, mean_comp, var_comp):
		if(mean_comp < 0.0001):
			return
		if(g == 100):
			self.setMutRate(0.001)
		if(g > 300 and (g+1)%20 == 0):
			self.setMutRate(2/(self.template.genome_size*self.template.seq_length*mean_comp))
		if(g > 500 and (g+1)%20 == 0):
			self.setMutRate(2/(self.template.genome_size*self.template.seq_length*mean_comp))
		if(var_comp == 0 and mean_comp > 0 and g > 10):
			self.setMutRate(2/(self.template.genome_size*self.template.seq_length*mean_comp))
		if(g > 300 and mean_comp < 0.5 and not self.changedK):
			self.setK(min(int(self.k + 0.5*self.k), self.pop_size - 1))
			self.changedK = True






### Syntax example:
### Simulation with default name
#	 python cellEvolver.py 
### Simulation with name, 4cell_mce1.cemod model:
# 	 python cellEvolver.py SimulationName
### Simulation with name and model file as input:
# 	 python cellEvolver.py SimulationName modelX.cemod
### Pass a pickle file and make a basic print of it:
# 	python cellEvolver.py newName -filename.pck
### Pass a pickle file and produce data frames (b for basic graph data, t for developmental time analysis. m for mutant TFs and s for TF sites in terminal features mutation):
# 	python cellEvolver.py newName -filename.pck btms
### Pass a pickle and execute another argument
#	python cellEvolver.py newname -filename.pck eval 'import XX;print(xy)'
### Pass a pickle and keep evolving:
# 	python cellEvolver.py newName filename.pck


def main():
	iname = sys.argv[1]
	if(iname is not None):
		print('Simulation name: ', iname, '\n')
	else:
		iname = 'mceX'
		print('Simulation name (default name): ', iname, '\n')
	if(len(sys.argv)>2):
		filenamein  = sys.argv[2]
		if(filenamein.endswith('.pck')):
			if(filenamein[0] == '-'):
				filenamein = filenamein[1:]
				mce = multicellEvolver.readMCE(filenamein, iname)
				if(len(sys.argv)>3):
					if(sys.argv[3] == 'eval'):
						for expression in sys.argv[4].split(';'):
							exec (expression)
					else:
						mce.produceDataFrames(sys.argv[3], filenamein)
				else:
					mce.printOrg()
				return mce			
			else:
				mce = multicellEvolver.readMCE(filenamein, iname)
		elif(filenamein.endswith('.cemod')):
			mce = multicellTournament(organismclass = Organism, template_org = OrganismTemplate, phenotype_eval = SimplePhenotypeEvaluator, competition = 6, sim_params = None, pop_size = 24, generations = 10000, error = 0.01, modelFile = filenamein, sexual = True, predefined_model = 4,  instance_name = iname, k = 12)
	else:
		mce = multicellTournament(organismclass = Organism, template_org = OrganismTemplate, phenotype_eval = SimplePhenotypeEvaluator, competition = 6, sim_params = None, pop_size = 24, generations = 10000, error = 0.01, modelFile = '4cell_mce1.cemod', sexual = True, predefined_model = 4,  instance_name = iname, k = 12)
	print(mce.instance_name, ': Initiating simulations...\n')
	mce.parrun(True, 200)
	mce.saveSelf()
	mce.produceDataFrames('btms', iname)
	return mce

if __name__ == "__main__":
    main()

