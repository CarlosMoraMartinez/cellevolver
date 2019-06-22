from collections import namedtuple
from Bio import motifs
from Bio.Seq import Seq
import numpy as np
import time

PSSMTuple = namedtuple("PSSMTuple", "pos neg")
interactionTuple = namedtuple("interactionTuple", "tf1 tf2 weight")

a=1+1
class MotifGenerator():
	""" produces random PWMs from dirichlet distributions with different alphas for each nucleotide. No instance of this class should be generated (lacks constructor);
		most of the variables are constants and are assumed to be parameters of the model. Uses the following algorithm to generate PSSMs:
		
		1) For i in MOTIF_SIZE, produces a probability distribution from a Dirichlet probability distribution with a given alpha. A diferent alpha is specified for each i, so that alphas in the 			   middle are smaller and yield more unbalanced distributions.
		2) Basing on this distributions, randomly generate a set of nucleotides (of size == NSITES) 
		3) Join the nucleotides across the i sites, in order to get a sample of NSITES 'functional' sites. Each of the sites has a particular nucleotide distribution
		4) Use the motifs library to generate pssms from your motif sample, using PSEUDOCOUNTS and DEFAULT_BACKGROUND

	METHODS:
		getPSSM: Optionally gets:
					-alphas: a list of alphas (a Dirichlet alpha for each nucleotide in the motif)
					-background: a dictionary with keys 'A', 'C', 'G', 'T'. Values should sum up to 1. Defaults all to 0.25
					@returns a motifs.pssm object randomly generated.

 """
	NSITES = 100	#number of random sites generated to calculate the PWM
	PSEUDOCOUNTS = 0.1*NSITES
	ALPHABET = ('A', 'C', 'G', 'T')
	MOTIF_SIZE = 8
	DEFAULT_ALPHAS = [0.4, 0.3, 0.2, 0.1, 0.1, 0.2, 0.3, 0.4]	#length must be equal to MOTIF_SIZE
	DEFAULT_BACKGROUND = {'A':0.25,'C':0.25,'G':0.25,'T':0.25}
	@classmethod
	def getPSSM(self, alphas = None, background = None):
		if(alphas == None):
			alphas = self.DEFAULT_ALPHAS
		if(background == None):
			background = self.DEFAULT_BACKGROUND
		probs = np.zeros([len(self.ALPHABET), len(alphas)]) ## Can't initialize a Motif object directly with this matrix
		counts = np.chararray([self.NSITES, len(alphas)])	## Therefore I should generate random sites
		instances = []
		for i in range(len(alphas)):
			probs[:, i] = np.random.dirichlet(np.full(len(self.ALPHABET), alphas[i]), size = 1)
			counts[:, i] = np.random.choice(self.ALPHABET, size = self.NSITES, replace = True, p = probs[:, i])
		for i in range(self.NSITES):
			instances.append(counts[i, :].tostring().decode())
		m = motifs.create(instances)
		#m.weblogo("logo.png")
		pssm = m.counts.normalize(pseudocounts=self.PSEUDOCOUNTS).log_odds(background = background)
		res = PSSMTuple(pssm, pssm.reverse_complement())
		return res


class SeqAnnotation():
	""" Holds the motifs found in a sequence, but not other features such as interaction matrix between sites etc. AnnotatedSeq objects have one of these in their sites field.
		 Also, several functions to easily get the TF binding sites data.
	FIELDS:
		-formats (constant): formats of the self.arr np.array columns
		-fields (constant): names of the self.arr np.array columns
		-arr: numpy array with defined names and formats containing the TF BINDING SITES with their score, position, strand, etc
	METHODS:
		-assign: Takes field, value. Assigns an array or list to the specified field. Fields can be specified with strings or numbers. value must be of shape = (n, 1), where n equals self.arr rows
		-append: Concatenates another SeqAnnotation object to self.
		-getTFs 
		-getLLR
		-getPercents
		-getStrands
		-getEfSite
		-getField: takes the name or number of a field of self.arr. @returns a vector with that field
		-length: returns the number of motifs
		-get: returns the full numpy array with all the fields.

	
	The core data structure is a numpy array with names and types defined for each column. Fields and types of this array are defined in the formats and fields constants.

	"""
	

	formats = ['int8', 'int16', 'float32', 'float32', 'int8', 'float32']
	fields = ['tf_ind', 'positions', 'LLR', 'percent', 'strand', 'efsite']
	def __init__(self, size = 0):
		self.arr = np.zeros(size, dtype = {'names':self.fields, 'formats':self.formats})
	def assign(self, field, value): 
		if(isinstance(field, str)):
			self.arr[field] = value
		elif(isinstance(field, int)):
			self.arr[self.fields[field]] = value
	def append(self, other):
		self.arr = np.append(self.arr, other.get())
	def getTFs(self):
		return self.arr['tf_ind']
	def getPos(self):
		return self.arr['positions']
	def getLLR(self):
		return self.arr['LLR']
	def getPercents(self):
		return self.arr['percent']
	def getStrands(self):
		return self.arr['strand']
	def getEfSite(self):
		return self.arr['efsite']
	def getField(self, field):
		if(isinstance(field, str)):
			return self.arr[field]
		elif(isinstance(field, int)):
			return self.arr[self.fields[field]]
	def length(self):
		return len(self.arr)
	def get(self):
		return self.arr

		

class AnnotatedSeq():
	"""Holds a SeqAnnotation object (the TF binding sites) and their interactions and overlaps"""
	MAX_INTERACTION_RANGE = 50
	def __init__(self, ind, sites, interactions = None, Qonpartial = None, Qoffpartial = None):
		self.ind = ind
		self.sites = sites
		self.nsites = self.sites.length()
		self.intmatrix = self.getInteractionMatrix(interactions)
		self.Qonpartial = Qonpartial
		self.Qoffpartial = Qoffpartial
	def getInteractionMatrix(self, interactions):
		if(interactions is None):
			mat = np.ones((self.nsites, self.nsites))
			dists = self.getDistances(self.sites.getPos())
			mat[dists < MotifGenerator.MOTIF_SIZE] = 0
			return mat
		elif(self.nsites is 0 or self.nsites is None):
			return None
		elif(self.nsites is 1):
			return np.ones((1, 1))
		mat = np.ones((self.nsites, self.nsites))
		dists = self.getDistances(self.sites.getPos())
		for pair in interactions:
			ind1 = np.where(self.sites.getTFs() == pair.tf1)
			ind2 = np.where(self.sites.getTFs() == pair.tf2)
			if(len(ind1[0])>0 and len(ind2[0])>0):
				for i in np.nditer(ind1):
					for j in np.nditer(ind2):
						if(dists[i, j] <= self.MAX_INTERACTION_RANGE):
							mat[i, j] = pair.weight
							mat[j, i] = pair.weight
		"""Motifs that are too close compete (cant bind at the same time)"""
		mat[dists < MotifGenerator.MOTIF_SIZE] = 0
		return mat
	@classmethod
	def getDistances(self, pos):
		d = np.ones((len(pos), len(pos)))
		for i in np.nditer(range(len(pos))):
			d[:,i] = np.absolute([p - pos[i] for p in pos])
		return d

class ExpressionCalculator():
	""" Calculates expression for given sequence annotations """
	BASAL_EXPRESSION = 0.0
	def __init__(self):
		self.current_seq  = None
		self.stack = None		
	def activate(self, seq, conTF = None, interactions = False, siteOpening = None):
		if(seq is None):
			return self.BASAL_EXPRESSION
		else:
			self.current_seq = seq

		if(conTF is None):
			concentrations = np.ones(len(seq.Qonpartial))
		else:
			concentrations = np.array([conTF[i] for i in seq.sites.getTFs()])
		if(siteOpening is not None):
			concentrations = concentrations*siteOpening
		if(interactions is False):
			""" Qon is exp(LLRi - LLRmax)*Kmax*alphasbtm. exp is calculated by SeqAnnotator"""
			""" Qoff is exp(LLRi - LLRmax)*Kmax"""
			indexes = np.where(concentrations > 0)  #### this has not been tested!!!!
			#qonvec = np.array(seq.Qonpartial*concentrations)
			#qoffvec = np.array(seq.Qoffpartial*concentrations)
			qonvec = seq.Qonpartial*concentrations
			qoffvec = seq.Qoffpartial*concentrations
			Zon = self.partitionGeneral(qonvec[indexes]) + self.BASAL_EXPRESSION
			Zoff = self.partitionGeneral(qoffvec[indexes]) + 1
		else:
			#self.current_seq = seq
			self.stack = []
			Zon = self.partitionInteraction(seq.Qonpartial*concentrations) + self.BASAL_EXPRESSION
			Zoff = self.partitionInteraction(seq.Qoffpartial*concentrations) + 1
		
		E = Zon/(Zon + Zoff)
		return E

	@classmethod
	def partitionGeneral(self, arr, i = 0):
		if(len(arr) == 1):
			return arr[0]
		elif(len(arr) == 0 or len(arr) > 40):
			return 0
		m = np.zeros(np.power(len(arr), 2)).reshape([len(arr), len(arr)])
		m[0, :] = arr
		for i in range(1, len(arr)):
			for j in range(i, len(arr)):
				m[i, j] = m[0, j]*np.sum(m[i-1, :j])
		return np.sum(m)

	def partitionInteraction(self, arr, i = 0):
		if(len(arr) > 25):
			#raise ValueError('partitionOff overflow (more than 50 motifs)\n')
			return 0
		elif(i >= len(arr)):
			return 0
		elif(i == (len(arr) - 1)):
			self.stack.append(i)
			inters  = self.getInteractions()
			self.stack.pop()
			return arr[i-1]*inters

		r = 0
		for j in range(i, len(arr)):
			self.stack.append(j)
			inters  = self.getInteractions()
			r += (arr[j]*inters + inters*arr[j]*self.partitionInteraction(arr, j+1))
			self.stack.pop()	
		return r
	
	def getInteractions(self):
		if(len(self.stack)<2):
			return 1
		else:
			last  = self.stack[len(self.stack)-1]
			others = self.stack[0:(len(self.stack)-1)]
			return np.prod(self.current_seq.intmatrix[last, others])
		

class ODERunner():
	""" Uses ExpressionCalculator class to run simulations over a period of time. Uses Euler integration to solve the following system:
		Concentration(t+1) = Concentration(t) + activate(annotation, TFconcentrations(t)) - Beta*Concentration(t)
		Uses a vectorized version of ExpressionCalculator.activate to get all the activation values at the same time.
	FIELDS:
		- tfs: a TFSet object. Usually initialized before passing it to the constructor
		- with_interactions: boolean variable. If True, will use ExpressionCalculator.activate() with interactions set to True, i.e., will use partitionInteraction()
		- h: Step size for numeric integration (Euler method). With 0.1 is not very slow, so used as default.
		- tmax: Will stop after tmax/h iterations (or when error < min_variation)
		- min_variation: if the maximum variation (among all sequences) in concentration in a step is smaller than min_variation*h, the loop will stop
		- betas: degradation coefficients for each gene. Either atomic or of len == num sequences
		- calculator: ExpressionCalculator object
		- vcalc: vectorized version of calculator.activate method


"""
	def __init__(self, tfs, betas, with_interactions = True, h = 0.1, tmax = 100, min_variation = 0.00001):
		self.tfs = tfs
		self.with_interactions = with_interactions
		self.h = h
		self.tmax = tmax
		self.min_variation = min_variation 
		self.betas = betas
		self.calculator =  ExpressionCalculator()
		self.vcalc = np.vectorize(self.calculator.activate, excluded = ['conTF', 'interactions', 'siteOpening'])

	def step(self, concentration, annotation, conTF):
		der = self.vcalc(annotation, conTF = conTF, interactions = self.with_interactions, siteOpening = None) - self.betas*concentration
		y_next = concentration + self.h*der
		y_next[y_next < 0] = 0.0
		return y_next
	def run(self, annotation, start, storeCon = True):
		max_time = int(self.tmax/self.h)
		if(storeCon):
			concentrations = np.zeros((len(start), max_time))
			concentrations[:, 0] = start
			for t in range(1, max_time):
				concentrations[:, t] = self.step(concentration = concentrations[:, t-1], annotation = annotation, conTF = concentrations[self.tfs.source_inds, t-1])
				if(np.max(np.abs(concentrations[:, t] - concentrations[:, t-1])) <= self.min_variation*self.h):
					concentrations = concentrations[:, range(0, t+1)]
					break
		else:
			y_last = start
			for t in range(1, max_time):
				concentrations = self.step(concentration = y_last, annotation = annotation, conTF = y_last[self.tfs.source_inds])
				if(np.max(np.abs(y_last - concentrations)) <= self.min_variation):
					break
				y_last = concentrations
		return concentrations


ChromatinResult = namedtuple("ChromatinResult", "concentrations chromatinState")
class ODERunnerChromatin(ODERunner):
	ACTIVATOR_STRENGTH = 1
	INHIBITOR_STRENGTH = -3
	def __init__(self, tfs, betas, with_interactions = True, h = 0.1, tmax = 100, min_variation = 0.00001, opening_start = 0.5, opening_basal_rate = -0.2, max_opening = 5, scale = 0.2, normalize = (True, False, False, False), amplitude = 30, activator_strength = 1, inhibitor_strength = -3, tf_strengths=None, min_opening = 0.1):
		super().__init__(tfs, betas, with_interactions, h, tmax, min_variation)
		self.opening_start = opening_start
		self.opening_basal_rate = opening_basal_rate
		self.max_opening = max_opening
		self.min_opening = min_opening
		self.scale = scale
		self.normalize_by_TFnum, self.normalize_by_influence, self.normalize_by_conc, self.normalize_total = normalize
		if(isinstance(amplitude, (int, float))):
			self.all_same_amplitude = True
			self.amplitude = np.array([amplitude for i in tfs.inds])
		else:
			self.amplitude = np.array(amplitude)
			self.all_same_amplitude = False
		self.open_precalc = []
		self.vcalc = np.vectorize(self.calculator.activate, excluded = ['conTF', 'interactions']) #siteOpening is vectorized in this class
		#self.gauss = lambda x, med, sd, alpha: alpha*np.exp(-np.power(x-med, 2)/(2*np.power(sd, 2))) #can't pickle lambda functions!
		self.site_size = len(tfs.getPSSMs()[0].pos['G'])
		if(tf_strengths is None):
			if(activator_strength is not None):
				self.activator_strength = activator_strength
			else:
				self.activator_strength = self.ACTIVATOR_STRENGTH
			if(inhibitor_strength is not None):
				self.inhibitor_strength = inhibitor_strength
			else:
				self.inhibitor_strength = self.INHIBITOR_STRENGTH
			self.chromatin_influence_strength = np.array([self.activator_strength if i==TFSet.ACTIVATOR_TYPE else self.inhibitor_strength for i in self.tfs.direction ])
		else:
			assert len(tf_strengths) == len(self.tfs.source_inds), "ODERunnerChromatin initialization Error: different number of TFs and tf_strengths"
			self.chromatin_influence_strength = np.array(tf_strengths)
	@staticmethod
	def gauss(x, med, sd, alpha):
		return alpha*np.exp(-np.power(x-med, 2)/(2*np.power(sd, 2)))	
	def precalculateChromatinMatrixes(self, annotation, fullSeq=False, seqlen=300):
		self.open_precalc = []
		if(not fullSeq):
			for ann in annotation:
				if(ann is None):
					mat = np.array([0])
				else:
					pos = ann.sites.getPos() + (self.site_size-1)/2
					mat = np.concatenate([self.gauss(pos, pos[i], self.amplitude[t], self.chromatin_influence_strength[t]).reshape([pos.shape[0], 1]) for i, t in enumerate(ann.sites.getTFs())], axis = 1)
				self.open_precalc.append(mat)
		else:
			for ann in annotation:
				if(ann is None):
					mat = np.array(range(max(seqlen, np.max(ann.sites.getPos()) + self.site_size)))
					mat = mat.reshape([mat.shape[0], 1])
				else:
					pos = ann.sites.getPos() + (self.site_size-1)/2	#center of the motif
					nts = np.array(range(max(seqlen, np.max(ann.sites.getPos()) + self.site_size)))
					mat = np.concatenate([self.gauss(nts, pos[i], self.amplitude[t], self.chromatin_influence_strength[t]).reshape([nts.shape[0], 1]) for i, t in enumerate(ann.sites.getTFs())], axis = 1)
				self.open_precalc.append(mat)
				
	def changeChromatin(self, seq, mat, conTF = None):
		if(seq is None):
			return np.array([0])
		if(conTF is None):
			concentrations = np.ones(len(seq.Qonpartial))
		else:
			concentrations = np.array([conTF[i] for i in seq.sites.getTFs()])
		#directions = self.chromatin_influence_strength[seq.sites.getTFs()]
		influence = (concentrations*seq.Qoffpartial).reshape(1, mat.shape[1])
		influence = np.sum(mat*influence, axis=1)
		if(self.normalize_by_TFnum):
			influence = influence/mat.shape[1]
		if(self.normalize_by_influence):
			influence = influence/np.sum(np.absolute(influence))
		if(self.normalize_by_conc):
			influence = influence/np.sum(concentrations)
		return influence

	def step(self, concentration, annotation, conTF, chromatinState, fullSeq=False):
		#1) calculate first change in gene concentrations
		if(fullSeq):		## necessary to compute activation (filter out nucleotides without motif center)
			chromatin_state_all=chromatinState
			chromatinState = [chrst[ann.sites.getPos() + int(self.site_size/2)] if ann is not None else np.array([0]) for chrst, ann in zip(chromatinState, annotation)]
		der = np.array([self.calculator.activate(seq=ann, conTF = conTF, interactions = self.with_interactions, siteOpening=chrst) for ann, chrst in zip(annotation, chromatinState)]) - self.betas*concentration
		y_next = concentration + self.h*der
		y_next[y_next < 0] = 0.0
		#2) calculate change in chromatin opening
		if(fullSeq):
			chromatinState = chromatin_state_all
		chr_der = [self.scale*(self.changeChromatin(seq, mat, conTF) + self.opening_basal_rate*chrst) for seq, mat, chrst in zip(annotation, self.open_precalc, chromatinState)]
		chr_next = [chrom_i + self.h*der_i for chrom_i, der_i in zip(chromatinState, chr_der)]
		for i in range(len(chr_next)):
			chr_next[i][chr_next[i] < self.min_opening] = self.min_opening
			chr_next[i][chr_next[i] > self.max_opening] = self.max_opening
			if(self.normalize_total):
				chr_next[i] = chr_next[i]/np.sum(chr_next[i])
		return (y_next, chr_next)
	def run(self, annotation, start, storeCon = True, start_chrom=None, fullSeq = False, seqlen=300):
		max_time = int(self.tmax/self.h)
		if(fullSeq and storeCon is False):
			fullSeq = False			#If you don't store change over time, you can't calculate all the sequence points either
		self.precalculateChromatinMatrixes(annotation, fullSeq, seqlen)
		if(start_chrom is None):
			if(not fullSeq):
				chromatin_state = [np.array([self.opening_start for i in range(ann.nsites)]) if ann is not None else np.array([self.opening_start]) for ann in annotation]
			else:
				chromatin_state = [np.array([self.opening_start for i in range(max(seqlen, np.max(ann.sites.getPos()) + self.site_size))]) if ann is not None else np.array([self.opening_start for i in range(seqlen)]) for ann in annotation]
		else:
			## uniform sequence for each gene, but a single starting value can be specified for each gene
			if(not fullSeq):
				chromatin_state = [np.array([start_chrom[gene] for i in range(ann.nsites)])  if ann is not None else np.array([start_chrom[gene]]) for gene, ann in enumerate(annotation)]
			else:
				chromatin_state = [np.array([start_chrom[gene] for i in range(max(seqlen, np.max(ann.sites.getPos()) + self.site_size))]) if ann is not None else np.array([start_chrom[gene] for i in range(seqlen)]) for gene, ann in enumerate(annotation)]
						
		if(storeCon):
			concentrations = np.zeros((len(start), max_time))
			concentrations[:, 0] = start
			chromatin_states = [chromatin_state]
			for t in range(1, max_time):
				a, chromatin_state = self.step(concentration = concentrations[:, t-1], annotation = annotation, conTF = concentrations[self.tfs.source_inds, t-1], chromatinState = chromatin_state, fullSeq = fullSeq)
				concentrations[:, t] = a
				chromatin_states.append(chromatin_state)
				if(np.max(np.abs(concentrations[:, t] - concentrations[:, t-1])) <= self.min_variation*self.h):
					concentrations = concentrations[:, range(0, t+1)]
					break
			if(fullSeq):
				return ChromatinResult(concentrations, chromatin_states)
			else:
				return concentrations
		else:
			y_last = start
			for t in range(1, max_time):
				concentrations, chromatin_state = self.step(concentration = y_last, annotation = annotation, conTF = y_last[self.tfs.source_inds], chromatinState = chromatin_state)
				if(np.max(np.abs(y_last - concentrations)) <= self.min_variation):
					break
				y_last = concentrations
			return concentrations	
		return concentrations


class TFSet():
	""" Holds a set of TF with all its data, including an Annotator object"""
	INHIBITOR_TYPE = 0
	ACTIVATOR_TYPE = 1
	PROB_INH = 0.25
	DEF_ALPHA_INH = 0.005
	DEF_ALPHA_ACT = 2#4
	DEF_KMAX = 1
	DEF_KMAX_ACT = DEF_KMAX
	DEF_KMAX_INH = 1*DEF_KMAX
	DEF_INT_PROB_ACTIVATORS = 0.2
	DEF_INT_PROB_INHIBITORS = 0.1
	DEF_INT_STRENGTH = 3
	def __init__(self, inds = None, interactions = None, alphasbtm = None, direction = None, kmax = None, annotator = None, difKmax = False, min_inhibitor_index = 4):
		if(annotator == None and inds != None):
			self.source_inds = inds
			self.inds = range(len(inds))
			self.annotator = SeqAnnotator(self.inds)
		elif(annotator == None and inds == None):
			self.annotator = None
			self.inds = None
		elif(annotator != None and inds != None):
			assert(annotator.inds == inds)
			self.annotator = annotator
			self.inds = self.annotator.inds
		elif(annotator != None and inds == None):
			self.annotator = annotator
			self.inds = self.annotator.inds
		self.min_inhibitor_index = min_inhibitor_index
		if(alphasbtm == None and direction == None):
			self.__setRandAlphas()
		elif(alphasbtm == None):
			self.direction = direction
			self.__setAlphasFromDir()
		elif(isinstance(alphasbtm, int)):
			self.alphasbtm = np.full(len(inds), alphasbtm)
			self.__setDirFromAlphas()
		else:
			assert(len(alphasbtm) == len(self.inds))
			self.alphasbtm = np.array(alphasbtm)
		if(kmax == None):
			if(difKmax):
				self.kmax = [self.DEF_KMAX_ACT if d==self.ACTIVATOR_TYPE else self.DEF_KMAX_INH for d in self.direction]
			else:
				self.kmax = [self.DEF_KMAX_ACT for d in self.direction]
			self.kmax = np.array(self.kmax)
		elif(isinstance(kmax, int)):
			self.kmax = [kmax for d in self.direction]
			self.kmax = np.array(self.kmax)
		else:
			assert(len(kmax) == len(self.inds))
			self.kmax = np.array(kmax)
		if(interactions is not None):
			self.originalInteractions = interactions
			self.interactions = self.interactionsToInnerIndex()
		else:
			self.originalInteractions = self.__generateRandomInteractions()
			self.interactions = self.originalInteractions
	def interactionsToInnerIndex(self):
		newInts = []
		for i in range(len(self.originalInteractions)):
			newInts.append(interactionTuple(self.getInternalIndex(self.originalInteractions[i].tf1), 
							self.getInternalIndex(self.originalInteractions[i].tf2),
							self.originalInteractions[i].weight))
		return newInts
	def setInteractions(self, new_interactions):
		self.originalInteractions = new_interactions
		self.interactions = self.interactionsToInnerIndex()	
	def annotate(self, seqstrvec, threshold = 0.7):
		bindingsites = self.annotator.searchAllpssms(seqstrvec, threshold)
		annotatedseqs = []
		for i in range(len(bindingsites)):
			if(bindingsites[i].length() > 0):
				#Q = sites[i].getEfSite()*self.kmax
				annseq = AnnotatedSeq(ind = i, sites = bindingsites[i], interactions = self.interactions,
						Qonpartial = self.alphasbtm[bindingsites[i].getTFs()] * bindingsites[i].getEfSite() * self.kmax[bindingsites[i].getTFs()], 
						Qoffpartial = bindingsites[i].getEfSite() * self.kmax[bindingsites[i].getTFs()])
				annotatedseqs.append(annseq)
			else:
				annotatedseqs.append(None)
		return annotatedseqs
	def addAlphaError(self, mean = 0, var = 2):
		self.alphasbtm += np.random.normal(mean, var, len(self.inds))
		self.alphasbtm = np.absolute(self.alphasbtm)
		self.alphasbtm = [self.DEF_ALPHA_INH if alpha == 0 else alpha for alpha in self.alphasbtm]
		self.direction = [self.INHIBITION_TYPE if alpha < 1 else self.ACTIVATOR_TYPE for alpha in self.alphasbtm]
	def getPSSMs(self):
		return self.annotator.pssms
	def getTrialSeq(self, n=1, strand = None, randomize = 0, inter = 'rand', max_inter=50):
		return self.annotator.generateSeqsWithConsensus(n, strand, randomize, inter, max_inter)
	def getSourceIndex(self, indexes):
		if(isinstance(indexes, int)):
			if(0 <= indexes < len(self.source_inds)):
				return self.source_inds[indexes]
			else:
				return None
		else:
			sinds = []
			for i in indexes:
				if(0 <= i < len(self.source_inds)):
					sinds.append(self.source_inds[i])
				else:
					sinds.append(None)
			return sinds
	def getInternalIndex(self, indexes):
		if(isinstance(indexes, int)):
			for j in self.inds:
				if(self.source_inds[j] == indexes):
					return j
			return None
		else:	
			intind = []
			for i in indexes:
				for j in self.inds:
					if(self.source_inds[j] == i):
						intind.append(j)
						break
				else:
					intind.append(None)
			return intind	
	def getPSSMs(self):
		return self.annotator.pssms	
	def getPSSMFromSourceIndex(self, indexes):
		return self.annotator.pssms[self.getInternalIndex(indexes)]
	def interTupleToInternal(self, inttuples):
		newtuples = []
		for t in inttuples:
			newtuples.append(interactionTuple(self.getInternalIndex(t.tf1), self.getInternalIndex(t.tf2), t.weight))
		return newtuples
	def interTupleToExternal(self, inttuples):
		newtuples = []
		for t in inttuples:
			newtuples.append(interactionTuple(self.getSourceIndex(t.tf1), self.getSourceIndex(t.tf2), t.weight))
		return newtuples
	def __setRandAlphas(self):
		self.direction = np.random.choice((self.INHIBITOR_TYPE, self.ACTIVATOR_TYPE), size = len(self.inds), replace = True, p = (self.PROB_INH, 1-self.PROB_INH))
		for i in range(self.min_inhibitor_index):		
			self.direction[i] = self.ACTIVATOR_TYPE
		self.__setAlphasFromDir()
	def __setAlphasFromDir(self):
		self.alphasbtm = np.array([self.DEF_ALPHA_INH if i == self.INHIBITOR_TYPE else self.DEF_ALPHA_ACT for i in self.direction]) #Make np.array?
	def __setDirFromAlphas(self):
		self.direction = np.array([self.INHIBITOR_TYPE if i < 1 else self.ACTIVATOR_TYPE for i in self.alphasbtm])
	def __generateRandomInteractions(self):
		inter_list = []
		act = np.where(self.direction == self.ACTIVATOR_TYPE)[0]
		inh = np.where(self.direction == self.INHIBITOR_TYPE)[0]
		act_n_inters = np.random.poisson(self.DEF_INT_PROB_ACTIVATORS*len(act), 1)[0]
		inh_n_inters = np.random.poisson(self.DEF_INT_PROB_INHIBITORS*len(inh), 1)[0]
		for i in range(act_n_inters):
			inter_list.append(interactionTuple(np.random.choice(act, 1)[0], np.random.choice(act, 1)[0],self.DEF_INT_STRENGTH))
		for i in range(inh_n_inters):
			inter_list.append(interactionTuple(np.random.choice(inh, 1)[0], np.random.choice(inh, 1)[0],self.DEF_INT_STRENGTH))
		return inter_list

class SeqAnnotator():
	""" 	Creates and holds a set of PSSMS".	
		Given sequences, annotates matches to its member TFs
		TFinds will hold the genome indexes of the TF
		pssms will hold a tuple with FW and RV PWMs (to avoid calculating rv each time)
		alphabtms are alphas for the expression thermodynamic model
		alphas and background are parameters to generate random pssms 
	FIELDS:
		- POS (constant). Is set to 0. Positive strand is equivalent to number 0
		- NEG (constant). Is set to 1. Negative strand is equivalent to number 1
		- MAX_NUMBER_OF_SITES: If finds a greater number of binding sites above threshold, selects the first MAX_NUMBER_OF_SITES and omits the rest
		- inds: TF indexes (usually as stored internally in a TFSet object, i.e, not necessarily equal to those stored in a genome object)
		- pssms: a vector of motifs.pssm type objects. Length should be equal to inds length. If they are not passed to the constructor, they are randomly generated through the MotifGenerator class
		- TFmeta: can contain other TF data, but it is not used for the moment

	METHODS:
		-CONSTRUCTOR: needs a LIST of indexes of TF. TFSet objects will convert and use indexes as 0, 1...n, so be careful when interacting with them.
			Optionally: 
					- pssms: list of PSSMs of the TFs. If None will be initialized randomly
					- alphas: Numeric List of length equal to TFinds length. Each element is the Dirichlet prior parameter to generate randomly one column of the PWM (the more close to the center of the motif, the higher the alpha should be). 
					- background: genomic background. Usually is set uniform
					- TFMeta: other TF data. Not implemented any use for it for the moment.
		- searchpssm: searches one PSSM in a sequence. First calculates the probability at each position and then filters out the positions with probability lower than a threshold. Takes:
				-pssm: a pssm object (motifs library)
				-seqstr: a string of characters 'ACGT' containing sequence to search motif
				-strand: strand where is searching. Positive = 0, Negative = 1. Only is used to format the output data
				-tf: TF index of the TF which PSSM is being matched against the sequence. Only used to format the output
				-threshold: float between 0 and 1. Defaults to 0.9. threshold for the energy of a pssm match with respect to the PSSM consensus energy. 
				@returns an object of class SeqAnnotation
		- searchAllpssms: searches all PSSMs in a list in all the sequences in a list of strings, both in the + and - strand. Takes:
				- seqstrvec: a vector of strings with only 'ACGT'
				- threshold: float number to pass to searchpssm function
				@returns a list of SeqAnnotation objects. In each object, matches for all PSSMs, both in the + and - strand, are concatenated
		- getConsensus: @returns the consensus for all PSSMs in self.pssms list
		- generateSeqsWithConsensus:Generates a list of sequences with the consensus for all the object PSSMs concatenated. Optionally it can randomize the sequence, each nucleotide with certain probability. Takes:
				- n: Number of sequences to generate. If randomize = 0, all the sequences will be identical
				- strand: If 0 or 1 will only concatenate the FW or RV motifs, respectively. Otherwise, both will be added, one after another.
				- randomize: the probability that each nucleotide is mutated. Defauts to 0
				@ returns: a list of strings which only contain 'ACGT'
"""

	POS = 0
	NEG = 1
	MAX_NUMBER_OF_SITES = 20
	def __init__(self, TFinds, pssms = None, alphas = None, background = None, TFmeta = None):  #Remove alphasbtm
		self.TFmeta = TFmeta
		self.inds = TFinds
		if(pssms == None):
			self.pssms = []
			for i in range(len(self.inds)):
				self.pssms.append(MotifGenerator.getPSSM(alphas, background))
		else:
			assert(len(pssms) == len(self.inds))
			self.pssms = pssms
	@classmethod
	def searchpssm(self, pssm, seqstr, strand = 0, tf = 0, threshold = 0.9):
		vals = np.array(pssm.calculate(Seq(seqstr, pssm.alphabet)), dtype = np.float)
		positions = np.where(vals/pssm.max >= threshold)
		if(len(positions)>self.MAX_NUMBER_OF_SITES):
			positions = positions[0:(self.MAX_NUMBER_OF_SITES-1)]
		if (len(positions[0])>0):
			vals = vals[positions]
			res = SeqAnnotation(len(vals))
			res.assign('tf_ind', np.full(len(vals), tf))
			res.assign('positions',  positions[0])
			res.assign('LLR', vals)
			res.assign('percent', vals/pssm.max)
			res.assign('strand', np.full(len(vals), strand))
			res.assign('efsite', np.exp(vals - pssm.max))
			return res
		else:
			return None	
	def searchAllpssms(self, seqstrvec, threshold = 0.9, restrict = True):
		allseqsAnnotations = []
		for i in range(len(seqstrvec)):
			seqsites = SeqAnnotation()
			for tf in range(len(self.pssms)):
				annotation = self.searchpssm(self.pssms[tf].pos, seqstrvec[i], strand = self.POS, tf = tf, threshold = threshold)
				if(annotation is not None):
					seqsites.append(annotation)
				annotation = self.searchpssm(self.pssms[tf].neg, seqstrvec[i], strand = self.NEG, tf = tf, threshold = threshold)
				if(annotation is not None):
					seqsites.append(annotation)				
			allseqsAnnotations.append(seqsites)
		return allseqsAnnotations
	def getConsensus(self):
		cons = []
		for i in range(len(self.pssms)):
			cons.append([str(self.pssms[i].pos.consensus), str(self.pssms[i].neg.consensus)])
		return cons
	def generateSeqsWithConsensus(self, n=1, strand = None, randomize = 0, inter = 'rand', max_inter=50):
		cons = self.getConsensus()
		strings = []
		for i in range(n):
			s = ''
			## Add motif consensus (in one strand or both strands)
			if(strand != self.POS and strand != self.NEG):
				for i in range(len(cons)):
					s += cons[i][self.POS] + cons[i][self.NEG] + self.__addSeqBetweenMotifs(inter, max_inter)
			else:
				for i in range(len(cons)):
					s += cons[i][strand] + self.__addSeqBetweenMotifs(inter, max_inter)
			## Finally, randomize nucleotides with certain probability, in case you want to get some degenerate sites.
			if(randomize > 0):
				s = np.array(list(s))
				ind = np.where(np.random.choice((True, False), size= len(s), p = np.array([randomize, 1-randomize])))[0]
				s[ind] = np.random.choice(MotifGenerator.ALPHABET, ind.shape[0], p = list(MotifGenerator.DEFAULT_BACKGROUND.values()))				
			strings.append(''.join(list(s)))
		return strings
	def __addSeqBetweenMotifs(self, inter, max_inter):
		## Add (probably) inactive sequence between motifs
		if(inter == 'rand'):
			# of random length with max_inter maximum length
			s = 'A'*np.random.choice(range(max_inter),1 )[0]
		elif(inter != ''):
			# if inter is numeric, add seq with length indicated by number
			if(inter.isnumeric()):
				s = 'A'*int(inter)
			# if inter is not numeric, the length of the string will be the length of nucleotides added
			else:
				s = 'A'*len(inter)
		return s









##########################

"""

def anntry():
	ann = SeqAnnotator(range(10))
	s1 = ''
	s2 = ''
	s3 = ann.pssms[0].pos.consensus + 'ATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT'
	s4 = 'ATCGACAGACGCGATACGTGCAGTAGATCGTATGACTAGCCTAG'
	for i in range(10):
		s1 += ann.pssms[i].pos.consensus
		s2 += ann.pssms[i].pos.consensus
	sv = [s1.tostring(), s2.tostring(), s3.tostring(), s4]
	print(sv)
	rr = ann.searchAllpssms(sv)
	print ('\n***\n')
	for i in range(len(rr)):
		print(sv[i] + "\n")
		print(rr[i].arr.view())
	return (ann, rr)



x = anntry()



tfs = gx.TFSet([0, 3, 5, 7, 9], alphasbtm = None, direction = [1,1,1,0,1],interactions = [gx.interactionTuple(3, 5, 1)])
ss=tfs.getTrialSeq(n=10, inter='10', randomize=0.1)
ann = tfs.annotate(ss)
ex = gx.ExpressionCalculator()
e1 = []
t = time.time()
times = [t]
for i in range(len(ann)):
	e1.append(ex.activate(ann[i]))
	times.append(time.time()-times[i-1])
finalt = time.time() - t


start1= [1, 0,0,0,0, 0, 0, 0, 0, 0]
start2= [1, 0,0,1,0, 0, 0, 0, 0, 0]
start3= [1, 0,0,1,0, 0, 0, 7, 0, 0]
oder = gx.ODERunner(tfs, np.array([0.1 for i in ss]),with_interactions = False, h=1, tmax=10)
final_expr1 = oder.run(ann, start1); final_expr1
final_expr2 = oder.run(ann, start2); final_expr2
final_expr3 = oder.run(ann, start3); final_expr3
oderX = gx.ODERunnerChromatin(tfs, np.array([0.1 for i in ss]),with_interactions = False, h=1, tmax=10, amplitude=5)
final_expr4 = oderX.run(ann, start1); final_expr4
final_expr5 = oderX.run(ann, start2); final_expr5
final_expr6 = oderX.run(ann, start3); final_expr6


print('Final time for ' + str(len(tfs.inds)+1) + ' sequences of about 130bp: ' + str(finalt))
print('Mean per seq: ' + str(finalt/len(ann)))
print('Individual times: ')
print (times)

a = tfs.getInternalIndex(3)
b = tfs.getInternalIndex([3, 3, 5, 10])
c = tfs.getInternalIndex(18)
d = tfs.getInternalIndex([3, 3, 18, 10])


e = tfs.getSourceIndex(0)
f = tfs.getSourceIndex([0, 0, 3, 4])
g = tfs.getSourceIndex(20)
h = tfs.getSourceIndex([0, 3, 20, 4])

print(a, b, c, d, e, f, g, h)
##### tests




#### test interactions
tfs = TFSet([3,5], alphasbtm = 2, interactions = [interactionTuple(3, 5, 2)])
ss=tfs.getTrialSeq(n=1, strand = 0, inter='10')
ex = ExpressionCalculator()

ann = tfs.annotate(ss)
ex1 = ex.activate(ann[0], conTF = [2, 1])
ex2 = ex.activate(ann[0], conTF = [2, 1],interactions = True)
print('expression 1 is: ', ex1, ' Expression with interactions = 2 is: ', ex2)
tfs.setInteractions([interactionTuple(3, 5, 5)])
ann = tfs.annotate(ss)
ex1 = ex.activate(ann[0])
ex2 = ex.activate(ann[0], interactions = True)
print('expression 1 is: ', ex1, ' Expression with interactions = 5 is: ', ex2)
print (tfs.interactions)

### test sequences without TFBS
ss=tfs.getTrialSeq(n=1, strand = 0, randomize = 1)
ex = ExpressionCalculator()
ann = tfs.annotate(ss)

"""

	
