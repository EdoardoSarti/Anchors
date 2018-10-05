#!/Users/Orfeo/anaconda3/bin/python3

import argparse, re, collections, time, subprocess, os
import numpy as np
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument('-i1', '--input_sequence_1', nargs=1)
parser.add_argument('-i2', '--input_sequence_2', nargs=1)
parser.add_argument('-s', '--simscorefile', nargs=1)
parser.add_argument('-o', '--outputfile', nargs='?', default='NW_score.txt')
parser.add_argument('-oo', '--outputfile_optimization', nargs='?', default='alignments.txt')
parser.add_argument('-Op', '--gap_opening_penalty', nargs='?')
parser.add_argument('-Ext', '--gap_extension_penalty', nargs='?')
parser.add_argument('-OpAbove', '--above_threshold_gap_opening_penalty', nargs='?')
parser.add_argument('-OpBelow', '--below_threshold_gap_opening_penalty', nargs='?')
parser.add_argument('-ExtAbove', '--above_threshold_gap_extension_penalty', nargs='?')
parser.add_argument('-ExtBelow', '--below_threshold_gap_extension_penalty', nargs='?')
parser.add_argument('-TermOp', '--termini_gap_opening_penalty', nargs=1)
parser.add_argument('-TermExt', '--termini_gap_extension_penalty', nargs=1)
parser.add_argument('-thr', '--threshold', nargs='?')
parser.add_argument('-optimize', dest='optimization', action='store_true')
parser.set_defaults(optimization=False)
parsed = parser.parse_args()

print("STOP!!! THERE PROBABLY IS A SIGN ERROR WHEN UPDATING THE POTENTIAL IN ALL THE ASSOCIATE_XXX ROUTINES!!! CHECK!!!")


class RegExp():
	def __init__(self):
		self.emptyline = re.compile('^\s*$')
		self.headerline = re.compile('^>.*$')

class Element():
	def __init__(self, value=-1000000000000):
		self.paths = np.ones(9)*value
		self.history = ''
		self.directions = [
		                    'D',
		                    'H',
		                    'V'
		                  ]
		self.path_names = {
		                    'D>D' : 0,
                                    'H>D' : 1,
                                    'V>D' : 2,
                                    'D>H' : 3,
                                    'H>H' : 4,
                                    'V>H' : 5,
                                    'D>V' : 6,
                                    'H>V' : 7,
                                    'V>V' : 8
		                  }

	def zeros(self):
		self.paths = np.zeros(9)
	
	def add_path_value(self, value, path_name=None, path_id=None):
		if path_id:
			self.paths[path_id] = value
		elif path_name:
#			print(path_name, value)
			self.paths[self.path_names[path_name]] = value
		else:
			raise NameError("ERROR: when adding a value to a NW element you should specify either path_name or path_id")

	def best_subpath_value(self, subpath_name=None, subpath_id=None):
		if subpath_id:
			return np.amax(self.paths[subpath_id*3 : subpath_id*3 + 3])
		elif subpath_name:
			i = self.directions.index(subpath_name)
#			path_ind = np.argwhere(self.paths[i*3 : i*3 + 3] == np.amax(self.paths[i*3 : i*3 + 3])).flatten()
#			if(len(path_ind)>1):
#				print(path_ind, np.amax(self.paths[i*3 : i*3 + 3]))
#			print(subpath_name, np.amax(self.paths[i*3 : i*3 + 3]), self.paths[i*3 : i*3 + 3])
			return np.amax(self.paths[i*3 : i*3 + 3])
		else:
			raise NameError("ERROR: when calculating the best subpath of a NW element you should specify either subpath_name or subpath_id")

	def best_path(self):
		return np.argmax(self.paths), np.amax(self.paths)

	def best_group_subpath(self):
		path_ind = np.argwhere(self.paths == np.amax(self.paths)).flatten()
		group = path_ind/3.0
		g = list(set(group.astype(int).tolist()))
		gg = [ self.directions[g[i]] for i in range(len(g)) ]
		return gg, np.amax(self.paths)

class Gaps():
	def __init__(self):
		self.gap = {
		             'Op' : 0.0,
		             'Ext' : 0.0,
		             'OpAbove' : 0.0,
		             'OpBelow' : 0.0,
		             'ExtAbove' : 0.0,
		             'ExtBelow' : 0.0,
		             'TermOp' : 0.0,
                             'TermExt' : 0.0
		           }

	def add_gap_value(self, gap_type, value):
		self.gap[gap_type] = value

def weighted_moving_average(x,y,bin_center,width,windowtype):

	def identity(x,wsize,mean,amp=1):
		y = np.zeros(len(x))
		y[mean] = 1
		return amp*y

	def triangular(x,wsize,mean,amp=1):
		y = np.zeros(len(x))
		for i in range(len(x)):
			if abs(mean - i) > (wsize-1)/2:
				y[i] = 0
			else:
				y[i] = ((wsize+1)/2 - abs(mean - i)) / ((wsize + 1)/2)
		return amp*y

	def rectangular(x,wsize,mean,amp=1):
		y = np.zeros(len(x))
		for i in range(len(x)):
			if abs(mean - i) > (wsize-1)/2:
				y[i] = 0
			else:
				y[i] = 1
		return amp*y

	wtypes = {
                   'none' : identity,
                   'triangular' : triangular,
                   'triangular_msa' : triangular,
                   'rectangular' : rectangular
                  }

	weights = wtypes[windowtype](x,width,bin_center)
	bin_avg = np.average(y,weights=weights)

	return bin_avg

def check_attributes(record_fields, nr, *args):
	for i in args:
		if not i in record_fields:
			raise NameError("ERROR: Record {0} does not have a \"{1}\" keyword".format(nr, i))

def plot_inputs(plottype, nr, record_fields, seqs, instr_1, instr_2, *args):
	filename = plottype + '_' + str(nr) + '.png'
	if (plottype == 'substitution_matrix' or plottype == 'pssm'):
		fig = plt.figure(nr, (15,6))
	else:
		fig = plt.figure(nr, (10,6))
	if (plottype == 'scale' or plottype == 'substitution_matrix'):
		fig.suptitle(plottype + ': ' + record_fields['file:'].split('/')[-1])
	else:
		fig.suptitle(plottype + ': similarity score file line ' + str(nr))

	if (plottype == 'substitution_matrix' or plottype == 'pssm'):
		aspect = 1
		ax = plt.subplot(211, adjustable='box', aspect=aspect)
	else:
		ax = plt.subplot(211)
	plt.title(seqs[0][0])
	if (plottype == 'scale'):
		ax.plot(instr_1[0], instr_1[1], ':', color='0.5')
		ax.plot(instr_1[2], instr_1[3], lw=2.0, color='b', label='window averaging, size='+str(args[0]))
#		print(list(zip(instr_1[2], instr_1[3])))
		ax.legend()
	elif (plottype == 'substitution_matrix' or plottype == 'pssm'):
		ax.pcolor(instr_1[0], cmap=plt.cm.Blues)
	else:
		ax.plot(instr_1[0], instr_1[1])
	ax.set_xlim(0, len(seqs[0][1]))
	
	if (plottype == 'substitution_matrix' or plottype == 'pssm'):
		aspect = 1
		ax = plt.subplot(212, adjustable='box', aspect=aspect)
	else:
		ax = plt.subplot(212)
	plt.title(seqs[1][0])
	if (plottype == 'scale'):
		ax.plot(instr_2[0], instr_2[1], ':', color='0.5')
		ax.plot(instr_2[2], instr_2[3], lw=2.0, color='b', label='window averaging, size='+str(args[0]))
#		print(list(zip(instr_2[2], instr_2[3])))
		ax.legend()
	elif (plottype == 'substitution_matrix' or plottype == 'pssm'):
		ax.pcolor(instr_2[0], cmap=plt.cm.Blues)
	else:
		ax.plot(instr_2[0], instr_2[1])
	ax.set_xlim(0, len(seqs[1][1]))

	plt.savefig(filename)
		

def associate_scale(potential, seq_1, seq_2, record_fields, nr):
	RE = RegExp()
	check_attributes(record_fields, nr, 'file:', 'windowtype:', 'windowsize:')

	f = open(record_fields['file:'], 'r')
	text = f.read().split('\n')
	scale = {}
	for nl in range(len(text)):
		if RE.emptyline.match(text[nl]):
			continue
		elif (len(text[nl].split()) == 2):
			scale[text[nl].split()[0]] = float(text[nl].split()[1])
	wsize = float(record_fields['windowsize:'])
	wtype = str(record_fields['windowtype:'])
	scaleprofile_1 = [ scale[seq_1[1][i]] for i in range(len(seq_1[1])) ]
	scaleprofile_2 = [ scale[seq_2[1][i]] for i in range(len(seq_2[1])) ]
	winav_1 = [ weighted_moving_average(range(len(scaleprofile_1)),scaleprofile_1,i,wsize,wtype) for i in range(len(seq_1[1])) ]
	winav_2 = [ weighted_moving_average(range(len(scaleprofile_2)),scaleprofile_2,i,wsize,wtype) for i in range(len(seq_2[1])) ]

	for na1 in range(len(seq_1[1])):
		for na2 in range(len(seq_2[1])):
			potential[na1, na2]['D'] -= float(record_fields['weight:']) * abs(winav_1[na1] - winav_2[na2])

	plot_inputs( 
	             'scale',
	             nr,
	             record_fields,
	             (seq_1, seq_2),
	             [range(len(scaleprofile_1)), scaleprofile_1, range(len(winav_1)), winav_1],
	             [range(len(scaleprofile_2)), scaleprofile_2, range(len(winav_2)), winav_2],
	             wsize
	           )
	return potential


def associate_seqsim(potential, seq_1, seq_2, record_fields, nr):
	RE = RegExp()
	check_attributes(record_fields, nr, 'file:')

	f = open(record_fields['file:'], 'r')
	text = f.read().split('\n')
	matrix = {}
	for nl in range(len(text)):
		fields = text[nl].split()
		if RE.emptyline.match(text[nl]):
			continue
		elif (len(fields) == 20):
			col_names = []
			for nf in range(len(fields)):
				col_names.append(fields[nf])
		elif (len(fields) == 21):
			row_name = str(fields[0])
			for nf in range(1, len(fields)):
				matrix[(row_name, col_names[nf-1])] = float(fields[nf])
	
	matrixprofile_1 = np.empty( (20, len(seq_1[1])) )
	for i in range(0, 20):
		for j in range(len(seq_1[1])):
			matrixprofile_1[i, j] = matrix[(seq_1[1][j], col_names[i])]
	matrixprofile_2 = np.empty( (20, len(seq_2[1])) ) 
	for i in range(0, 20):
		for j in range(len(seq_2[1])):
			matrixprofile_2[i, j] = matrix[(seq_2[1][j], col_names[i])]

	for na1 in range(len(seq_1[1])):
		for na2 in range(len(seq_2[1])):
			potential[na1, na2]['D'] -= float(record_fields['weight:']) * matrix[(seq_1[1][na1], seq_2[1][na2])] 

	plot_inputs(
	             'substitution_matrix',
	             nr,
	             record_fields,
	             (seq_1, seq_2),
	             [matrixprofile_1],
                     [matrixprofile_2]
                   )
	return potential

def associate_profile(potential, seq_1, seq_2, record_fields, nr, threshold=None):
	RE = RegExp()
	check_attributes(record_fields, nr, 'column:', 'headerlines:', 'profile1:', 'profile2:')

	c = int(record_fields['column:'])-1

	f = open(record_fields['profile1:'], 'r')
	text = f.read().split('\n')
	il, hl = 0, 0
	profile_1 = []
	for nl in range(len(text)):
		if RE.emptyline.match(text[nl]):
			continue
		elif hl < int(record_fields['headerlines:']):
			hl += 1
			continue
		else:
			profile_1.append(float(text[nl].split()[c]))
			il += 1

	f = open(record_fields['profile2:'], 'r')
	text = f.read().split('\n')
	il, hl = 0, 0
	profile_2 = []
	for nl in range(len(text)):
		if RE.emptyline.match(text[nl]):
			continue
		elif hl < int(record_fields['headerlines:']):
			hl += 1
			continue
		else:
			profile_2.append(float(text[nl].split()[c]))
			il += 1

	for na1 in range(len(seq_1[1])):
		for na2 in range(len(seq_2[1])):
			potential[na1, na2]['D'] -= float(record_fields['weight:']) * abs(profile_1[na1] - profile_2[na2]) 
	
	plot_inputs(
	             'profile',
	             nr,
	             record_fields,
	             (seq_1, seq_2),
	             [range(len(profile_1)), profile_1],
	             [range(len(profile_2)), profile_2]
	           )

	if threshold:
		binary_profile = {}
		for na1 in range(len(seq_1[1])):
			binary_profile[na1] = {}
			for na2 in range(len(seq_2[1])):
				if abs(profile_1[na1] - profile_2[na2]) >= threshold:
					binary_profile[na1][na2] = 'Above'
				else:
					binary_profile[na1][na2] = 'Below'

		return potential, binary_profile

	return potential


def associate_PSSM(potential, seq_1, seq_2, record_fields, nr):
	RE = RegExp()
	check_attributes(record_fields, nr, 'pssm1:', 'pssm2:')

	naa = None
	f = open(record_fields['pssm1:'], 'r')
	text = f.read().split('\n')
	matrixprofile_1 = np.empty( (20, len(seq_1[1])) )
	for nl in range(len(text)):
		fields = text[nl].split()
		if RE.emptyline.match(text[nl]):
			continue
		elif (naa == None and (len(fields) == 20 or len(fields) == 40) and max([len(fields[i]) for i in range(len(fields))])==1):
			col_names = []
			for nf in range(20):
				col_names.append(fields[nf])
			naa = 0
		elif (naa != None and naa >=0 and int(fields[0])==naa+1):
			for nf in range(2, 22):
#				print(nf-2, naa, float(fields[nf]))
				matrixprofile_1[nf-2, naa] = float(fields[nf])
			naa += 1
		if naa != None and naa >= len(seq_1[1]):
			break
	naa = None
	f = open(record_fields['pssm2:'], 'r')
	text = f.read().split('\n')
	matrixprofile_2 = np.empty( (20, len(seq_2[1])) )
	for nl in range(len(text)):
		fields = text[nl].split()
		if RE.emptyline.match(text[nl]):
			continue
		elif ((len(fields) == 20 or len(fields) == 40) and max([len(fields[i]) for i in range(len(fields))])==1):
			col_names = []
			for nf in range(20):
				col_names.append(fields[nf])
			naa = 0
		elif (naa != None and naa >=0 and int(fields[0])==naa+1):
			for nf in range(2, 22):
				matrixprofile_2[nf-2, naa] = float(fields[nf])
			naa += 1
		if naa != None and naa >= len(seq_2[1]):
			break



	for na1 in range(len(seq_1[1])):
		for na2 in range(len(seq_2[1])):
#			print(na1, seq_1[1][na1], col_names.index(seq_1[1][na1]), na2, matrixprofile_2[col_names.index(seq_1[1][na1]), na2])
			potential[na1, na2]['D'] += float(record_fields['weight:']) * 0.5*(matrixprofile_1[col_names.index(seq_2[1][na2]), na1] + matrixprofile_2[col_names.index(seq_1[1][na1]), na2])

	plot_inputs(
	             'pssm',
	             nr,
	             record_fields,
	             (seq_1, seq_2),
	             [matrixprofile_1],
                     [matrixprofile_2]
                   )
	return potential

def impose_anchors(potential, seq_1, seq_2, record_fields, nr, readonly=False):
	RE = RegExp()
	check_attributes(record_fields, nr, 'file:')

	f = open(record_fields['file:'], 'r')
	text = f.read().split('\n')
	anchors = collections.OrderedDict()
	anch_boundaries = collections.OrderedDict()
	for nl in range(len(text)):
		fields = re.sub('[\s+]', ' ', text[nl]).split()
#		print(fields)
		if nl == 0:
			n_anchors = int(fields[0])
		elif nl <= n_anchors:
			if len(fields) == 2:
				anchors[(int(fields[0]), int(fields[1]))] = None
			elif len(fields) == 3:
				anchors[(int(fields[0]), int(fields[1]))] = float(fields[2])
			elif len(fields) == 7:
				anchors[(int(fields[0]), int(fields[1]))] = float(fields[2])
				anch_boundaries[(int(fields[0]), int(fields[1]))] = ((int(fields[3]), int(fields[4])), (int(fields[5]), int(fields[6])))

	if readonly:
		return anchors, anch_boundaries

	el = Element()
	for na1 in range(len(seq_1[1])):
		for na2 in range(len(seq_2[1])):
			for a in anchors:
				if (na1 < a[0] and na2 >= a[1]) or (na1 >= a[0] and na2 < a[1]):
					for d in el.directions:
						potential[na1, na2][d] -= float(record_fields['weight:']) * anchors[a]
	return potential

def to_fasta(sequences, alignment, anchors=None, enable_DHV=False):
	if alignment[0] != (0, 0):
		raise NameError("ERROR: The alignment to be translated must start at (0,0)")

	aligned_seqs = ([], [])
	matches = []
	nD, nH, nV = 0, 0, 0
	i_old, j_old = alignment[0][0], alignment[0][1]
	for k in range(1, len(alignment)):
		i, j = alignment[k][0], alignment[k][1]
		if i == i_old and j != j_old:
			aligned_seqs[0].append('-')
			aligned_seqs[1].append(sequences[1][j-1])
			nH += 1
		elif i != i_old and j == j_old:
			aligned_seqs[0].append(sequences[0][i-1])
			aligned_seqs[1].append('-')
			nV += 1
		elif i != i_old and j != j_old:
			aligned_seqs[0].append(sequences[0][i-1])
			aligned_seqs[1].append(sequences[1][j-1])
			nD += 1
		i_old, j_old = i, j
		if aligned_seqs[0][-1] == aligned_seqs[1][-1]:
			matches.append(':')
		else:
			matches.append(' ')

	if enable_DHV:
		return nD, nH, nV

	signatures = ['', '']
	if anchors == None:
		return aligned_seqs, matches, signatures

	na = 0
	pos_anch = []
	count = [0, 0]
	old = [0, 0]	
	for a in anchors:
		pos_anch.append([-1, -1])
		dcount = [0, 0]
		for s in [0, 1]:
			for i in range(old[s], len(aligned_seqs[s])):
				if aligned_seqs[s][i] != '-':
					if count[s] == a[s]:
						pos_anch[na][s] = dcount[s]
						count[s] += 1
						old[s] = i+1
						break
					else:
						count[s] += 1
						dcount[s] += 1
				else:
					dcount[s] += 1
				old[s] = i+1
		for s in [0, 1]:
			if pos_anch[na][s] == -1:
				pos_anch[na][s] = dcount[s]
		na += 1

	na = 0
	for a in anchors:
		for s in [0, 1]:
			signatures[s] = signatures[s] + ' '*pos_anch[na][s] + str(na+1)
		na += 1

	return aligned_seqs, matches, signatures

def needleman_wunsch(sequence_1, sequence_2, potential, profile, gap):
#	NWM = []
	NWM = np.empty((len(sequence_1)+1, len(sequence_2)+1), dtype=object)

#	print(profile)

	# First line
#	NWline0 = []
	el = Element(0)
#	NWline0.append(el)
	NWM[0, 0] = el
	for j in range(1, len(sequence_2)+1):
		el = Element()
		el.add_path_value(gap.gap['TermOp'] + (j-1)*gap.gap['TermExt'] + sum([ potential[0, k-1]['H'] for k in range(1, j+1)]), path_name='H>H')
		NWM[0, j] = el
#		NWline0.append(el)
#	NWM.append(NWline0)

	# Main loop
	for i in range(1, len(sequence_1)):
#		NWline = []
		el = Element()

		# First column elements
		el.add_path_value(gap.gap['TermOp'] + (i-1)*gap.gap['TermExt'] + sum([ potential[k-1, 0]['V'] for k in range(1, i+1)]), path_name='V>V')
#		NWline.append(el)
		NWM[i, 0] = el
		# Core elements
		for j in range(1, len(sequence_2)):
			el = Element(0)

			# Diagonal values
			d2, el_past = 'D', NWM[i-1, j-1]
			for d1 in el.directions:
				path_name = d1 + '>' + d2
#				print(i,j,path_name, el_past.best_subpath_value(d1), potential[i-1, j-1][d2])
				value = el_past.best_subpath_value(d1) + potential[i-1, j-1][d2]
				el.add_path_value(value, path_name=path_name)

			# Horizontal and vertical values
#			for d2, el_past in [('H', NWline[j-1]), ('V', NWM[i-1][j])]:
			for d2, el_past in [('H', NWM[i, j-1]), ('V', NWM[i-1, j])]:
				# When d1 and d2 are different...
				for d1 in list(set(el.directions) - set(d2)):
					path_name = d1 + '>' + d2
#					print(i,j,path_name, el_past.best_subpath_value(d1), gap.gap[ 'Op' + profile[i-1][j-1] ], potential[i-1, j-1][d2])
					value = el_past.best_subpath_value(d1) + gap.gap[ 'Op' + profile[i-1][j-1] ] + potential[i-1, j-1][d2]
					el.add_path_value(value, path_name=path_name)
				# ...and when they are equal.
				d1 = d2
				path_name = d1 + '>' + d2
#				print(i,j,path_name, el_past.best_subpath_value(d1), gap.gap[ 'Ext' + profile[i-1][j-1] ], potential[i-1, j-1][d2])
				value = el_past.best_subpath_value(d1) + gap.gap[ 'Ext' + profile[i-1][j-1] ] + potential[i-1, j-1][d2]
				el.add_path_value(value, path_name=path_name)
	
#			NWline.append(el)
			NWM[i, j] = el
#			if(i<10 and  j <10):
#				print(i, j, NWM[i, j].paths[0])
#		NWM.append(NWline)
	
	# Last line
#	NWlineT = []
	el = Element()
	i = len(sequence_1)
	el.add_path_value(gap.gap['TermOp'] + (i-1)*gap.gap['TermExt'] + sum([potential[k-1, 0]['V'] for k in range(1, i+1)]), path_name='V>V')
#	NWlineT.append(el)
	NWM[i,0] = el
	for j in range(1, len(sequence_2)):
		el = Element(0)
		
		# Diagonal values
		d2, el_past = 'D', NWM[i-1, j-1]
		for d1 in el.directions:
			path_name = d1 + '>' + d2
			value = el_past.best_subpath_value(d1) + potential[i-1, j-1][d2]
			el.add_path_value(value, path_name=path_name)

		# Horizontal values (terminal gaps)
#		d2, el_past = 'H', NWlineT[j-1]
		d2, el_past = 'H', NWM[i, j-1]
		for d1 in list(set(el.directions) - set(d2)):
			path_name = d1 + '>' + d2
			value = el_past.best_subpath_value(d1) + gap.gap[ 'TermOp' ] + potential[i-1, j-1][d2]
			el.add_path_value(value, path_name=path_name)
		d1 = d2
		path_name = d1 + '>' + d2
		value = el_past.best_subpath_value(d1) + gap.gap[ 'TermExt' ] + potential[i-1, j-1][d2]
		el.add_path_value(value, path_name=path_name)

		# Vertical values (normal gaps)
		d2, el_past = 'V', NWM[i-1, j]
		for d1 in list(set(el.directions) - set(d2)):
			path_name = d1 + '>' + d2
			value = el_past.best_subpath_value(d1) + gap.gap[ 'Op' + profile[i-1][j-1] ] + potential[i-1, j-1][d2]
			el.add_path_value(value, path_name=path_name)
		d1 = d2
		path_name = d1 + '>' + d2
		value = el_past.best_subpath_value(d1) + gap.gap[ 'Ext' + profile[i-1][j-1] ] + potential[i-1, j-1][d2]
		el.add_path_value(value, path_name=path_name)

		NWM[i, j] = el
#		NWlineT.append(el)
#	NWM.append(NWlineT)

	# Last row (last element of each line)
	j = len(sequence_2)
	for i in range(1, len(sequence_1)): 
		el = Element(0)

		# Diagonal values
		d2, el_past = 'D', NWM[i-1, j-1]
		for d1 in el.directions:
			path_name = d1 + '>' + d2
			value = el_past.best_subpath_value(d1) + potential[i-1, j-1][d2]
			el.add_path_value(value, path_name=path_name)

		# Horizontal values (normal gaps)
		d2, el_past = 'H', NWM[i, j-1]
		for d1 in list(set(el.directions) - set(d2)):
			path_name = d1 + '>' + d2
			value = el_past.best_subpath_value(d1) + gap.gap[ 'Op' + profile[i-1][j-1] ] + potential[i-1, j-1][d2]
			el.add_path_value(value, path_name=path_name)
		d1 = d2
		path_name = d1 + '>' + d2
		value = el_past.best_subpath_value(d1) + gap.gap[ 'Ext' + profile[i-1][j-1] ] + potential[i-1, j-1][d2]
		el.add_path_value(value, path_name=path_name)

		# Vertical values (terminal gaps)
		d2, el_past = 'V', NWM[i-1, j]
		for d1 in list(set(el.directions) - set(d2)):
			path_name = d1 + '>' + d2
			value = el_past.best_subpath_value(d1) + gap.gap[ 'TermOp' ] + potential[i-1, j-1][d2]
			el.add_path_value(value, path_name=path_name)
		d1 = d2
		path_name = d1 + '>' + d2
		value = el_past.best_subpath_value(d1) + gap.gap[ 'TermExt' ] + potential[i-1, j-1][d2]
		el.add_path_value(value, path_name=path_name)

#		NWM[i].append(el)
		NWM[i, j] = el

	# Last element (terminal gaps either way)
	i, j = len(sequence_1), len(sequence_2)
	el = Element(0)

	# Diagonal values
	d2, el_past = 'D', NWM[i-1, j-1]
	for d1 in el.directions:
		path_name = d1 + '>' + d2
		value = el_past.best_subpath_value(d1) + potential[i-1, j-1][d2]
		el.add_path_value(value, path_name=path_name)
	
	# Horizontal and vertical values (terminal gaps)
	for d2, el_past in [('H', NWM[i, j-1]), ('V', NWM[i-1, j])]:
		# When d1 and d2 are diffrent...
		for d1 in list(set(el.directions) - set(d2)):
			path_name = d1 + '>' + d2
			value = el_past.best_subpath_value(d1) + gap.gap[ 'TermOp' ] + potential[i-1, j-1][d2]
			el.add_path_value(value, path_name=path_name)
		# ...and when they are equal
		d1 = d2
		path_name = d1 + '>' + d2
		value = el_past.best_subpath_value(d1) + gap.gap[ 'TermExt' ] + potential[i-1, j-1][d2]
		el.add_path_value(value, path_name=path_name)

#	NWM[i].append(el)
	NWM[i, j] = el

#	print('length matrix ', len(NWM[0]))

	return NWM


def traceback(sequence_1, sequence_2, NWM, alignments=None, scores=None, code='0', end_row=None, end_column=None):

	if ((alignments or scores or end_row or end_column) and (not (alignments and scores)) and (end_row == None or end_column == None)):
#		print(alignments[code], scores[code], end_row, end_column)
		raise NameError("ERROR: If traceback is called on a partially completed alignment, you have to provide the partial alignment, score and end-row and -column")

	if (alignments and (len(alignments.keys()) > max([len(sequence_1), len(sequence_2)]))):
		raise NameError("ERROR: The number of bifurcations is very large ({0}). Check flag \"-override_bifurcation_limit\" to override, and prepare to be patient.".format(len(alignments.keys())))

	if code == '0':
		alignments = {}
		alignments['0'] = []
		scores = {}
		scores['0'] = 0
		i, j = NWM.shape[0]-1, NWM.shape[1]-1
	else:
		i = end_row
		j = end_column

#	print(code, i, j)
	alignments[code].append((i,j))

	branched = False
	while i>0 and j>0 and branched == False:
		group, value = NWM[i, j].best_group_subpath()
#		print(i, j, group, value)
		scores[code] += value
		if (len(group) > 1):
#			print('branch!', code, i, j)
			for branch in group:
				ii, jj = i, j
				if branch == 'D':
					ii = i-1
					jj = j-1
				elif branch == 'H':
					jj = j-1
				elif branch == 'V':
					ii = i-1

				newcode = code + str(branch)
				alignments[newcode] = []
				scores[newcode] =  scores[code]
				alignments[newcode].extend(alignments[code][:])
				print('branching now!', code, newcode, ii, jj)
				alignments, scores = traceback(sequence_1, sequence_2, NWM, alignments=alignments, scores=scores, code=newcode, end_row=ii, end_column=jj)
			branched = True
		else:
			if group[0] == 'D':
				i -= 1
				j -= 1
			elif group[0] == 'H':
				j -= 1
			elif group[0] == 'V':
				i -= 1
			alignments[code].append((i,j))
		
	if i>0 and j==0 and branched == False:
		while i > 0:
			i -= 1
			alignments[code].append((i,j))
			scores[code] += NWM[i, j].paths[NWM[i, j].path_names['V>V']]
	elif i==0 and j>0 and branched == False:
		while j > 0:
			j -= 1
			alignments[code].append((i,j))
			scores[code] += NWM[i, j].paths[NWM[i, j].path_names['H>H']]
#	else:
#		raise NameError("Fuck")

#	if branched == False:
#		alignments[code].append((i,j))

	return alignments, scores

def limits(potential_0, seq_1, seq_2, windowsize=21):
	d = []
	for na1 in range(len(seq_1[1])-windowsize):
		for na2 in range(len(seq_2[1])-windowsize):
			pot = 0.0
			for k in range(windowsize):
				pot += potential_0[na1+k, na2+k]['D']
			d.append(pot)

	dnp = np.array(d)
	fig = plt.figure(10, (15,6))
	ax = plt.subplot(111)
	M = np.amax(dnp)
	m = np.amin(dnp)
#	print(d,M,m)
#	fig, ax = plt.subplot(111)
	hist, bins = np.histogram(dnp, bins=np.arange(int(m),int(M)+1), normed=True)
	width = 0.8 * (bins[1] - bins[0])
	center = (bins[:-1] + bins[1:]) / 2
#	ax.bar(center, hist, align='center', width=width)

	return m, M, hist, bins
	
def relax_anchors(anchors, anch_boundaries, potential, sequence_1, sequence_2):
	print("Relax anchors")
	new_anchors = collections.OrderedDict()
	new_anch_boundaries = collections.OrderedDict()
	for a in anchors:
		coords = (a[0], a[1])
		maxvalue = -1000000000000
		for i in range(anch_boundaries[a][0][0], anch_boundaries[a][0][1]+1):
			for j in range(anch_boundaries[a][1][0], anch_boundaries[a][1][1]+1):
				p = potential_area(potential, sequence_1, sequence_2, (i, j), radius=3)
				if p > maxvalue:
					coords = (i, j)
					maxvalue = p
		new_anchors[coords] = anchors[a]
		new_anch_boundaries[coords] = anch_boundaries[a]
		print("anchor {0} shifted to {1}".format(a, coords))

	return new_anchors, new_anch_boundaries

def potential_area(potential, sequence_1, sequence_2, center, radius=3, weight_function='uniform'):

	def uniform(center, radius):
		return np.ones((radius*2+1, radius*2+1))

	def pyramid(center,radius):
		y = np.zeros((radius*2+1, radius*2+1))
		for i in range(center[0]-radius, center[0]+radius+1):
			for j in range(center[1]-radius, center[1]+radius+1):
				if abs(center - i) > (radius-1)/2 or abs(center - j) > (radius-1)/2:
					y[i,j] = 0
				else:
					y[i,j] = ((radius+1)/2 - abs(center - max(i,j))) / ((radius + 1)/2)
		return y
	
	allowed_weight_functions = {
	                             'uniform' : uniform,
	                             'pyramid' : pyramid
	                           }

	if weight_function not in allowed_weight_functions:
		raise NameError("ERROR: weight function {0} is not implemented".format(weight_function))

	potential_area = np.zeros((radius*2+1, radius*2+1))
	weights = allowed_weight_functions[weight_function](center, radius)
	for i in range(radius*2+1):
		ii = center[0]-radius+i
		for j in range(radius*2+1):
			jj = center[1]-radius+j
			potential_area[i,j] = potential[ii,jj]['D']

	return np.average(potential_area,weights=weights)

def enhance_diagonal(potential, anchors, anch_boundaries, weight):
	for a in anchors:
		i, j = a[0], a[1]
		while i >= anch_boundaries[a][0][0] and j >= anch_boundaries[a][1][0]:
			potential[i,j]['H'] -= weight
			potential[i,j]['V'] -= weight
			i -= 1
			j -= 1
		while i <= anch_boundaries[a][0][1] and j <= anch_boundaries[a][1][1]:
			potential[i,j]['H'] -= weight
			potential[i,j]['V'] -= weight
			i += 1
			j += 1
	return potential

def parse_inputs(parsed):
	RE = RegExp()

	# Parsing 1st sequence
	f = open(parsed.input_sequence_1[0], 'r')
	text = f.read().split('\n')
	for nl in range(len(text)):
		if RE.headerline.match(text[nl]):
			if len(text[nl]) > 1:
				seqname_1 = text[nl][1:]
			else:
				seqname_1 = 'sequence 1'
		elif not RE.emptyline.match(text[nl]):
			sequence_1 = text[nl]
			break

	# Parsing 2nd sequence
	f = open(parsed.input_sequence_2[0], 'r')
	text = f.read().split('\n')
	for nl in range(len(text)):
		if RE.headerline.match(text[nl]):
			if len(text[nl]) > 1:
				seqname_2 = text[nl][1:]
			else:
				seqname_2 = 'sequence 2'
		elif not RE.emptyline.match(text[nl]):
			sequence_2 = text[nl]
			break

	# Reading gaps
	gap = Gaps()
	if (parsed.above_threshold_gap_opening_penalty or
	    parsed.above_threshold_gap_extension_penalty or
	    parsed.below_threshold_gap_opening_penalty or
	    parsed.below_threshold_gap_extension_penalty) and ((not (parsed.above_threshold_gap_opening_penalty and
	    parsed.above_threshold_gap_extension_penalty and
	    parsed.below_threshold_gap_opening_penalty and
	    parsed.below_threshold_gap_extension_penalty)) or not parsed.threshold):
		raise NameError("ERROR: thresholding includes putting six gap-penalty flags and the -threshold flag")
	elif parsed.above_threshold_gap_opening_penalty:
		gap.gap['OpAbove'] = -float(parsed.above_threshold_gap_opening_penalty)
		gap.gap['OpBelow'] = -float(parsed.below_threshold_gap_opening_penalty)
		gap.gap['ExtAbove'] = -float(parsed.above_threshold_gap_extension_penalty)
		gap.gap['ExtBelow'] = -float(parsed.below_threshold_gap_extension_penalty)

	if not (parsed.termini_gap_opening_penalty and
	       	parsed.termini_gap_extension_penalty):
		raise NameError("ERROR: you have to put the two termini gap penalty flags")
	else:
		gap.gap['TermOp'] = -float(parsed.termini_gap_opening_penalty[0])
		gap.gap['TermExt'] = -float(parsed.termini_gap_extension_penalty[0])

	if (parsed.gap_opening_penalty or parsed.gap_extension_penalty) and not (parsed.gap_opening_penalty and parsed.gap_extension_penalty):
		raise NameError("ERROR: you should put at least four gap-penalty flags")
	elif parsed.gap_opening_penalty:
		gap.gap['Op'] = -float(parsed.gap_opening_penalty)
		gap.gap['Ext'] = -float(parsed.gap_extension_penalty)

	# Parsing similarity score file
	f = open(parsed.simscorefile[0], 'r')
	text = f.read().split('\n')
	nr = 0
	record_fields = []
	for line in text:
		if not RE.emptyline.match(line):
			fields = line.split()
			record_fields.append({})
			for nf in range(0, len(fields), 2):
				record_fields[nr][fields[nf]] = fields[nf+1]
			nr = nr + 1

	# Check option consistency
	#  Optimization
	if parsed.optimization:
		found_it = False
		for nr in range(len(record_fields)):
			if record_fields[nr]['type:'] == 'Anchors':
				found_it = True
		if found_it == False:
			raise NameError("ERROR: anchor weights optimization was enabled, but simscore file does not contain an \"Anchors\" restraint")

	return (seqname_1, sequence_1), (seqname_2, sequence_2), record_fields, gap

def build_potential(parsed, seq_1, seq_2, record_fields):

	defined_types = {
	                   'SequenceSimilarity' : associate_seqsim,
	                   'ScaleSimilarity' : associate_scale,
	                   'UniversalProfileSimilarity' : associate_profile,
	                   'PositionSpecificSimilarity' : associate_PSSM,
	                   'Anchors' : impose_anchors
	}


	sequence_1, sequence_2 = seq_1[1], seq_2[1]
	# Initialization: potential, profile
	el = Element()
	potential = np.zeros((len(sequence_1), len(sequence_2)), dtype=[('D', 'f'),('H', 'f'), ('V', 'f')])
	potential_0 = np.zeros((len(sequence_1), len(sequence_2)), dtype=[('D', 'f'),('H', 'f'), ('V', 'f')])
	profile = []

	# Build potential, profile
	for nr in range(len(record_fields)):
		if 'weight:' not in record_fields[nr]:
			raise NameError("ERROR: Record {0} does not have a \"weight:\" keyword".format(nr))
		if 'type:' in record_fields[nr]:
			inputtype = record_fields[nr]['type:']
		else:
			raise NameError("ERROR: Record {0} does not have a \"type:\" keyword".format(nr))
		if inputtype not in defined_types:
			raise NameError("ERROR: Record {0} has type = {1}, which is not one of the defined types".format(nr, inputtype))
		if inputtype == 'Anchors':
			potential = defined_types[inputtype](potential, seq_1, seq_2, record_fields[nr], nr)
			continue
		if nr == 0 and parsed.threshold:
			if 'Profile' in inputtype:
				potential_0, profile = defined_types[inputtype](potential_0, seq_1, seq_2, record_fields[nr], nr, threshold=float(parsed.threshold))
				continue
			else:
				raise NameError("ERROR: threshold is enabled, but record 0 is not a profile")
		potential_0 = defined_types[inputtype](potential_0, seq_1, seq_2, record_fields[nr], nr)

	# potential is the total one and goes into the NW, potential_0 is the potential without anchors, to optimize their weight.
	# Write potential
	f = open(str(parsed.outputfile), 'w')
	for na1 in range(len(sequence_1)):
		for na2 in range(len(sequence_2)):
			potential[na1,na2]['D'] += potential_0[na1,na2]['D']
			f.write("{0}  {1}  {2}\n".format(na1, na2, potential_0[na1,na2]['D']))

	# Initialization empty profile
	if not profile:
		for i in range(len(sequence_1)):
			profile.append([])
			for j in range(len(sequence_2)):
				profile[i].append('')

	anchors = collections.OrderedDict()

	for nr in range(len(record_fields)):
		if record_fields[nr]['type:'] == 'Anchors':
			anchors_weight = float(record_fields[nr]['weight:'])
			anchors, anch_boundaries = defined_types['Anchors'](None, (None, None), (None, None), record_fields[nr], nr, readonly=True)

	return potential, potential_0, profile, anchors, anch_boundaries, anchors_weight


def alignment(potential, profile, seq_1, seq_2):
	# Needleman-Wunsch
	t0 = time.time()
	print("Needleman-Wunsching...")
	NWM = needleman_wunsch(sequence_1, sequence_2, potential, profile, gap)
	print("\t{0} seconds".format(time.time() - t0))

	# Traceback
	t0 = time.time()
	print("Backtracing...")
	alignments, scores = traceback(sequence_1, sequence_2, NWM)
	# Choose only the alignments that reach (0,0)
	full_alignments = []
	for p in alignments:
		if alignments[p][-1] == (0,0):
			full_alignments.append((scores[p], list(reversed(alignments[p]))))
	print("\t{0} seconds".format(time.time() - t0))

	full_alignments = sorted(full_alignments, key=lambda fa: -fa[0])
	best_aln = full_alignments[0][1]

	print("\tnumber of alignments: ", len(full_alignments))

	return full_alignments


def optimization(potential, potential_0, profile, seq_1, seq_2, anchors, anch_boundaries, anch_weight):

	ncycles = 10

	def call_correct_pl(seqname_1, seqname_2, aligned_seqs):
		f = open('prov.tmp', 'w')
		f.write(">{0}:\n{1}\n>{2}:\n{3}\n".format(seqname_1, ''.join(str(p) for p in aligned_seqs[0]), seqname_2, ''.join(str(p) for p in aligned_seqs[1])))
		f.flush()
		os.fsync(f.fileno())
		f.close()
		subprocess.call(['perl', 'correct.pl', '1ZOY_D_2H88_C.fa', 'prov.tmp'])

	seqname_1, seqname_2 = seq_1[0], seq_2[0]
	sequence_1, sequence_2 = seq_1[1], seq_2[1]

	free_alignments = alignment(potential_0, profile, seq_1, seq_2)

	best_score, best_aln = free_alignments[0][0], free_alignments[0][1]
	aligned_seqs, matches, signatures = to_fasta((sequence_1, sequence_2), best_aln, anchors=anchors)
	print("Cycle 0: unrestrained")
	print("Score {0}:\n> {1}:\n{2}\n{3}\n> {4}:\n{5}\n{6}\n".format(best_score, seqname_1, ''.join(str(p) for p in aligned_seqs[0]), signatures[0], seqname_2, ''.join(str(p) for p in aligned_seqs[1]), signatures[1]))
	call_correct_pl(seqname_1, seqname_2, aligned_seqs)

#	potential_0 = modify_penalties(potential_0, profile, seq_1, seq_2, anchors)

	for i in range(ncycles):
		anchors = anchors_weights(anchors, anch_boundaries, seq_1, seq_2, best_aln)

		el = Element()
		for na1 in range(len(sequence_1)):
			for na2 in range(len(sequence_2)):
				for a in anchors:
					if (na1 < a[0] and na2 >= a[1]) or (na1 >= a[0] and na2 < a[1]):
						for d in el.directions:
							potential[na1,na2][d] = potential_0[na1,na2][d] - anch_weight * anchors[a]	

		restrained_alignments = alignment(potential, profile, seq_1, seq_2)
		best_score, best_aln = restrained_alignments[0][0], restrained_alignments[0][1]
		aligned_seqs, matches, signatures = to_fasta((sequence_1, sequence_2), best_aln, anchors=anchors)
		print("Cycle {0}: restrained.\n\tG-Weight: {1}".format(i+1, anch_weight))
		for a in anchors:
			print("\tanchor: {0}\t\tR-Weight: {1}".format(a, anchors[a]))
		print("Score {0}:\n> {1}:\n{2}\n{3}\n> {4}:\n{5}\n{6}\n".format(best_score, seqname_1, ''.join(str(p) for p in aligned_seqs[0]), signatures[0], seqname_2, ''.join(str(p) for p in aligned_seqs[1]), signatures[1]))
		call_correct_pl(seqname_1, seqname_2, aligned_seqs)


def modify_penalties(potential, profile, seq_1, seq_2, anchors):

	enh_factor = 0.2
	def DHV_count(block, potential, enh_factor, best_d):
		count = {
		          'D' : 0,
		          'H' : 0,
		          'V' : 0
		        }

		prev = {
		         'D' : lambda x, y: y-1,
		         'H' : lambda x, y: y if x==0 else y-1,
		         'V' : lambda x, y: y-1 if x==0 else y,
		       }
		
		for ii in range(blocks[i][0][0], blocks[i][1][0]): 
			for jj in range(blocks[i][0][1], blocks[i][1][1]):
				best_d[ii, jj]['enh'] = 0
				best_d[ii, jj]['val'] = -1000000000000
				for d in ['D', 'V', 'H']:
					p_ii, p_jj = prev[d](0, ii), prev[d](1, jj)
					if best_d[ii, jj]['val'] <  potential[ii, jj][d]:
						best_d[ii, jj]['dir'] = d
						best_d[ii, jj]['val'] = potential[ii, jj][d]
						if p_ii >= 0 and p_jj >= 0 and best_d[p_ii, p_jj]['dir'] == d:
							best_d[ii, jj]['enh'] = best_d[p_ii, p_jj]['enh'] + enh_factor
		for ii in range(blocks[i][0][0], blocks[i][1][0]):
			for jj in range(blocks[i][0][1], blocks[i][1][1]):
				count[best_d[ii, jj]['dir']] += 1 + count[best_d[ii, jj]['enh']]
		norm = 0
		for d in ['D', 'V', 'H']:
			norm += count[d]
		for d in ['D', 'V', 'H']:
			count[d] /= norm

		return count, best_d

	blocks = []
	blocks.append([(0,0)])
	na = 0
	for a in anchors:
		blocks[na].append((a[0],a[1]))
		na += 1
		blocks.append([(a[0],a[1])])
	blocks[na].append((len(seq_1[1]),len(seq_2[1])))

	free_alignments = alignment(potential, profile, seq_1, seq_2)
	best_score, best_aln = free_alignments[0][0], free_alignments[0][1]

	DHV_counts = []
	best_d = np.empty(potential.size, dtype=[('dir', np.str_, 1), ('val', np.float64), ('enh', np.float64)])
	for i in range(len(blocks)):
		count, best_d = DHV_count(blocks[i], potential, enh_factor, best_d)
		DHV_counts.append(count)

def anchors_weights(anchors, anch_boundaries, seq_1, seq_2, best_aln):
	windowsize = 21
	lmargin, rmargin = 0, 0

	seqname_1, seqname_2 = seq_1[0], seq_2[0]
	sequence_1, sequence_2 = seq_1[1], seq_2[1]

	if windowsize % 2 == 0:
		raise NameError("ERROR: window size must be an odd number\n")
	elif windowsize > min(len(sequence_1), len(sequence_2)):
		raise NameError("ERROR: window size must be less than sequence lengths\n")
	else:
		mg = int(windowsize/2)
		lmargin, rmargin = mg, mg

	new_anchors = collections.OrderedDict()
	for a in anchors:
		magic_number = a[0] + a[1]
		for couple in best_aln:
			if couple[0] + couple[1] == magic_number or couple[0] + couple[1] == magic_number - 1:
				middle_point = [couple[0], couple[1]]

		check_extrema = False
		while check_extrema == False and lmargin > int(windowsize/4) and lmargin > int(windowsize/4) and rmargin > int(windowsize/4) and rmargin > int(windowsize/4):
			if middle_point[0] - lmargin >= 0 and middle_point[0] + rmargin < len(sequence_1) and middle_point[1] - lmargin >= 0 and middle_point[1] + rmargin < len(sequence_2):
				check_extrema = True
			else:
				lmargin -= 1
				rmargin-= 1

		while check_extrema == False:
			if middle_point[0] - lmargin >= 0 and middle_point[0] + rmargin < len(sequence_1) and middle_point[1] - lmargin >= 0 and middle_point[1] + rmargin < len(sequence_2):
				check_extrema = True
			else:
				if middle_point[0] - lmargin < 0:
					middle_point[0] += 1
				elif middle_point[0] + rmargin >= len(sequence_1):
					middle_point[0] -= 1
				if middle_point[1] - lmargin < 0:
					middle_point[1] += 1
				elif middle_point[1] + rmargin >= len(sequence_2):
					middle_point[1] -= 1

		pot = 0.0
		for k in range(rmargin+lmargin+1):
			i = middle_point[0] - lmargin + k
			j = middle_point[1] - lmargin + k
			pot += potential_0[i,j]['D']

		min_score, max_score, score_histogram, score_h_edges = limits(potential_0, (seqname_1, sequence_1), (seqname_2, sequence_2), windowsize=rmargin+lmargin+1)
		new_anchors[a] = np.cumsum(score_histogram)[np.amin(score_h_edges[score_h_edges > pot])]

	for a in new_anchors:
		anchors[a] = new_anchors[a]
#	print(anchors)

	return anchors



t00 = time.time()
t0 = time.time()
print("Initializing...", end=" ")

seq_1, seq_2, record_fields, gap = parse_inputs(parsed)
sequence_1, sequence_2 = seq_1[1], seq_2[1]

potential, potential_0, profile, anchors, anch_boundaries, anch_weight = build_potential(parsed, seq_1, seq_2, record_fields)
	
print("{0} seconds".format(time.time() - t0))

anchors, anch_boundaries = relax_anchors(anchors, anch_boundaries, potential_0, sequence_1, sequence_2)

optimization(potential, potential_0, profile, seq_1, seq_2, anchors, anch_boundaries, anch_weight)

#f = open('alignments.txt', 'w')
#for i in range(len(full_alignments)):
#	aligned_seq_1, aligned_seq_2, matches = to_fasta(sequence_1, sequence_2, full_alignments[i][1])
#	f.write("Score {0}:\n> {1}:\n{2}\n> {3}:\n{4}\n".format(full_alignments[i][0], seqname_1, ''.join(str(p) for p in aligned_seq_1), seqname_2, ''.join(str(p) for p in aligned_seq_2)))

tf = time.time()
print("{0} seconds".format(tf - t0))
print("Total time: {0} seconds".format(tf - t00))
print("Done!")
