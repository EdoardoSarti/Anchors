#!/Users/Orfeo/anaconda3/bin/python3

import argparse, re, collections, time
import numpy as np
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument('-i1', '--input_sequence_1', nargs=1)
parser.add_argument('-i2', '--input_sequence_2', nargs=1)
parser.add_argument('-s', '--simscorefile', nargs=1)
parser.add_argument('-o', '--outputfile', nargs='?', default='NW_score.txt')
parser.add_argument('-oo', '--outputfile_optimization', nargs='?', default='alignments.txt')
parser.add_argument('-thr', '--threshold', nargs='?')
parser.add_argument('-optimize', dest='optimization', action='store_true')
parser.set_defaults(optimization=False)
parsed = parser.parse_args()

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
		             'Op' : -15.36,
		             'Ext' : -0.88,
		             'OpAbove' : 0.0,
		             'OpBelow' : 0.0,
		             'ExtAbove' : 0.0,
		             'ExtBelow' : 0.0,
		             'TermOp' : -1.69,
                             'TermExt' : -0.25
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
	if (plottype == 'substitution_matrix'):
		fig = plt.figure(nr, (15,6))
	else:
		fig = plt.figure(nr, (10,6))
	if (plottype == 'scale' or plottype == 'substitution_matrix'):
		fig.suptitle(plottype + ': ' + record_fields['file:'].split('/')[-1])
	else:
		fig.suptitle(plottype + ': similarity score file line ' + str(nr))

	if (plottype == 'substitution_matrix'):
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
	elif (plottype == 'substitution_matrix'):
		ax.pcolor(instr_1[0], cmap=plt.cm.Blues)
	else:
		ax.plot(instr_1[0], instr_1[1])
	ax.set_xlim(0, len(seqs[0][1]))
	
	if (plottype == 'substitution_matrix'):
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
	elif (plottype == 'substitution_matrix'):
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
			potential[((na1, seq_1[1][na1]), (na2, seq_2[1][na2]), 'D')] -= float(record_fields['weight:']) * abs(winav_1[na1] - winav_2[na2])

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
			potential[((na1, seq_1[1][na1]), (na2, seq_2[1][na2]), 'D')] -= float(record_fields['weight:']) * matrix[(seq_1[1][na1], seq_2[1][na2])] 

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
			potential[((na1, seq_1[1][na1]), (na2, seq_2[1][na2]), 'D')] -= float(record_fields['weight:']) * abs(profile_1[na1] - profile_2[na2]) 
	
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
	pass

def impose_anchors(potential, seq_1, seq_2, record_fields, nr):
	RE = RegExp()
	check_attributes(record_fields, nr, 'file:')

	f = open(record_fields['file:'], 'r')
	text = f.read().split('\n')
	f.close()
	anchors = {}
	for nl in range(len(text)):
		fields = text[nl].split()
		if nl == 0:
			n_anchors = int(fields[0])
		elif nl <= n_anchors:
			if len(fields) == 2:
				anchors[(int(fields[0]), int(fields[1]))] = None
			if len(fields) == 3:
				anchors[(int(fields[0]), int(fields[1]))] = float(fields[2])

	el = Element()
	for na1 in range(len(seq_1[1])):
		for na2 in range(len(seq_2[1])):
			for a in anchors:
				if (na1 < a[0] and na2 >= a[1]) or (na1 >= a[0] and na2 < a[1]):
					for d in el.directions:
						potential[((na1, seq_1[1][na1]), (na2, seq_2[1][na2]), d)] -= anchors[a]
	return potential

def to_fasta(sequence_1, sequence_2, alignment):
	if alignment[0] != (0, 0):
		raise NameError("ERROR: The alignment to be translated must start at (0,0)")

	aligned_seq_1 = []
	aligned_seq_2 = []
	matches = []
	i_old, j_old = alignment[0][0], alignment[0][1]
	for k in range(1, len(alignment)):
		i, j = alignment[k][0], alignment[k][1]
		if i == i_old:
			aligned_seq_1.append('-')
			aligned_seq_2.append(sequence_2[j])
		elif j == j_old:
			aligned_seq_1.append(sequence_1[i])
			aligned_seq_2.append('-')
		else:
			aligned_seq_1.append(sequence_1[i])
			aligned_seq_2.append(sequence_2[j])
		i_old, j_old = i, j
		if aligned_seq_1[-1] == aligned_seq_2[-1]:
			matches.append(':')
		else:
			matches.append(' ')
	return aligned_seq_1, aligned_seq_2, matches

def needleman_wunsch(sequence_1, sequence_2, potential, profile):
	NWM = []
	gap = Gaps()

	# First line
	NWline0 = []
	el = Element(0)
	NWline0.append(el)
	for j in range(1, len(sequence_2)+1):
		el = Element()
		el.add_path_value(gap.gap['TermOp'] + (j-1)*gap.gap['TermExt'] + sum([ potential[((0, sequence_1[0]), (k-1, sequence_2[k-1]), 'H')] for k in range(1, j+1)]), path_name='H>H')
		NWline0.append(el)
	NWM.append(NWline0)

	# Main loop
	for i in range(1, len(sequence_1)):
		NWline = []
		el = Element()

		# First column elements
		el.add_path_value(gap.gap['TermOp'] + (i-1)*gap.gap['TermExt'] + sum([ potential[((k-1, sequence_1[k-1]), (0, sequence_2[0]), 'V')] for k in range(1, i+1)]), path_name='V>V')
		NWline.append(el)

		# Core elements
		for j in range(1, len(sequence_2)):
			el = Element(0)

			# Diagonal values
			d2, el_past = 'D', NWM[i-1][j-1]
			for d1 in el.directions:
				path_name = d1 + '>' + d2
#				print(i,j,path_name, el_past.best_subpath_value(d1), potential[((i-1, sequence_1[i-1]), (j-1, sequence_2[j-1]), d2)])
				value = el_past.best_subpath_value(d1) + potential[((i-1, sequence_1[i-1]), (j-1, sequence_2[j-1]), d2)]
				el.add_path_value(value, path_name=path_name)

			# Horizontal and vertical values
			for d2, el_past in [('H', NWline[j-1]), ('V', NWM[i-1][j])]:
				# When d1 and d2 are different...
				for d1 in list(set(el.directions) - set(d2)):
					path_name = d1 + '>' + d2
#					print(i,j,path_name, el_past.best_subpath_value(d1), gap.gap[ 'Op' + profile[i-1][j-1] ], potential[((i-1, sequence_1[i-1]), (j-1, sequence_2[j-1]), d2)])
					value = el_past.best_subpath_value(d1) + gap.gap[ 'Op' + profile[i-1][j-1] ] + potential[((i-1, sequence_1[i-1]), (j-1, sequence_2[j-1]), d2)]
					el.add_path_value(value, path_name=path_name)
				# ...and when they are equal.
				d1 = d2
				path_name = d1 + '>' + d2
#				print(i,j,path_name, el_past.best_subpath_value(d1), gap.gap[ 'Ext' + profile[i-1][j-1] ], potential[((i-1, sequence_1[i-1]), (j-1, sequence_2[j-1]), d2)])
				value = el_past.best_subpath_value(d1) + gap.gap[ 'Ext' + profile[i-1][j-1] ] + potential[((i-1, sequence_1[i-1]), (j-1, sequence_2[j-1]), d2)]
				el.add_path_value(value, path_name=path_name)
	
			NWline.append(el)
	
		NWM.append(NWline)
	
	# Last line
	NWlineT = []
	el = Element()
	i = len(sequence_1)
	el.add_path_value(gap.gap['TermOp'] + (i-1)*gap.gap['TermExt'] + sum([ potential[((k-1, sequence_1[k-1]), (0, sequence_2[0]), 'H')] for k in range(1, i+1)]), path_name='V>V')
	NWlineT.append(el)
	for j in range(1, len(sequence_2)):
		el = Element(0)

		# Diagonal values
		d2, el_past = 'D', NWM[i-1][j-1]
		for d1 in el.directions:
			path_name = d1 + '>' + d2
			value = el_past.best_subpath_value(d1) + potential[((i-1, sequence_1[i-1]), (j-1, sequence_2[j-1]), d2)]
			el.add_path_value(value, path_name=path_name)

		# Horizontal values (terminal gaps)
		d2, el_past = 'H', NWlineT[j-1]
		for d1 in list(set(el.directions) - set(d2)):
			path_name = d1 + '>' + d2
			value = el_past.best_subpath_value(d1) + gap.gap[ 'TermOp' ] + potential[((i-1, sequence_1[i-1]), (j-1, sequence_2[j-1]), d2)]
			el.add_path_value(value, path_name=path_name)
		d1 = d2
		path_name = d1 + '>' + d2
		value = el_past.best_subpath_value(d1) + gap.gap[ 'TermExt' ] + potential[((i-1, sequence_1[i-1]), (j-1, sequence_2[j-1]), d2)]
		el.add_path_value(value, path_name=path_name)

		# Vertical values (normal gaps)
		d2, el_past = 'V', NWM[i-1][j]
		for d1 in list(set(el.directions) - set(d2)):
			path_name = d1 + '>' + d2
			value = el_past.best_subpath_value(d1) + gap.gap[ 'Op' + profile[i-1][j-1] ] + potential[((i-1, sequence_1[i-1]), (j-1, sequence_2[j-1]), d2)]
			el.add_path_value(value, path_name=path_name)
		d1 = d2
		path_name = d1 + '>' + d2
		value = el_past.best_subpath_value(d1) + gap.gap[ 'Ext' + profile[i-1][j-1] ] + potential[((i-1, sequence_1[i-1]), (j-1, sequence_2[j-1]), d2)]
		el.add_path_value(value, path_name=path_name)

		NWlineT.append(el)
	NWM.append(NWlineT)

	# Last row (last element of each line)
	j = len(sequence_2)
	for i in range(1, len(sequence_1)): 
		el = Element(0)

		# Diagonal values
		d2, el_past = 'D', NWM[i-1][j-1]
		for d1 in el.directions:
			path_name = d1 + '>' + d2
			value = el_past.best_subpath_value(d1) + potential[((i-1, sequence_1[i-1]), (j-1, sequence_2[j-1]), d2)]
			el.add_path_value(value, path_name=path_name)

		# Horizontal values (normal gaps)
		d2, el_past = 'H', NWM[i][j-1]
		for d1 in list(set(el.directions) - set(d2)):
			path_name = d1 + '>' + d2
			value = el_past.best_subpath_value(d1) + gap.gap[ 'Op' + profile[i-1][j-1] ] + potential[((i-1, sequence_1[i-1]), (j-1, sequence_2[j-1]), d2)]
			el.add_path_value(value, path_name=path_name)
		d1 = d2
		path_name = d1 + '>' + d2
		value = el_past.best_subpath_value(d1) + gap.gap[ 'Ext' + profile[i-1][j-1] ] + potential[((i-1, sequence_1[i-1]), (j-1, sequence_2[j-1]), d2)]
		el.add_path_value(value, path_name=path_name)

		# Vertical values (terminal gaps)
		d2, el_past = 'V', NWM[i-1][j]
		for d1 in list(set(el.directions) - set(d2)):
			path_name = d1 + '>' + d2
			value = el_past.best_subpath_value(d1) + gap.gap[ 'TermOp' ] + potential[((i-1, sequence_1[i-1]), (j-1, sequence_2[j-1]), d2)]
			el.add_path_value(value, path_name=path_name)
		d1 = d2
		path_name = d1 + '>' + d2
		value = el_past.best_subpath_value(d1) + gap.gap[ 'TermExt' ] + potential[((i-1, sequence_1[i-1]), (j-1, sequence_2[j-1]), d2)]
		el.add_path_value(value, path_name=path_name)

		NWM[i].append(el)

	# Last element (terminal gaps either way)
	i, j = len(sequence_1), len(sequence_2)
	el = Element(0)

	# Diagonal values
	d2, el_past = 'D', NWM[i-1][j-1]
	for d1 in el.directions:
		path_name = d1 + '>' + d2
		value = el_past.best_subpath_value(d1) + potential[((i-1, sequence_1[i-1]), (j-1, sequence_2[j-1]), d2)]
		el.add_path_value(value, path_name=path_name)
	
	# Horizontal and vertical values (terminal gaps)
	for d2, el_past in [('H', NWM[i][j-1]), ('V', NWM[i-1][j])]:
		# When d1 and d2 are diffrent...
		for d1 in list(set(el.directions) - set(d2)):
			path_name = d1 + '>' + d2
			value = el_past.best_subpath_value(d1) + gap.gap[ 'TermOp' ] + potential[((i-1, sequence_1[i-1]), (j-1, sequence_2[j-1]), d2)]
			el.add_path_value(value, path_name=path_name)
		# ...and when they are equal
		d1 = d2
		path_name = d1 + '>' + d2
		value = el_past.best_subpath_value(d1) + gap.gap[ 'TermExt' ] + potential[((i-1, sequence_1[i-1]), (j-1, sequence_2[j-1]), d2)]
		el.add_path_value(value, path_name=path_name)

	NWM[i].append(el)

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
		i = len(sequence_1)-1
		j = len(sequence_2)-1
	else:
		i = end_row
		j = end_column

	alignments[code].append((i,j))

	while i>0 and j>0:
		group, value = NWM[i][j].best_group_subpath()
#		print(i, j, group, value)
		scores[code] += value
		if (len(group) > 1):
#			print('branch', i, j, group, value)
			for branch in range(len(group)):
				if group[branch] == 'D':
					i -= 1
					j -= 1
				elif group[branch] == 'H':
					j -= 1
				elif group[branch] == 'V':
					i -= 1

				newcode = code + str(branch)
				scores[newcode] =  scores[code]
				alignments[newcode] = alignments[code]
				alignments, scores = traceback(sequence_1, sequence_2, NWM, alignments=alignments, scores=scores, code=newcode, end_row=i, end_column=j)
		else:
			if group[0] == 'D':
				i -= 1
				j -= 1
			elif group[0] == 'H':
				j -= 1
			elif group[0] == 'V':
				i -= 1
			alignments[code].append((i,j))
		
	if i>0 and j==0:
		while i > 0:
			alignments[code].append((i,j))
			scores[code] += NWM[i][j].paths[NWM[i][j].path_names['V>V']]
			i -= 1
	elif i==0 and j>0:
		while j > 0:
			alignments[code].append((i,j))
			scores[code] += NWM[i][j].paths[NWM[i][j].path_names['H>H']]
			j -= 1
	else:
		raise NameError("Fuck")

	alignments[code].append((i,j))

	return alignments, scores

def optimization():
	pass




defined_types = {
                   'SequenceSimilarity' : associate_seqsim,
                   'ScaleSimilarity' : associate_scale,
                   'UniversalProfileSimilarity' : associate_profile,
                   'PositionSpecificSimilarity' : associate_PSSM,
                   'Anchors' : impose_anchors
}


t00 = time.time()
t0 = time.time()
print("Initializing...", end=" ")

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

# Initialization: potential, profile
el = Element()
potential = {}
for na1 in range(len(sequence_1)):
	for na2 in range(len(sequence_2)):
		for d in el.directions:
			potential[((na1, sequence_1[na1]), (na2, sequence_2[na2]), d)] = 0.0
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
	if nr == 0 and parsed.threshold:
		if inputtype == 'Profile':
			potential, profile = defined_types[inputtype](potential, (seqname_1, sequence_1), (seqname_2, sequence_2), record_fields[nr], nr, threshold=float(parsed.threshold))
			continue
		else:
			raise NameError("ERROR: threshold is enabled, but record 0 is not a profile")
	potential = defined_types[inputtype](potential, (seqname_1, sequence_1), (seqname_2, sequence_2), record_fields[nr], nr)

	'''
	f = open('pot_'+str(nr)+'_'+inputtype+'.txt', 'w')
	for na1 in range(len(sequence_1)):
		for na2 in range(len(sequence_2)):
			f.write("{0}  {1}  {2}\n".format(na1, na2, potential[((na1, sequence_1[na1]), (na2, sequence_2[na2]), 'D')]))
	'''

# Write potential
f = open(str(parsed.outputfile), 'w')
for na1 in range(len(sequence_1)):
	for na2 in range(len(sequence_2)):
		f.write("{0}  {1}  {2}\n".format(na1, na2, potential[((na1, sequence_1[na1]), (na2, sequence_2[na2]), 'D')]))

# Initialization empty profile
if not profile:
	for i in range(len(sequence_1)):
		profile.append([])
		for j in range(len(sequence_2)):
			profile[i].append('')

print("{0} seconds".format(time.time() - t0))


# Needleman-Wunsch
t0 = time.time()
print("Needleman-Wunsching...", end=" ")
NWM = needleman_wunsch(sequence_1, sequence_2, potential, profile)
print("{0} seconds".format(time.time() - t0))


# Traceback
t0 = time.time()
print("Backtracing...", end=" ")
alignments, scores = traceback(sequence_1, sequence_2, NWM)

# Choose only the alignments that reach (0,0)
full_alignments = []
for p in alignments:
	if alignments[p][-1] == (0,0):
		full_alignments.append((scores[p], list(reversed(alignments[p]))))

print("{0} seconds".format(time.time() - t0))


# Write alignments in fasta format
t0 = time.time()
print("Writing...", end=" ")
full_alignments = sorted(full_alignments, key=lambda fa: fa[0])
f = open('alignments.txt', 'w')
for i in range(len(full_alignments)):
	aligned_seq_1, aligned_seq_2, matches = to_fasta(sequence_1, sequence_2, full_alignments[i][1])
	f.write("Score {0}:\n> {1}:\n{2}\n> {3}:\n{4}\n\n\n\n".format(full_alignments[i][0], seqname_1, ''.join(str(p) for p in aligned_seq_1), seqname_2, ''.join(str(p) for p in aligned_seq_2)))

tf = time.time()
print("{0} seconds".format(tf - t0))
print("Total time: {0} seconds".format(tf - t00))
print("Done!")
