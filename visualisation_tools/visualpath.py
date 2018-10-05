#!/Users/Orfeo/anaconda3/bin/python3
import sys, argparse, operator
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--fasta', nargs='*')
parser.add_argument('-s', '--score', nargs='*')
parser.add_argument('-cs', '--colorscore', nargs='*')
parser.add_argument('--normalize_score', dest='normalization', action='store_true')
parser.add_argument('-ev', '--every', nargs='?')
parser.add_argument('-co', '--cutoff', nargs='?')
parser.add_argument('-aa', '--anchors', nargs='?')
parser.set_defaults(every=20)
parser.set_defaults(cutoff=-1000)
parser.set_defaults(normalization=False)
parsed = parser.parse_args()

title1 = {} 
seq1 = {}
title2 = {}
seq2 = {}
chars1 = {}
chars2 = {}

if parsed.anchors:
	f = open(parsed.anchors, 'r')
	text = f.read().split('\n')
	anchors = {}
	for nl in range(len(text)):
		fields = text[nl].split()
		if nl == 0:
 			n_anchors = int(fields[0])
		elif nl <= n_anchors:
 			anchors[(int(fields[0]), int(fields[1]))] = float(fields[2])


def plot_score(sf, ax, norm, step=10, start=0, color='0.4'):
	scorefile = open(sf)
	text = scorefile.read().split('\n')
	score = {}
	for line in range(len(text)):
		if(text[line] != ''):
			if(len(text[line].split()) != 3):
				raise NameError("corrupted score file")
			prov = text[line].split()
			if float(prov[2]) > float(parsed.cutoff) or parsed.normalization:
				score[(int(prov[0])), int(prov[1])] = float(prov[2])
			else:
				score[(int(prov[0])), int(prov[1])] = float(parsed.cutoff)

	del score[(0, 0)]
	maxind = max(score.items(), key=operator.itemgetter(0))[0]
	del score[maxind]
	del score[(maxind[0]-1, maxind[1])]
	del score[(maxind[0], maxind[1]-1)]

	if norm:
		maxscore = max(score.items(), key=operator.itemgetter(1))[1]
		minscore = min(score.items(), key=operator.itemgetter(1))[1]
#		print(minscore, maxscore)
		for i in score:
			score[i] = (score[i] - minscore)/(maxscore - minscore)

	if step == 1:
		ax.scatter(list(zip(*score.keys()))[0][::100], list(zip(*score.keys()))[1][::100], list(score.values())[::100])
		"""
		provx = []
		provy = []
		provs = []
		for i, j in score:
			if score[i,j] < -1.8:
				provx.append(i)
				provy.append(j)
				provs.append(score[(i,j)])
		ax.scatter(provx, provy, provs)
		"""
		return ax, score

	coord_x = {}
	coord_y = {}
	score_x = {}
	score_y = {}

	j = 0
	while((1, j) in score):
		coord_y[j] = []
		score_y[j] = []
		j += step
	coord_y[maxind[1]] = []
	score_y[maxind[1]] = []

	i = 0
	while ((i, 1) in score):
		if i == 0 and start == 0:
			j = step
		else:
			j = start
		coord_x[i] = []
		score_x[i] = []
		while((i, j) in score):
			coord_x[i].append((i, j))
			score_x[i].append(score[(i, j)])
			coord_y[j].append((i, j))
			score_y[j].append(score[(i, j)])
			j += step
		if(i, maxind[1]) in score:
			coord_x[i].append((i, maxind[1]))
			score_x[i].append(score[(i, maxind[1])])
			coord_y[maxind[1]].append((i, maxind[1]))
			score_y[maxind[1]].append(score[(i, maxind[1])])
		i += step
	coord_x[maxind[0]] = []
	score_x[maxind[0]] = []
	j = start
	while((maxind[0], j) in score):
#		print(maxind[0], j)
		coord_x[maxind[0]].append((maxind[0], j))
		score_x[maxind[0]].append(score[(maxind[0], j)])
		coord_y[j].append((maxind[0], j))
		score_y[j].append(score[(maxind[0], j)])
		j += step

	for i in coord_x:
		px, py = zip(*coord_x[i])
		ax.plot(px, py, score_x[i], c=color)
	for j in coord_y:
#		print(j, coord_y[j])
		px, py = zip(*coord_y[j])
		ax.plot(px, py, score_y[j], c=color)

#	print(coord_x[0], score_x[0])

	return ax, score

def plot_with_m(score_filename, parsed, step=1, color=0.4):
	score_file = open(score_filename)
	text = score_file.read().split('\n')
	score_file.close()
	score = {}
	for line in range(len(text)):
		if(text[line] != ''):
			if(len(text[line].split()) != 3):
				raise NameError("corrupted score file")
			prov = text[line].split()
			if float(prov[2]) > float(parsed.cutoff) or parsed.normalization:
				score[(int(prov[0])), int(prov[1])] = float(prov[2])
			else:
				score[(int(prov[0])), int(prov[1])] = float(parsed.cutoff)

	if parsed.normalization:
		maxscore = max(score.items(), key=operator.itemgetter(1))[1]
		minscore = min(score.items(), key=operator.itemgetter(1))[1]
		for i in score:
			score[i] = (score[i] - minscore)/(maxscore - minscore)

	maxind = (max([x[0] for x in score]), max([x[1] for x in score]))
	score_mx = np.zeros((maxind[0]+1, maxind[1]+1))
	for coords in score:
		score_mx[coords[0], coords[1]] = score[coords]

	x = []
	y = []
	z = []
	s = []
	connections = []
	N = 300
	index = 0
#	score_mx = score_mx[:10,:10]
	for i in range(score_mx.shape[0]):
		x.append(np.ones(score_mx.shape[1])*i)
		y.append(np.linspace(0, score_mx.shape[1]-1, score_mx.shape[1]))
		z.append(score_mx[i])
		s.append(np.ones(score_mx.shape[1])*np.pi)
#		print(x[-1].shape, y[-1].shape, z[-1].shape)
		connections.append(np.vstack([np.arange(index, index + score_mx.shape[1] - 1.5),np.arange(index + 1, index + score_mx.shape[1] - .5)]).T)
		index += score_mx.shape[1]
	for j in range(score_mx.shape[1]):
		y.append(np.ones(score_mx.shape[0])*j)
		x.append(np.linspace(0, score_mx.shape[0]-1, score_mx.shape[0]))
		z.append(score_mx.T[j])
		s.append(np.zeros(score_mx.shape[0]))
		connections.append(np.vstack([np.arange(index, index + score_mx.shape[0] - 1.5),np.arange(index + 1, index + score_mx.shape[0] - .5)]).T)
		index += score_mx.shape[0]

	x = np.hstack(x)
	y = np.hstack(y)
	z = np.hstack(z)
	s = np.hstack(z)
	connections = np.vstack(connections)
#	print(x,y,z,connections)
#	print(x.shape, y.shape, z.shape, connections.shape)

	return x,y,z,s,connections,score_mx,index

"""
fig = plt.figure(figsize=(15,10))
ax = fig.add_subplot(111, projection='3d')
ax.view_init(elev=90., azim=0.)
"""

if parsed.score:
	"""
	if parsed.colorscore and len(parsed.score) == len(parsed.colorscore):
		scoreandcolor = zip(parsed.score, parsed.colorscore)
		for name, color in scoreandcolor:
			ax, score = plot_score(name, ax, parsed.normalization, step=int(parsed.every), color=color)
	"""
	for name in parsed.score:
		x,y,z,s,connections,score_mx,index = plot_with_m(name, parsed)
#		ax, score = plot_score(name, ax, parsed.normalization, step=int(parsed.every))
	"""
	if parsed.anchors:
		c_x = []
		c_y = []
		c_z = []
		for a in anchors:
			c_x.append(a[0])
			c_y.append(a[1])
			c_z.append(score[(a[0],a[1])])
		ax.scatter(c_x, c_y, c_z, s=40, c='k')
	"""

#if (parsed.normalization):
#	ax.auto_scale_xyz([0, max(score.items(), key=operator.itemgetter(0))[0][0]], [0, max(score.items(), key=operator.itemgetter(0))[0][1]], [-1, 2])

xp = []
yp = []
zp = []
sp = []
connectionsp = []
if parsed.fasta:
	for name in parsed.fasta:
		fastafile = open(name, 'r+')
		text = fastafile.read().split('\n')
		title1[name] = text[0]
		seq1[name] = text[1]
		title2[name] = text[2]
		seq2[name] = text[3]
	
		# Checks
		if (title1[name][0] != '>' or title2[name][0] != '>'):
			raise NameError("File {0}: \nSequences must be preceded by a title line starting with \">\"".format(name))
	
		if (len(seq1[name]) != len(seq2[name])):
			raise NameError("File {0}: \nThe two sequences must have the same length. \nSeq1: {0} \nSeq2: {1}".format(name, len(seq1[name]), len(seq2[name])))
	
		chars1[''.join(filter(str.isalpha, seq1[name]))] = name
		chars2[''.join(filter(str.isalpha, seq2[name]))] = name
		if (len(chars1.keys()) != 1):
			raise NameError("File {0}: \nSequence 1 does not coincide with other files' analogous")
		if (len(chars2.keys()) != 1):
			raise NameError("File {0}: \nSequence 2 does not coincide with other files' analogous")
	
		nseq1 = 0
		nseq2 = 0
		color = 2*score_mx.max()
		xp = [ nseq1 ]
		yp = [ nseq2 ]
		zp = [ score_mx[nseq1,nseq2] ]
		sp = [ color ]
#		print("{0}\n{1}".format(seq1[name], seq2[name]))
		for i in range(1,len(seq1[name])):
			a = seq1[name][i-1]
			b = seq2[name][i-1]
			if (a != '-'):
				nseq1 += 1
			if (b != '-'):
				nseq2 += 1
			xp.append(nseq1)
			yp.append(nseq2)
			zp.append(score_mx[xp[-1],yp[-1]])
			sp.append(color)
#		ax.plot(x, y, s, label=name, linewidth=2)

#	plt.xlabel(title1[parsed.fasta[0]][1:])
#	plt.ylabel(title2[parsed.fasta[0]][1:])
#	ax.legend()	
	x = np.hstack((x,xp))
	y = np.hstack((y,yp))
	z = np.hstack((z,zp))
	s = np.hstack((s,sp))
	connectionsp = np.vstack([np.arange(index, index + len(xp) - 1.5),np.arange(index + 1, index + len(xp) - .5)]).T
	connections = np.vstack((connections, connectionsp))

#print(connections)

from mayavi import mlab
mlab.figure(1, size=(400, 400), bgcolor=(0, 0, 0))
mlab.clf()	

# Create the points
src = mlab.pipeline.scalar_scatter(x, y, z, s)
src2 = mlab.pipeline.scalar_scatter([xp], [yp], [zp], [sp])

# Connect them
src.mlab_source.dataset.lines = connections
src.update()
src2.mlab_source.dataset.lines = connectionsp
src2.update()

# The stripper filter cleans up connected lines
#lines = mlab.pipeline.stripper(src)

# Finally, display the set of lines
mlab.pipeline.surface(src, colormap='Blues', line_width=1, opacity=.4)
mlab.pipeline.surface(src2, colormap='Reds', line_width=3, opacity=.4)
# And choose a nice view
#mlab.move(int(score_mx.shape[0]/2), int(score_mx.shape[0]/2), score_mx[int(score_mx.shape[0]/2)][int(score_mx.shape[1]/2)])
mlab.view(33.6, 106, 5.5, [int(score_mx.shape[0]/2), int(score_mx.shape[0]/2), score_mx[int(score_mx.shape[0]/2)][int(score_mx.shape[1]/2)]])
#mlab.roll(125)
mlab.show()

cmap = plt.cm.Blues
norm = plt.Normalize(score_mx.min(), score_mx.max())
rgba = cmap(norm(score_mx))

# Set the diagonal to red...
if xp and yp and len(xp)==len(yp):
	for i in range(len(xp)):
		rgba[xp[i], yp[i], :3] = 1, 0, 0

heat_fig, heat_ax = plt.subplots()
im = heat_ax.imshow(rgba, interpolation='nearest')
heat_fig.savefig('heat.png')

#plt.show()
