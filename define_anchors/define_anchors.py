#!/Users/Orfeo/anaconda3/bin/python3

import sys, argparse, re
import numpy as np
import pylab as pl
import matplotlib.pyplot as plt
from matplotlib import collections  as mc
from scipy.interpolate import spline
from scipy.signal import argrelextrema
#from windowaverage import weighted_moving_average

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

parser = argparse.ArgumentParser()
parser.add_argument('-a', '--alignment', nargs='?')
parser.add_argument('-p', '--profiles', nargs='?')
parser.add_argument('-c', '--columns', nargs='*')
parser.add_argument('-x', '--xcoord', nargs=2)
parser.add_argument('-o', '--output', nargs='?')
parser.add_argument('-fasta', dest='fasta', action='store_true')
parser.set_defaults(output='aligned_profiles.png')
parser.set_defaults(fasta=False)
parsed = parser.parse_args()

colorcodes = ['#08306B', '0.5']
heights = [1.15, 1.12]
#heights = [0.1, 0.05]

emptyline = re.compile('^\s*$')

f = open(parsed.alignment, 'r')
text = f.read().split('\n')
f.close()

if parsed.fasta:
	ch_idx = -1
	chain = {}
	for nl in range(len(text)):
		if emptyline.match(text[nl]):
			continue
		elif text[nl][0] == ">":
			ch_idx = ch_idx+1
			continue
		elif ch_idx in chain:
			chain[ch_idx].extend(text[nl])
		else:
			chain[ch_idx] = list(text[nl])
else:
	for nl in range(len(text)):
		if nl == 0:
			ch_idx = 0
			chain = {}
			continue
		if emptyline.match(text[nl]):
			ch_idx = 0
			continue
		if text[nl] != "":
			if "*" in text[nl] or text[nl].split()[0] == "Sequence" or text[nl].split()[0] == "Matched":
				continue
			if ch_idx in chain:
				chain[ch_idx].extend(text[nl].split()[1])
			else:
				chain[ch_idx] = list(text[nl].split()[1])
			ch_idx = ch_idx+1

#for i in chain.keys():
#	print(chain[i], len(chain[i]))

f = open(parsed.profiles, 'r')
text = f.read().split('\n')
f.close()

pr_lines = {}
templine = {}
tempax = {}
for c in parsed.columns:
	c = int(c)-1
	templine[c] = []
	pr_lines[c] = []
	tempax[c] = []

product = []
conversion1 = {}
conversion2 = {}
back_conversion = {}
e = {}
n_prod = 0
c1_count, c2_count = 0, 0
for nl in range(len(text)):
	if emptyline.match(text[nl]) or text[nl][0] == "#":
		continue
	else:
		for c in parsed.columns:
			c = int(c)-1
			if text[nl].split()[c] != "?0":
				templine[c].append(float(text[nl].split()[c]))
				tempax[c].append(int(text[nl].split()[0]))
				e[c] = 'yes'
			else:
				if templine[c]:
					pr_lines[c].append((tempax[c], templine[c]))
					templine[c] = []
					tempax[c] = []
				e[c] = ''
		if len(parsed.columns) == 2:
			c1 = int(parsed.columns[0]) - 1
			c2 = int(parsed.columns[1]) - 1
			back_conversion[n_prod] = (0,0)
			if text[nl].split()[c1] == "?0":
				a1 = -2.0
			else:
				conversion1[c1_count] = n_prod
				back_conversion[n_prod] = (c1_count, back_conversion[n_prod][1])
				c1_count += 1
				a1 = float(text[nl].split()[c1])
			if text[nl].split()[c2] == "?0":
				a2 = -2.0
			else:
				conversion2[c2_count] = n_prod
				back_conversion[n_prod] = (back_conversion[n_prod][0], c2_count)
				c2_count += 1
				a2 = float(text[nl].split()[c2])
			product.append( (a1 + 2.0) * (a2 + 2.0) - 4.0)
			n_prod += 1
for c in parsed.columns:
	c = int(c)-1
	if e[c]:
		pr_lines[c].append((tempax[c], templine[c]))


#if parsed.anchors:
#	f = open(str(parsed.anchors), 'r')
#	text = f.read().split('\n')
#	anchors = {}
#	for nl in range(len(text)):
#		fields = text[nl].split()
#		if nl == 0:
#			n_anchors = int(fields[0])
#		elif nl <= n_anchors:
#			anchors[(int(fields[0]), int(fields[1]))] = float(fields[2])

#print(parsed.columns, pr_lines)

minval = 1000000
maxval = -1000000
for c in parsed.columns:
	c = int(c)-1
	for x, y in pr_lines[c]:
		if minval > min(y):
			minval = min(y)
		if maxval < max(y):
			maxval = max(y)
	for i in range(len(pr_lines[c])):
		for j in range(len(pr_lines[c][i][1])):
			pr_lines[c][i][1][j] = (pr_lines[c][i][1][j] - minval)/(maxval - minval)


fig = plt.figure(1, (10,3))
ax = plt.subplot(111)
c_idx = 0
#ax2 = plt.axes([0,0,1,1], axisbg=(1,1,1,0))
for c in parsed.columns:
	c = int(c)-1
	for x, y in pr_lines[c]:
		if len(x) < 2:
			x = [x[0]-0.5] + x
			x.append(x[-1]+0.5)
			y = [y[0]] + y
			y.append(y[-1])
			ax.plot(x, y, lw=2.5, color=colorcodes[c_idx])
			ax.fill_between(x, 0, y, facecolor=colorcodes[c_idx], alpha=0.5)
		else:
			x = list(np.arange(x[0]-0.5, x[0]-0.05, 0.1)) + x + list(np.arange(x[-1]+0.05, x[-1]+0.5, 0.1))
			y = list(np.ones(5)*y[0]) + y + list(np.ones(5)*y[-1])
#			print(x,y)
			xnew = np.arange(x[0], x[-1]+0.05 ,0.05)
			smooth = spline(x, y, xnew)
#			print(xnew[-5:],smooth[-5:])
			ax.plot(xnew, smooth, lw=2.5, color=colorcodes[c_idx])
			ax.fill_between(xnew, 0, smooth, facecolor=colorcodes[c_idx], alpha=0.5)
		h = [ heights[c_idx] for i in range(len(x)) ]
#		ax2.plot(x, h, lw=2.5, color=colorcodes[c_idx])
		ax.plot(x, h, lw=2.5, color=colorcodes[c_idx])
	c_idx = c_idx+1

if (int(parsed.xcoord[1]) - int(parsed.xcoord[0]) <= 60):
	tag = []
	for i in range(len(chain[0])):
		tag.append(str(chain[0][i])+'\n'+str(chain[1][i]))
	ax.set_xticks(np.arange(len(chain[0])), minor=False)
	ax.set_xticklabels(tag, minor=False)

#ax.set_xticks(np.arange(len(chain[0])), minor=False)
#ax.set_xticklabels(tag, minor=False)
ax.set_xlim(int(parsed.xcoord[0]), int(parsed.xcoord[1]))
ax.set_ylim(0,1.2)

fig.subplots_adjust(bottom=0.15)

plt.savefig(parsed.output)

fig = plt.figure(2, (10,3))
ax = plt.subplot(111)

d_thr = 0.15
if len(parsed.columns) == 2:
	product = [ weighted_moving_average(range(len(product)),product,i,15,'triangular') for i in range(len(product)) ]
	np_product = np.asarray(product)
	locmax = argrelextrema(np_product, np.greater)[0]
	d_product = np.gradient(np_product)
	peaks = []
	limits = []
	list_locmax = locmax.tolist()
	for i in range(len(list_locmax)):
		if product[list_locmax[i]] > -2:
#			print(i, (back_conversion[list_locmax[i]][0], back_conversion[list_locmax[i]][1]))
			l_lim, r_lim = list_locmax[i], list_locmax[i]
			for pos in range(list_locmax[i], len(product)):
#				print(pos, list_locmax[i], d_product[pos])
				if abs(d_product[pos]) > d_thr:
					break
				else:
					r_lim = pos
			for pos in reversed(range(list_locmax[i])):
				if abs(d_product[pos]) > d_thr:
					break
				else:
					l_lim = pos
			limits.append(((back_conversion[l_lim][0], back_conversion[r_lim][0]), (back_conversion[l_lim][1], back_conversion[r_lim][1])))
			peaks.append((back_conversion[list_locmax[i]][0], back_conversion[list_locmax[i]][1], list_locmax[i]))
	ax.plot(range(len(product)), product)
#	for a in anchors:
#		plt.axvline(x=conversion1[a[0]], color='r')
#	for a in anchors:
#		plt.axvline(x=conversion2[a[1]], color='k')
	f = open('new_anchors_tmp.txt', 'w')
	f.write("{0}\n".format(len(peaks)))
	for np in range(len(peaks)):
		f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format(peaks[np][0], peaks[np][1], 0, limits[np][0][0], limits[np][0][1], limits[np][1][0], limits[np][1][1]))
		plt.axvline(x=peaks[np][2], color='g')


plt.savefig('product_and_anchors.png')
