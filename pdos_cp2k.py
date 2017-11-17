""" class for reading, fitting, saving and plotting pdos files from CP2K .pdos files.
	written by Claudio Padilha @ University of York, Oct 2017
	version 1.0.0 """

import re, math
import numpy as np
import matplotlib.pyplot as plt
# import matplotlib.lines as mline

# from cp2k output, convert the energies from au to eV
au2ev = 2.72113838565563E+01


# gaussian function we want to fit, w is the weight (used for the projections)
def g(x, m, s, w):
	return w * np.exp(- 0.5 * ((x - m) / s) ** 2) / (np.sqrt(2 * math.pi) * s ** 2)


class Dos:
	"""A class that stores and smooths pdos data from Cp2k"""
	# the constructor

	def __init__(self, file=None):
		""" default constructor, it receives a filename in the format file.pdos
			if no file name is provided, an empty object is created"""
		if file is not None:
			# save the filename used to get the data
			self.file = file
			# read the data from the file and close it
			with open(file, "r") as data:
				# get a list with the lines of data
				lines = data.readlines()

			# number of dos used to create this object
			self.ndos = 1

			# find the atomic kind at the header (line 1)
			m = re.search('kind (.*) at', lines[0])
			if m:
				self.kind = str(m.group(1).strip())
			# find Ef at the header, convert it to eV
			n = re.search('E\(Fermi\) = (.*) a.u.', lines[0])
			if n:
				self.ef = float(n.group(1).strip()) * au2ev

			# depending on the atomic kind, you will have different orbitals. Extract them from line 2
			# and stored them in a dictionary. Set each orbital as a np array of zeros.
			self.orbs = dict((o, np.zeros(len(lines[2:]))) for o in re.sub('\s+', ' ', re.sub('\[a.u.\]', '', re.sub('#', ' ', lines[1].strip()).strip()).strip()).split(" ")[3:])

			# array to store the eigenvalues
			self.eig = np.zeros(len(lines[2:]))

			# read the file from line 3 on
			for i in range(2, len(lines)):
				# we have to turn many spaces ('\s+') into a single one and get a list from each line
				line = re.sub('\s+', ' ', lines[i].strip()).split(" ")
				# the eigenvalues are on the second column, convert them to eV and store
				self.eig[i - 2] = float(line[1]) * au2ev
				# the orbitals are on the forth column on
				for j in range(3, 3 + len(self.orbs)):
					self.orbs[list(self.orbs.keys())[j - 3]][i - 2] = float(line[j])
		else:
			self.file = None
			self.kind = None
			self.ef = None
			self.ndos = 0
			self.orbs = {}
			self.eig = []

		# the dos is not fitted when the constructor is called, you should call Dos.fit()
		self.fitted = False

	def copy(self):
		""" make a copy of the Dos object """
		cdos = Dos()
		cdos.file = self.file
		cdos.kind = self.kind
		cdos.ef = self.ef
		cdos.ndos = self.ndos
		cdos.eig = self.eig.copy()
		cdos.orbs = self.orbs.copy()
		if self.fitted:
			cdos.fitted = self.fitted
			cdos.sigma = self.sigma
			cdos.e = self.e
			cdos.dos = self.dos.copy()
		return cdos

	def fit(self, sigma=0.05, ngrid=100):
		""" this function fits gaussians of half-width sigma and mean values equal to the
		eigenvalues stored in the self.eig vector. ngrid is the multiplication factor used
		to derive the smooth energy grid from the self.eig grid"""

		# you can only fit a Dos object that contains any information
		if self.ndos > 0:
			# weight for the fitting
			w = 1 / self.ndos
			# save sigma into the object
			self.sigma = sigma
			# grid for the energies
			self.e = np.linspace(self.eig[0] - 4 * sigma, self.eig[-1] + 4 * sigma, len(self.eig) * ngrid)

			# the fitted projections are stores in a dictionary called dos
			# first we have all the ml contributions (s, px, py, pz, d-2, d-1, d0, etc.)
			self.dos = dict((o, np.zeros_like(self.e)) for o in self.orbs)

			# then we want the total projection per angular momenta
			#  s = s
			#  p = px + py + pz,
			#  d = d-2 + d-1 + d0 + d+1 + d+2
			# and so on...
			t = []
			for i in self.dos:
				if i[0:1] not in t:
					t.append(i[0:1])
			for i in t:
				self.dos[i] = np.zeros_like(self.e)

			# the total dos is stored here
			self.dos['total'] = np.zeros_like(self.e)

			# fit the dos using the gaussian function on the whole vectors
			for i in range(len(self.eig)):
				# the total dos is fitted here (weight = 1 / ndos)
				self.dos['total'] += g(self.e, self.eig[i], sigma, w)
				# the projected dos per orbital (ml) is fitted here (weights = orbital contribution / ndos)
				for orb in self.orbs:
					# print(orb, len(self.dos[orb]), i)
					self.dos[orb] += g(self.e, self.eig[i], sigma, w * self.orbs[orb][i])

			# sum the angular momenta (p = px + py + pz, d = d-2 + d-1 + d0 + d+1 + d+2, and so on)
			t = [o for o in self.dos if len(o) > 1 and o != 'total']
			for o in t:
				self.dos[o[0:1]] += self.dos[o]

			# if everything goes well, then this dos is fitted
			self.fitted = True
		else:
			print("You can only fit a Dos object which is not empty!")

		# returns true if the fitting process was ok, or false if not
		return self.fitted

	def save(self, file=None):
		""" saves the data after fitting into a text file """
		# check if the Dos object was fitted before saving it
		if not self.fitted:
			go = self.fit()
		else:
			go = True

		# if the fitting was ok
		if go:
			# build header
			header = '# fitted pdos from file ' + self.file + '\n# parameters: kind = ' + self.kind + ' Fermi-level = ' + str(self.ef) + ' sigma = ' + str(self.sigma) + '\n# energy[eV]'
			for o in self.dos:
				for i in range(12 - len(o)):
					header += " "
				header += o
			# user provided filename or use standard
			if file is None:
				file = self.file.split("/")[-1][:-4] + 'dat'
			# open the file where the data will be save into
			with open(file, "w") as data:
				# write header
				data.write(header)
				# for each energy value, write all the data
				for i in range(len(self.e)):
					data.write("\n%12f" % self.e[i])
					for o in self.dos:
						data.write('%12f' % self.dos[o][i])
		else:
			print("Problem with data fitting prior to saving. Check if Dos object is not empty.")


# plot Dos objects in subplots
def plot(dos, proj=None, xmin=-10, xmax=10, ymin=0, ymax=1000, show=False, hide_yticks=True, file=None, align=False):
	""" plots one or more pdos files into subplots
		dos should be a list of Dos objects
		proj should be a list of lists containing the projections in the same order
		as in the dos list. For example, if you call this function using

		>>> plot([dos1, dos2])

		in order to get the projections you should do something like this

		>>> plot([dos1, dos2], proj=[['s', 'p'], ['p', 'd']])

		and you'll get s and p orbitals for dos1 and p and d for dos2.

		if align == True you will align all graphs by the first eigenvalue of the first Dos """

	# The number of Dos objects in the list dos
	ndos = len(dos)

	# loop over all Dos objects
	for i in range(ndos):
		# if they are not fitted, fit using default parameters
		if not dos[i].fitted:
			go = dos[i].fit()
		else:
			go = True

		# if user wants to align all pdos by the core level
		if align:
			shift = dos[i].eig[0] - dos[0].eig[0]
		else:
			shift = dos[i].ef

		# if the fitting is ok
		if go:
			# add subplots in sequence
			plt.subplot(ndos, 1, i + 1)
			# plot total dos using a black thicker line
			plt.plot(dos[i].e - shift, dos[i].dos['total'], color='black', linewidth=2)
			# if you want projections...
			if proj is not None:
				if proj[i] is not None:
					# each element of the list proj shoud be a list with the projections
					for orb in proj[i]:
						# look for the projection in each Dos object
						if orb.lower() in dos[i].dos:
							# add plot for each projection found
							plt.plot(dos[i].e - shift, dos[i].dos[orb], label=dos[i].kind + '$_{' + orb + '}$')
						else:
							# if not found, tells the user and goes on
							print(' warning: orbital ' + orb + ' not found!')
				# add legend with the kind and projection info
				plt.legend()
			else:
				# if you did not ask for projections, just add the kind
				if dos[i].kind is not None:
					plt.legend(dos[i].kind)
			# add a vertical line at E_F
			plt.plot([dos[i].ef - shift, dos[i].ef - shift], [ymin, ymax], linestyle='-', color='blue')
			# show ylabel and title
			plt.ylabel('pDOS')

			# hide the x ticks and corresponding tick labels for all plots
			# except the last one at the bottom
			if i < ndos - 1:
				plt.tick_params(
					axis='x',         	# changes apply to the x-axis
					which='both',     	# both major and minor ticks are affected
					bottom='off',       # ticks along the bottom edge are off
					top='off',      	# ticks along the top edge are off
					labelbottom='off')  # labels along the bottom edge are off
			# if user requested, hide yticks
			if hide_yticks:
				plt.tick_params(
					axis='y',         # changes apply to the y-axis
					which='both',     # both major and minor ticks are affected
					left='off',       # ticks along the left edge are off
					right='off',      # ticks along the right edge are off
					labelleft='off')  # labels along the left edge are off
			# set boundaries
			plt.xlim(xmin, xmax)
			plt.ylim(ymin, ymax)

		else:
			print("Problem with data fitting prior to plotting. Check if all Dos objects are not empty.")
			return False

	# add x label to last plot at bottom
	plt.xlabel('$E-E_F$(eV)')

	# make it look good
	plt.tight_layout()

	# if user wants just to show it
	if show:
		plt.show()
	# otherwise save to file
	else:
		# if user provided a filename, use it
		if file is not None:
			plt.savefig(file, dpi=300)
		# otherwise use default name 'pdos_cp2k.png'
		else:
			plt.savefig('pdos_cp2k', dpi=300)

	# close it all
	plt.clf()
	plt.cla()
	plt.close()


# combines Dos objects into a single Dos (total dos only)
def combine(dos, align=0):
	""" combines two or more Dos objects and returns another Dos object

		** Remark: it outputs total dos only

		align stores the indices of the core levels used to
		align the energies (first eigenvalue if none is provided)."""

	if isinstance(align, list):
		if len(dos) == len(align):
			shift = [dos[i].eig[align[i]] - dos[0].eig[align[0]] for i in range(len(dos))]
		else:
			print("dos and align should have same length!")
			return None
	else:
		shift = [dos[i].eig[align] - dos[0].eig[align] for i in range(len(dos))]

	# copy dos[0] eigenvalues into a new list (we want to keep the original)
	tdos = [i for i in dos[0].eig]

	# we will shift and combine all Fermi energies in order to get the average
	ef = [dos[0].ef]

	# stack all aligned eigenvalues and all aligned EF's
	for i in range(1, 1 + len(dos[1:])):
		temp = dos[i].eig - shift[i]
		for j in temp:
			tdos.append(j)
		ef.append(dos[i].ef - shift[i])

	# we will return a Dos object, empty at start
	cdos = Dos()
	# number of Dos objects used to create this new Dos object
	cdos.ndos = len(dos)
	# the Ef will be the average
	cdos.ef = np.mean(ef)
	# the eigenvalues are the sorted stacked values from before
	cdos.eig = np.sort(np.array(tdos, dtype=float))
	# there is no kind
	cdos.kind = None
	# return the Dos
	return cdos
