#! /usr/bin/env python 
import numpy, csv, random, math, scipy, scipy.interpolate
import matplotlib.pyplot as plt

class Flux(object):
	"""docstring for Flux"""
	def __init__(self):
		self._tempangles = numpy.array([1,0.6,0.4,0.3,0.2,0.1,0.05,0.0])
		self._angles = numpy.arccos(self._tempangles)
		self._datafile = list(csv.reader(open("muonflux2.csv", "rb" ), delimiter = '\t'))	# reads each row from the csv into a list, stores each list as an element of a list
		self._fluxtable = [item for sublist in self._datafile for item in sublist]	# flatten the list
		self._concatflux = []
		for i in xrange(len(self._fluxtable)/2):
			x = float(self._fluxtable.pop(0))
			y = float(self._fluxtable.pop(0))
			self._concatflux.append(x*10**y)
		
		self._fluxdata = numpy.array(self._concatflux).reshape(-1,9)	# takes flattened array and reshapes to a 9 column matrix with arbitrary row size
		self._cumsums = self._fluxdata
		# this loop iterates over sliced rows to replace them with cumulative sums of the slice. The slice is the whole row minus the first element
		for i in xrange(numpy.alen(self._cumsums)):
			self._cumsums[i,1:] = numpy.cumsum(self._cumsums[i,1:])
		
		
		self.interpgrid=scipy.interpolate.RectBivariateSpline(self._fluxdata[:,0],self._angles,self._fluxdata[:,1:],kx=1,ky=1)
		self.precision = 30
		self.newx = scipy.logspace(math.log10(self._fluxdata[0,0]),math.log10(self._fluxdata[-1,0]),self.precision)
		#print self.newx # the interpolated energies
		self.interpedgrid = self.interpgrid(self.newx,self._angles)
		
		self.interpcumsumgrid = scipy.interpolate.RectBivariateSpline(self._fluxdata[:,0],self._angles,self._cumsums[:,1:],kx=1,ky=1)
		self.interpedcumsums = self.interpcumsumgrid(self.newx,self._angles)
		
		
		self._ncs = numpy.cumsum(self.interpedcumsums[:,-1]) # takes the cumulative sum of the last column of cumsums[]
		
		self._ncs = self._ncs/self._ncs[-1] # normalizes to the last element which is the largest
		
		
		# this loop does the same thing, but for the rows
		for i in xrange(numpy.alen(self.interpedcumsums)):
			self.interpedcumsums[i,:] /= self.interpedcumsums[i,-1]
	
	def randFlux(self):
		"""docstring for randFlux"""
		
		self._randnum1 = random.random()
		self._randnum2 = random.random()
		self._energyindex = numpy.nonzero(self._ncs > self._randnum1)[0][0] # finds the index from the 1D array, ncs, which will be used to look up an energy value
		self._angleindex = numpy.nonzero(self.interpedcumsums[self._energyindex,:] > self._randnum2)[0][0] # same, but to look up an angle value
		self.phi = 2*numpy.pi*random.random()
		self.x = random.randint(-250000,250000)
		self.y = random.randint(-250000,250000)
		self.theta = self._angles[self._angleindex]
		self.energy = self.newx[self._energyindex] # self.newx is actually the interpolated energies
		return (self.energy, self.theta, self.phi, self.x,self.y) # returns energy, theta, phi, x, y as a tuple
	
	def tableSum(self):
		"""docstring for tableSum"""
		self.fsum = 0
		for i in xrange(numpy.alen(self.newx)):
			for j in xrange(numpy.alen(self._angles)):
				self.fsum += self.newx[i]*1.15*self._angles[j]*math.fabs(self.interpedgrid[i,j])
		
		
		return self.fsum

class Muon(object):
	"""docstring for Muon"""
	def __init__(self, energy, theta, phi, x, y):
		self.z = 250000
		self.dt = 1
		self.energy = energy
		self.theta = theta
		self.phi = phi
		self.x = x
		self.y = y
		self.__makeCartesian__()
		# self.t = numpy.arange(0.0, math.floor(self.r)+self.dt, self.dt)		
		# 		originarametric = numpy.empty((numpy.alen(self.t), 3))
		# 		self.__makeParametric__()
		# 	
	def __makeCartesian__(self):
		"""docstring for __makeCartesian__"""
		self.origin = numpy.array([self.x,self.y, self.z])
		self.r = math.sqrt(math.fabs(self.x)**2 + math.fabs(self.y)**2 + self.z**2)
		self.dx = self.r*math.sin(self.theta)*math.cos(self.phi)
		self.dy = self.r*math.sin(self.theta)*math.sin(self.phi)
		self.dz = (-1)*self.r*math.cos(self.theta)
		self.d = numpy.array([self.dx,self.dy,self.dz]) # direction vector in Cartesian coordinates
		self.dirvec = self.d/numpy.linalg.norm(self.d) # divide by its norm to turn into a unit vector
	
	def __makeSpherical__(self):
		"""docstring for __makeSpherical__"""
		self.ds = numpy.array([self.r,self.theta,originhi]) # direction vector in spherical coordinates
	
	# def __makeParametric__(self):
	# 		"""docstring for __makeParametric__"""
	# 		for i in xrange(numpy.alen(self.t)):
	# 			originarametric[i,:] = origin + self.t[i]*dirvec

class Detector(object):
	"""docstring for Detector"""
	def __init__(self):
		self.bmin = numpy.array([-50000, -50000,0])
		self.bmax = numpy.array([50000, 50000, 100000])
		self.tmin = 0
		self.tmax = 250000
		self.Ax = numpy.array([-50000, -50000, 0])
		self.Bx = numpy.array([50000, -50000, 0])
		self.Cx = numpy.array([-50000, -50000, 100000])
		self.Ay = numpy.array([50000, 50000,0])
		self.By = self.Bx
		self.Cy = numpy.array([50000, 50000, 100000])
		self.__getNormal__()
		
	def getIntersection(self, dirvec, origin):
		"""docstring for getIntersection"""
		if (numpy.inner(self.normalxz,dirvec) == 0 or numpy.inner(self.normalyz,dirvec) == 0):
			if(origin[0] > -50000 and origin[0] < 50000 and origin[1] > -50000 and origin[1] < 50000): 
				return 1
			else:
				return 0
		else:
			for i in xrange(numpy.alen(origin)):
				self.t1 = (self.bmin[i]-origin[i])/dirvec[i]
				self.t2 = (self.bmax[i]-origin[i])/dirvec[i]
				if (self.t1 > self.t2): self.t1,self.t2 = self.t2,self.t1
				if (self.t1 > self.tmin): self.tmin = self.t1
				if (self.t2 > self.tmax): self.tmax = self.t2
				if (self.tmin > self.tmax): return 0
			return 1
			
	
	def __getNormal__(self):
		self.normalxz = numpy.cross((self.Bx-self.Ax),(self.Cx-self.Ax))
		self.normalyz = numpy.cross((self.By-self.Ay),(self.Cy-self.Ay))


def getEE(gamma, mumass):
	"""docstring for getEE"""
	xstep = 0.001
	xpdf = numpy.arange(0.0, 0.053+xstep,xstep)
	G = 1.166*10**(-5)
	eepdf = ((G**2)/(12*math.pi**3))*(mumass**2)*(xpdf**2)*(3- ( (4*xpdf) / (mumass) ))
	eecdf = numpy.cumsum(eepdf)
	eecdf /= eecdf[-1]
	randnum = random.random()
	cdfindex = numpy.nonzero(eecdf > randnum)[0][0]
	electronenergy = xpdf[cdfindex]
	paramomentum = random.uniform(-1,1)*electronenergy
	#electronmomentum = 
	return electronenergy,paramomentum

def stepper(muon, detector):
	"""docstring for stepper"""
	decayed = False
	dt = 0.2 # in microseconds
	c = 29980 # speed of light in cm/us
	mumass = 105.66/1000 # muon mass in GeV
	energy = muon.energy
	gamma = energy/mumass
	lamb = 1/(gamma*2.2) # decay constant, inverse of the average lifetime of the muon
	beta = math.sqrt(gamma**2 -1)/gamma
	origin = muon.origin
	dirvec = muon.dirvec
	currpoint = origin
	#print 'Point origin:', currpoint
	#print 'Direction unit vector:', dirvec
	while (decayed == False):
		r = random.random()
		if r > lamb*dt and energy > 1.12:
			dist = beta*dt*c
			deltenergy = (0.002*0.93 + 0.000007*energy*0.93)*dist
			#print 'Energy is', energy, 'dist is', dist, 'Delta energy is ', deltenergy
			energy -= deltenergy
			#print 'energy is ', energy
			gamma = energy/mumass
			if (gamma < 1) : gamma = 1
			#print 'gamma is ', gamma
			beta = math.sqrt(gamma**2 -1)/gamma
			#print 'beta is ', beta
			currpoint += dist*dirvec
			#print currpoint
		else:
			electronenergy, paramomentum = getEE(gamma,mumass)
			energyelectronlab = gamma*electronenergy + beta*gamma*paramomentum
			#print 'Muon energy at decay time is ', energy, 
			#print 'Decay coordinates are ', currpoint[2]
			#print 'Beta is ', beta
			#print 'gamma is ', gamma
			decayed = True
	return energyelectronlab, currpoint[2]

def gettime(tablesum,numgen):
	"""docstring for gettime"""
	time = numgen/((2.5*10**(11))*tablesum )
	return time

def main():
	"""docstring for main"""
	histogram = []
	histoallenergy = []
	height = []
	height100 = []
	num = 1*10**6
	index = 0
	flux = Flux()
	detector = Detector()
	for i in xrange(num):
		energy, theta, phi, x, y = flux.randFlux()
		muon = Muon(energy,theta,phi,x,y)
		if (detector.getIntersection(muon.dirvec,muon.origin)):
			energyelectronlab, depth = stepper(muon, detector)
			height.append(depth)
			histoallenergy.append(energyelectronlab)
			#print 'Electron energy in the lab frame is ', energyelectronlab
			if (energyelectronlab >= 100):
				#print depth, energyelectronlab, energy, theta
				histogram.append(energyelectronlab)
				height100.append(depth)
				#print energyelectronlab
			index += 1
		else:
			pass
	print 'Total number out of events that should intersect the detector: %d out of %d' % (index,num)
	tablesum = flux.tableSum()
	gentime = gettime(tablesum,num)
	height = numpy.array(height)
	histo = numpy.array(histogram)
	historate = histo/gentime
	histoallenergy = numpy.array(histoallenergy)
	height100 = numpy.array(height100)
	print gentime
	#histoallenergy /= gentime
	gentime=numpy.array([gentime])
	
	
	numpy.savetxt("histoallenergy.txt",histoallenergy)
	numpy.savetxt("height.txt",height)
	numpy.savetxt("histo100gev.txt", histo)
	numpy.savetxt("height100gev.txt", height100)
	numpy.savetxt("gentime.txt", gentime)
	 

	

if __name__ == "__main__":
	main()
