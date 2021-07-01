import random, sys
from sympy import exp, prod, log, var, symbols, N, factorial
from scipy.optimize import root

class CVModel():


	def poiss(self, lamb):
		# knuth's method for random poisson variate (ideal up to ~30)
		l = exp(-lamb)
		k = 1
		p = random.random()
		while p > l:
			k += 1
			p *= random.random()
		return k-1


	def h0(self, ivl):
		model = self.model
		params = self.params
		i = ivl + 1
		if model == "GM":
			b = params[0]
			return b 		# eqn 22

		elif model == "NB":
			b = params[0]
			return (i*b**2)/(1 + b*(i - 1))

		elif model == "DW2":
			b = params[0]
			return 1 - b**(i**2 - (i - 1)**2)

		elif model == "DW3":
			b, c = params
			return 1 - exp(-c * i**b)

		elif model == "S":
			p, pi = params
			return p * (1 - pi**i)

		elif model == "TL":
			c, d = params
			return (1 - exp(-1/d))/(1 + exp(- (i - c)/d))

		elif model == "IFRSB":
			c = params[0]
			return 1 - c / i

		elif model == "IFRGSB":
			c, alpha = params
			return 1 - c / ((i - 1) * alpha + 1)

		raise Exception("model not implemented")


	def g(self, i):	# eqn 15, verified
		return exp( sum([ self.betas[j]*self.covariates[j][i] for j in range(len(self.betas))] ) )# , evaluate=False)


	def pr(self, i):	# eqn 19, verified
		Sser = [(1 - self.h0(k))**self.g(k) for k in range(i)]	# note k goes from 0 to i-1
		return (1 - (1 - self.h0(i))**self.g(i)) * prod(Sser)


	def LL(self, eval = True):		# evaluate LL

		if self.fcs == None:
			raise Exception("non-existent failures")
		n = self.length
		return N(sum([	-self.omega * self.pr(i) 	+ \
						self.fcs[i] * log(self.omega) + \
						self.fcs[i] * log( self.pr(i) ) + \
						- log( factorial( self.fcs[i] ) )
				for i in range(n)]))


	def gen_defects(self, count = None, override = False):

		if count == None:
			count = self.length
		fcs = [ int(self.poiss( self.omega * self.pr(i) + 0.5) ) for i in range(count)] 

		if override or (self.fcs == None):
			self.fcs = fcs

		return fcs


	def __init__(self, model, params, covariates, omega, betas = None, fcs = None):
		self.model = model
		self.covariates = covariates
		self.params = params
		self.omega = omega
		self.betas = betas
		self.fcs = fcs

		self.length = len(covariates[0])	# store number of intervals

