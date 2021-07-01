import random, sys, matplotlib.pyplot as plt
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


	def sim_pr_fn(self, i, b = 5, tot = False):	# basic exponential shape from 0 to 1
		return (exp(b*i) - exp(-b)) / (1 - exp(-b))

	def sim_pr_solve(self, x, count):
		samples = [self.sim_pr_fn(i/count - 1)*x for i in range(count)]
		return float(sum(samples) - 1)

	def sim_pr(self, count):
		val = root(self.sim_pr_solve, 1, method='hybr', tol=1e-10, args=(count))
		samples = [self.sim_pr_fn(i/count - 1)*val.x for i in range(count)]
		plt.scatter(range(1, count+1),samples)



	def __init__(self, model, params, covariates, omega, betas = None, fcs = None):
		self.model = model
		self.covariates = covariates
		self.params = params
		self.omega = omega
		self.betas = betas
		self.fcs = fcs

		self.length = len(covariates[0])	# store number of intervals




#ys = 	[1,1,2,1,8,9,6,7,4,3,0,4,1,0,2,2,3]	# ds1 data set fcs in intervals
if True:
	xs = [	[0.0531,0.0619,0.1580,0.0810,1.0460,1.7500,2.9600,4.97,0.42,4.7,0.9,1.5,2.0,1.2,1.2,2.2,7.6],
			[4,20,1,1,32,32,24,24,24,30,0,8,8,12,20,32,24],
			[1,0,0.5,0.5,2.0,5,4.5,2.5,4,2,0,4,6,4,6,10,8]	]	# covariates DS1
else:
	xs = [	[0.05, 1, 0.19, 0.41, 0.32, 0.61, 0.32, 1.83, 3.01, 1.79, 3.17,3.4, 4.2, 1.2],
			[1.3, 17.8, 5.0, 1.5, 1.5, 3.0, 3.0, 8, 30, 9, 25, 15, 15, 2],
			[0.5, 2.8, 1, 0.5, 0.5, 1, 0.5, 2.5, 3.0, 3.0, 6, 4, 4, 1]	]

b, w, beta = var('b'), var('w'), symbols('beta:3')

random.seed(1)



covar = CVModel(	model = "NB",
					params = [0.0262],
					covariates = xs,
					omega = 55.1228,
					betas = [0.2177, 0.0313, 0.1770])
covar.sim_pr(17)

'''
covar = CVModel(	model = "GM",
					params = [0.01],
					covariates = xs,
					omega = 55.1228,
					betas = [0.2177, 0.0313, 0.1770])'''

plt.scatter(range(1,covar.length+1), [covar.pr(i) for i in range(covar.length)])
plt.show()

#covar.gen_defects()
#print(covar.LL())
#print(covar.fcs)
#print(covar.LL().subs(b, 0.5))