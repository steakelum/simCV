from covar_sym import CVModel


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