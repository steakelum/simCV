from covar_sim import CVModel
from random import seed
from csv import writer

#ys = 	[1,1,2,1,8,9,6,7,4,3,0,4,1,0,2,2,3]	# ds1 data set fcs in intervals
if True:
	xs = [	[0.0531,0.0619,0.1580,0.0810,1.0460,1.7500,2.9600,4.97,0.42,4.7,0.9,1.5,2.0,1.2,1.2,2.2,7.6], # eVec
			[4,20,1,1,32,32,24,24,24,30,0,8,8,12,20,32,24], # fVec
			[1,0,0.5,0.5,2.0,5,4.5,2.5,4,2,0,4,6,4,6,10,8]		] # cVec
else:
	xs = [	[0.05, 1, 0.19, 0.41, 0.32, 0.61, 0.32, 1.83, 3.01, 1.79, 3.17,3.4, 4.2, 1.2], # eVec
			[1.3, 17.8, 5.0, 1.5, 1.5, 3.0, 3.0, 8, 30, 9, 25, 15, 15, 2], # fVec
			[0.5, 2.8, 1, 0.5, 0.5, 1, 0.5, 2.5, 3.0, 3.0, 6, 4, 4, 1]	] # cVec


seed(1)



covar = CVModel(	model = "GM",
					params = [ 0.0261919972473019 ],
					covariates = xs,
					omega = 55.1228111682891,
					betas = [ 0.217742237478512, 0.0313272153762373, 0.177020429713746 ])


print(sum([covar.pr(i) for i in range(covar.length)]))
covar.gen_defects()
print(covar.LL())
print(covar.fcs)
print(f'Cummulative FC is: {sum(covar.fcs)}')
#print(covar.LL().subs(b, 0.5))



with open('poiss_sym.csv', 'w') as myfile:
	wr = writer(myfile)
	wr.writerow(['T', 'FC'] + [f'x{i+1}' for i in range(len(covar.covariates))])
	wr.writerows(zip(
		range(1,covar.length+1),
		covar.fcs,
		*covar.covariates))
myfile.close()