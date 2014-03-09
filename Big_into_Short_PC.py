#!/usr/bin/python
#testing changes for git

import numpy as np
from numpy import log,exp
import scipy.optimize as opt
import matplotlib.pyplot as pl
import matplotlib

# Things will go wrong a lot in the code below
# Ignore most of the warnings
#np.seterr(all='ignore')

alpha = 0.49
Ns = 9
Nb = 79
dfs = [1.0, 2.0, 3.0, 4.0, 5.0]
phi_ss = np.linspace(0.05, 0.50, 50)
phibb = np.linspace(0.50, 0.05, 50)
phi_water_min = 0.01

def f(p, phi, phib, df):
	""" Implicit definition of P(phi, dF) """
	return - p + exp( - df + Ns*(log((1 - p*phi)/(1 - phi - phib)) + \
		(p - 1)*phi - phib + (9./4)*alpha*((phi + phib)**(5./4) - (p*phi)**(5./4))))

def P(phi, phib, df):
	""" Numerically solve for partition coefficient as a
	    function of \phi_s """
	if f(0,phi,phib,df)*f(1,phi,phib,df) < 0:
		return opt.bisect(f, 0, 1, args=(phi,phib,df), maxiter=500) # Bisection method
	else:
		return opt.newton(f, 1.0, args=(phi,phib,df), maxiter=5000) # Newton-Raphson

lts = [':','-.','--','-',':','-.','--','-',':','-.','--','-']
colors = ['k','r','b','g','c','m','y']

pl.figure()
#phibb = np.linspace(0.95, 0.05, 9)
#for i,phi_s in enumerate(phi_ss):
#	for k,phib in enumerate(phibb):
for j,df in enumerate(dfs):
	try: 
		ps = [P(phi_s,phib,df) for phi_s,phib in zip(phi_ss,phibb)]
	except Exception as err: 
		print 'tried and failed for this df=%.3f:ps=this combo bad:phi_s=%.3f:phib=%.3f' % (df,phi_s,phib)
		print err
		continue
	resids = [abs(f(p,phi_s,phib,df)) for p,phib in zip(ps, phibb)]	
	max_resid = max(resids)
	print 'Largest residual for df=%f: %e' % (df, max_resid)
#	if j==0:
#		labels=r'$\Phi_{PEG400} = %.2f$' % phi_s
#	else:
#		labels=None
	#pps = [phi_s*P(phi_s,phib,df) for phib in phibb]
	print 'Last phi tried for df=%f: phis=%.2f: phib=%.1f: total phi=%.2f' % (df, phi_s,phib,phi_s+phib)
	pl.plot(phibb,ps)#,color=colors[j],linestyle = lts[i],label=labels)#color = colors[j],
#		pl.axis([0.0, 1.0, -0.2, 10])

pl.text(0.02,5.50,r'$\Delta f = 0.0$', color='k')
pl.text(0.02,5.10,r'$\Delta f = 5.0$', color='r')
pl.text(0.02,4.70,r'$\Delta f = 10.0$', color='b')
pl.text(0.02,4.30,r'$\Delta F = 15.0$', color= 'g')

pl.legend(loc='upper left')
pl.title('Partition Coefficient for Fixed $\Phi_{PEG400}$(Out) and Varied $\Phi_{PEG3500}$')
pl.xlabel('Number Fraction of PEG3500')
pl.ylabel('Partition Coefficient for PEG400')

pl.savefig('110930_Part_Coeff_VaryBig_Complete_Generated.eps', dpi=600)


pl.show()




