#!/usr/bin/python

import numpy as np
from numpy import log,exp
import scipy.optimize as opt
import matplotlib.pyplot as pl
import matplotlib

Vbar = 1.0
alpha = 0.49
Ns = 9
Nb = 79
dfs = [0.0, 5.0,10.0, 15.0]
phi_bs = [0.05, 0.10, 0.15, 0.20]
phi_water_min = 0.01

def f(p, phi, phib, df):
	""" Implicit definition of P(phi, dF) """
	return - p + exp( - df + Ns*(log((1 - p*phi)/(1 - phi - phib)) + \
		(p - 1)*phi - phib + (9./4)*alpha*((phi + phib)**(5./4) - (p*phi)**(5./4))))


def Mu_s(phi,phib): 
	return phi + exp(1 - phi - (1-phi-phib)*Ns - phib*(Ns / Nb) +  \
	(1./2)*Ns*(1-phi-phib)**(2) - (5./4)*alpha*Ns*(phi + phib)**(9./4) + \
	(9./4)*alpha*Ns*(phi + phib)**(5./4))

#Optimization algorithms, if opposite sign use brent, else use newton
def P(phi, phib, df):
	""" Numerically solve for partition coefficient as a
	    function of \phi_s """
	if f(0,phi,phib,df)*f(1,phi,phib,df) < 0:
		return opt.bisect(f, 0, 1, args=(phi,phib,df), maxiter=500) # Bisection method
	else:
		return opt.newton(f, 1.0, args=(phi,phib,df), maxiter=5000) # Newton-Raphson

lts = [':','-.','--','-']
colors = ['k','r','b','g','c','m','y']

pl.figure()
for i,phi_b in enumerate(phi_bs):
	phis = np.linspace(0.0, 1.0-phi_b-phi_water_min, 100)
	cp = [Mu_s(phi,phi_b) for phi in phis]
	for j,df in enumerate(dfs):
		try: ps = [P(phi,phi_b,df) for phi in phis]
	       	except: continue
		resids = [abs(f(p,phi,phi_b,df)) for p,phi in zip(ps, phis)]
		max_resid = max(resids)
		print 'Largest residual for df=%f: %e' % (df, max_resid)
		if j==0:
			labels=r'$\Phi_{PEG3500} = %.2f$' % phi_b
		else:
			labels=None
		print 'Last phi tried for df=%f: phis=%.2f: phib=%.1f: total phi=%.2f' % (df, phi,phi_b,phi+phi_b)
		pl.plot(cp,ps,color=colors[j],linestyle = lts[i],label=labels)#color = colors[j],

pl.axis([0.0, 100, 0.0, 2.0])


pl.text(3.5,1.76,r'$\Delta f = 0.0$', color='k')
pl.text(3.5,1.64,r'$\Delta f = 5.0$', color='r')
pl.text(3.5,1.52,r'$\Delta f = 10.0$', color='b')
pl.text(3.5,1.40,r'$\Delta f = 15.0$', color= 'g')

pl.legend(loc='upper right')
pl.title('Partition Coefficient vs Chemical Potential, Fixed $\Phi_{PEG3500}$(out)')
pl.xlabel('Chemical Potential')
pl.ylabel('Partition Coefficient for PEG400')

pl.savefig('VaryShort_PC_vs_CP.eps', dpi=600)


pl.show()




