#!/usr/bin/python

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
dfs = np.linspace(0.0,10,60)# 100)
phi_bs = [0.05, 0.10, 0.15, 0.20]
phi_water_min = 0.01
phi = 0.20
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

lts = [':','-.','--','-']
#markers = ['x' , '^' , 'd' ]
#colors = ['k','r','b','g','c','m','y']

pl.figure()
for i,phi_b in enumerate(phi_bs):
	for j,df in enumerate(dfs):
		try: ps = [P(phi,phi_b,df) for df in dfs]
       		except: continue
		resids = [abs(f(p,phi,phi_b,df)) for p,df in zip(ps, dfs)]
		max_resid = max(resids)
#		print 'Largest residual for df=%f: %e' % (df, max_resid)
		if j==0:
			labels=r'$\phi_{b} = %.2f$' % phi_b
		else:
			labels=None
		#pps = [phi*P(phi,phi_b,df) for phi in phis]
		#print 'Last phi tried for df=%f: phis=%.2f: phib=%.1f: total phi=%.2f' % (df, phi,phi_b,phi+phi_b)
		#pl.plot(phis,ps,color=colors[j],linestyle = lts[i],label=labels)#color = colors[j],
	pl.plot(dfs,ps,'k',linestyle = lts[i], linewidth= 0.9,label=labels)#, color = colors[i])
#		pl.text(0.9,1.95,r'$\phi_{PEG3500} = %.2f$' % j(0), color='k')
#		pl.text(0.9,1.95,r'$\phi_{PEG3500} = %.2f$' % j(1), color='r')
#		pl.text(0.9,1.95,r'$\phi_{PEG3500} = %.2f$' % j(2), color='b')


	pl.axis([0, 10.0, 0.0, 1.8])

#pl.text(0.85,3.75,r'$\Delta f = 0.0$', color='k')
#pl.text(0.85,3.45,r'$\Delta f = 5.0$', color='r')
#pl.text(0.85,3.15,r'$\Delta f = 10.0$', color='b')
#pl.text(0.85,2.85,r'$\Delta f = 15.0$', color= 'g')

pl.legend(loc='upper right')
#pl.title('Partition Coefficient vs. Free Energy of Pore')
pl.xlabel(r'$\Delta f$', size = 'x-large')
pl.ylabel(r'$p  ( \phi_{s}  , \phi_{b} , \Delta f )$', size = 'x-large')

pl.savefig('VaryDeltaF_BW_PC.eps', dpi=600)


pl.show()



