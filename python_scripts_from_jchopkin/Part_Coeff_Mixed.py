#!/usr/bin/python

import numpy as np
from numpy import log
import scipy.optimize as opt
import matplotlib.pyplot as pl

# Things will go wrong a lot in the code below
# Ignore most of the warnings
np.seterr(all='ignore')

alpha = 0.49
Ns = 9
Nb = 79
dfs = [ -6.625, -6.4, -6.2, -5.0, -2.0, 0.0, 2.0, 4.0, 7.0, 10.0 ]
phi_bs = [0.0, 0.25, 0.5, 0.75, 1.0]
 
def f(p, phi, phib, df):
	""" Implicit definition of P(phi, dF) """
	return -log(p) - df + (p-1)*phi + \
	    ( phi*(1-p) + \
	      5./4*alpha*(phi*p)**(5./4)*(p*phi-9./5) - \
	      (phib/Nb)-5./4*alpha*(phi+phib)**(5./4)*(phi+phib-9./5) - \
	      1./2*((1-p*phi)**2 - (1-phi-phib)**2) ) * Ns

def P(phi, phib, df):
	""" Numerically solve for partition coefficient as a
	    function of \phi_s """
	if f(0,phi,phib,df)*f(1,phi,phib,df) < 0:
		#return opt.brentq(f, 0, 1, args=(phi,df)) # Brent's method
		#return opt.brenth(f, 0, 1, args=(phi,df)) # Brent's method
		return opt.bisect(f, 0, 1, args=(phi,phib,df)) # Bisection method
		#x,r = opt.bisect(f, 0, 1, args=(phi,df), full_output=True) # Bisection method
		#print r.iterations
		#return x
	else:
		return opt.newton(f, 0.5, args=(phi,phib,df)) # Newton-Raphson

phis = np.linspace(0, .85, 85)

# lts = [':', '--', '-', '-.'], put in the i loop, and put linetype/style=lt into pl.plot(...,)
pl.figure()
for i,phi_b in enumerate(phi_bs):
	pl.subplot(1,len(phi_bs),i+1)
	for df in dfs:
		try: ps = [P(phi,phi_b,df) for phi in phis]
		except: continue
		resids = [abs(f(p,phi,phi_b,df)) for p,phi in zip(ps, phis)]
		max_resid = max(resids)
		print 'Largest residual for df=%f: %e' % (df, max_resid)
		pl.plot(phis,ps,label=r'$\Delta F = %.3f$' % df)
	pl.axis([0, 0.87, -0.3, 6.0])


pl.legend(loc='upper right')
pl.title('Partition Coefficent Peg 3500')
pl.xlabel('Number fraction of short polymer')
pl.ylabel('Partition Coefficent')

pl.savefig('Part_Coeff_Mixed.png')


pl.show()

