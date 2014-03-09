#!/usr/bin/python

import numpy as np
from numpy import log
from scipy.optimize import newton, brentq
import matplotlib.pyplot as pl

alpha = 0.49
Ns = 9
df = +0.3

def f(p, phi):
	""" Implicit definition of P(phi, dF) """
	return -log(p) - df + (p-1)*phi + \
	       ( phi*(1-p) + \
 		 5./4*alpha*(phi*p)**(5./4)*(p*phi-9./5) - \
		 5./4*alpha*phi**(5./4)*(phi-9./5) - \
		 1./2*((1-p*phi)**2 - (1-phi)**2) ) * Ns

def P(phi):
	""" Numerically solve for partition coefficient as a
	    function of \phi_s """
	return newton(f, 0.5, args=(phi,)) # Newton-Raphson
	#return brentq(f, 0, 1, args=(phi,)) # Brent's method

phis = np.linspace(0, 1, 100)
ps = [P(phi) for phi in phis]

print 'Largest residual', max(abs(f(p,phi)) for p,phi in zip(ps,phis))

pl.plot(phis, ps)
pl.xlabel(r'$\phi_s$')
pl.ylabel(r'$P(\phi_s)$')
pl.show()

