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
dfs = [ -6.625, -6.4, -6.2, -5.0, -2.0, 0.0, 2.0, 4.0, 7.0, 10.0 ]
 
def f(p, phi, df):
	""" Implicit definition of P(phi, dF) """
	return -log(p) - df + (p-1)*phi + \
	    ( phi*(1-p) + \
	      5./4*alpha*(phi*p)**(5./4)*(p*phi-9./5) - \
	      5./4*alpha*phi**(5./4)*(phi-9./5) - \
	      1./2*((1-p*phi)**2 - (1-phi)**2) ) * Ns

def P(phi, df):
	""" Numerically solve for partition coefficient as a
	    function of \phi_s """
	if f(0,phi,df)*f(1,phi,df) < 0:
		#return opt.brentq(f, 0, 1, args=(phi,df)) # Brent's method
		#return opt.brenth(f, 0, 1, args=(phi,df)) # Brent's method
		return opt.bisect(f, 0, 1, args=(phi,df)) # Bisection method
		#x,r = opt.bisect(f, 0, 1, args=(phi,df), full_output=True) # Bisection method
		#print r.iterations
		#return x
	else:
		return opt.newton(f, 0.5, args=(phi,df)) # Newton-Raphson

phis = np.linspace(0, 1, 100)
pl.figure()
for df in dfs:
	try: ps = [P(phi,df) for phi in phis]
	except: continue

	resids = [abs(f(p,phi,df)) for p,phi in zip(ps, phis)]
	max_resid = max(resids)
	print 'Largest residual for df=%f: %e' % (df, max_resid)
# note: can drag to size fig then save fig() in terminal, saves last touched fig as is
#	pl.plot(phis, ps, label=r'$\Delta F = %g$' % df) # self selects number of digits
#	pl.plot(phis, ps, label=r'$\Delta F = %.3f$' % df) #give 3 floating digits
      	pl.semilogy(phis, ps, label=r'$\Delta F = %.3f$' % df) #give 3 floating digits

pl.axis([0, 1.0, -0.3, 6.0])
pl.legend(loc='lower right')
pl.title('Partition Coefficent as function of Number Fraction')
pl.xlabel(r'Number Fraction of Polymer $\phi_s$')
pl.ylabel(r'Partition Coefficent $P(\phi_s)$')
pl.savefig('Part_Coeff.png')


#pl.show()

