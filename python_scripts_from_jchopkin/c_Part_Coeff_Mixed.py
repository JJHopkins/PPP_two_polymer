#!/usr/bin/python

import numpy as np
from numpy import log,exp
import scipy.optimize as opt
import matplotlib.pyplot as pl
import matplotlib

# Things will go wrong a lot in the code below
# Ignore most of the warnings
np.seterr(all='ignore')

alpha = 0.49
Ns = 9
Nb = 79
dfs = [0.0, 2.0, 6.0, 10.0]
phi_bs = [0.0, 0.01, 0.05, 0.10]

# Define eqn that is to be solved numerically
#Variables: p=partition coefficent, phi=vol fract of short, phib=vol frac big
#df=Free energy cost by pore for entering pore, Ns and Nb are number of subunits
# for short and big polymers
def f(p, phi, phib, df):
	""" Implicit definition of P(phi, dF) """
	#return -log(p) - df + (p-1)*phi + \
	#    ( phi*(1-p) + \
	#      5./4*alpha*(phi*p)**(5./4)*(p*phi-9./5) - \
	#      (phib/Nb)-5./4*alpha*(phi+phib)**(5./4)*(phi+phib-9./5) - \
	#      1./2*((1-p*phi)**2 - (1-phi-phib)**2) ) * Ns
        return -log(p) - df + (p-1)*phi + \
	    ( phi*(1-p) + phib + \
	      (5./4)*alpha*phi**(9./4)*p**(9./4)- \
	      (9./4)*alpha*phi**(5./4)*p**(5./4) - \
	      (1./2)*(1 - p*phi)**2 - phib/Nb - (5./4)*alpha*(phi + phib)**(9./4) + \
	      (9./4)*alpha*(phi+phib)**(5./4)+ (1./2)*(1-phi-phib)**2) * Ns 


	#return -p + exp( - df + (p-1)*phi + \
	#    ( phi*(1-p) + phib + \
	#      (5./4)*alpha*phi**(9./4)*p**(9./4)- \
	#      (9./4)*alpha*phi**(5./4)*p**(5./4) - \
	#      (1./2)*(1 - p*phi)**2 - phib/Nb - (5./4)*alpha*(phi + phib)**(9./4) + \
	#      (9./4)*alpha*(phi+phib)**(5./4)+ (1./2)*(1-phi-phib)**2) * Ns) 


#Optimization algorithms, if opposite sign use brent, else use newton
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

# linestyles for different phi_big, colors for different deltaf's
lts = [':','-.','--','-']
colors = ['r','g','b','k','c','m','y']

pl.figure()
for i,phi_b in enumerate(phi_bs):
	phis = np.linspace(0, 1.0-phi_b, 100)
#	pl.subplot(len(phi_bs),1,i+1)
#	for df in dfs:
	for j,df in enumerate(dfs):
		try: ps = [P(phi,phi_b,df) for phi in phis]
	       	except: continue
		resids = [abs(f(p,phi,phi_b,df)) for p,phi in zip(ps, phis)]
		max_resid = max(resids)
		print 'Largest residual for df=%f: %e' % (df, max_resid)
		if j==3:
			labels=r'$\Phi_{PEG3500} = %.2f$' % phi_b
		else:
			labels=None
		pl.plot(phis,ps,color=colors[j],linestyle = lts[i],label=labels)#color = colors[j],
	pl.axis([0, 1.0, -0.1, 2.0])

#font = {'size': 22}
#matplotlib.rc('font',**font)
pl.text(0.02,1.07,r'$\Delta F = 0.0$', color='r')
pl.text(0.23,0.54,r'$\Delta F = 2.0$', color='g')
pl.text(0.38,.06,r'$\Delta F = 6.0$', color='b')
pl.text(.6,.02,r'$\Delta F = 10.0$', color= 'k')

pl.legend(loc='upper right')
pl.title('Partition Coefficent for Varied $\Phi_{PEG400}$ and Fixed $\Phi_{PEG3500}$')
pl.xlabel('Number fraction of PEG400 $\Phi_{PEG400}$')
pl.ylabel('Partition Coefficent  P( $\Phi_{PEG400}$, $\Phi_{PEG3500}$, $\Delta$F)')

pl.savefig('c_Part_Coeff_Mixed.png')


#pl.show()

