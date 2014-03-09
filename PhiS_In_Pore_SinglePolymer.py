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
dfs = [0.0, 5.0, 10.0, 15.0]

# Define eqn that is to be solved numerically
#Variables: p=partition coefficent, phi=vol fract of short, phib=vol frac big
#df=Free energy cost by pore for entering pore, Ns and Nb are number of subunits
# for short and big polymers
def f(p, phi, df):
	""" Implicit definition of P(phi, dF) """
	#return -log(p) - df + (p-1)*phi + \
	#    ( phi*(1-p) + \
	#      5./4*alpha*(phi*p)**(5./4)*(p*phi-9./5) - \
	#      (phib/Nb)-5./4*alpha*(phi+phib)**(5./4)*(phi+phib-9./5) - \
	#      1./2*((1-p*phi)**2 - (1-phi-phib)**2) ) * Ns

#	return -p + exp(- df + (p - 1)*phi + Ns*(phi*(1 - p) + phib + (5./4)*alpha*p**(5./4)*phi**(5./4)*(p*phi - (9./5))- \
#		phib/Nb - (5./4)*alpha*(phi + phib)**(5./4)*(phi + phib-(9./5))-\
#		(1./2)*((1 - p*phi)**(2) - (1 - phi - phib)**(2)))) 
	return -p + exp(- df + Ns*(log((1 - p*phi)/(1 - phi)) + (p - 1)*phi + (9./4)*alpha*(1 - p**(5./4))*phi**(5./4)))
#Optimization algorithms, if opposite sign use brent, else use newton
def P(phi, df):
	""" Numerically solve for partition coefficient as a
	    function of \phi_s """
	if f(0,phi,df)*f(1,phi,df) < 0:
		return opt.brentq(f, 0, 1, args=(phi,df), maxiter=500) # Brent's method
		#return opt.brenth(f, 0, 1, args=(phi,df)) # Brent's method
		#return opt.bisect(f, 0, 1, args=(phi,df)) # Bisection method
		#x,r = opt.bisect(f, 0, 1, args=(phi,df), full_output=True) # Bisection method
		#print r.iterations
		#return x
	else:
		return opt.newton(f, 1.0, args=(phi,df), maxiter=500) # Newton-Raphson

# linestyles for different phi_big, colors for different deltaf's
lts = [':','-.','--','-']
colors = ['r','g','b','k','c','m','y']

pl.figure()
#for i,phi_b in enumerate(phi_bs):

#	pl.subplot(len(phi_bs),1,i+1)
#	for df in dfs:
for j,df in enumerate(dfs):
	phis = np.linspace(0, 0.99, 100)
	try: ps = [P(phi,df) for phi in phis]
	except: continue
	resids = [abs(f(p,phi,df)) for p,phi in zip(ps, phis)]
	max_resid = max(resids)
	print 'Largest residual for df=%f: %e' % (df, max_resid)
		#if j==3:
		#	labels=r'$\Phi_{PEG3500} = %.2f$' % phi_b
		#else:
		#	labels=None
	pps = [phi*P(phi,df) for phi in phis]
	pl.plot(phis,pps,color=colors[j])#,linestyle = lts[i],label=labels)#color = colors[j],
	pl.axis([0, 1.0, -0.1, 1.2])

#font = {'size': 22}
#matplotlib.rc('font',**font)
pl.text(0.02,0.98,r'$\Delta F = 0.0$', color='r')
pl.text(0.02,0.94,r'$\Delta F = 5.0$', color='g')
pl.text(0.02,0.89,r'$\Delta F = 10.0$', color='b')
pl.text(0.02,0.85,r'$\Delta F = 15.0$', color= 'k')

pl.legend(loc='upper left')
pl.ylabel('Number Fraction Of PEG400 In Pore   $\Phi_{PEG400}$(In)')
pl.xlabel('Number Fraction of PEG400 in Bulk   $\Phi_{PEG400}$(Out)')
pl.title('Number Fraction $\Phi_{PEG400}$(In) as a function of $\Phi_{PEG400}$(Out), $\Delta$F')

pl.savefig('PhiS_In_Pore_SinglePolymer_Generated.eps', dpi=600)


pl.show()



