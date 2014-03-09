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
dfs = [0.0, 5.0,10.0, 15.0]
phi_bs = [0.0, 0.10, 0.20, 0.70]
phi_water_min = 0.01

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

	#return -p + exp(- df + (p - 1)*phi + Ns*(phi*(1 - p) + phib + (5./4)*alpha*p**(5./4)*phi**(5./4)*(p*phi - (9./5))- \
	#	phib/Nb - (5./4)*alpha*(phi + phib)**(5./4)*(phi + phib-(9./5))-\
	#	(1./2)*((1 - p*phi)**(2) - (1 - phi - phib)**(2)))) 
	##return -p + exp(- df + Ns*(log((1 - p*phi)/(1 - phi)) + (p - 1)*phi + (9./4)*alpha_tilda*(1 - p**(5./4))*phi**(5./4))
	return - p + exp( - df + Ns*(log((1 - p*phi)/(1 - phi - phib)) + \
		(p - 1)*phi - phib + (9./4)*alpha*((phi + phib)**(5./4) - (p*phi)**(5./4))))

#Optimization algorithms, if opposite sign use brent, else use newton
def P(phi, phib, df):
	""" Numerically solve for partition coefficient as a
	    function of \phi_s """
	if f(0,phi,phib,df)*f(1,phi,phib,df) < 0:
		#print 'bisect'
		#return opt.brentq(f, 0, 1, args=(phi,df)) # Brent's method
		#return opt.brenth(f, 0, 1, args=(phi,df)) # Brent's method
		return opt.bisect(f, 0, 1, args=(phi,phib,df), maxiter=500) # Bisection method
		#x,r = opt.bisect(f, 0, 1, args=(phi,df), full_output=True) # Bisection method
		#print r.iterations
		#return x
	else:
		#print 'newton'
		return opt.newton(f, 1.0, args=(phi,phib,df), maxiter=5000) # Newton-Raphson
		#return opt.bisect(f, 0, 1, args=(phi,phib,df)) # Bisection method
		# linestyles for different phi_big, colors for different deltaf's

lts = [':','-.','--','-']
colors = ['r','g','b','k','c','m','y']

pl.figure()
for i,phi_b in enumerate(phi_bs):
	phis = np.linspace(0.0, 1.0-phi_b-phi_water_min, 100)
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
		pps = [phi*P(phi,phi_b,df) for phi in phis]
		print 'Last phi tried for df=%f: phis=%.2f: phib=%.1f: total phi=%.2f' % (df, phi,phi_b,phi+phi_b)
		pl.plot(phis,pps,color=colors[j],linestyle = lts[i],label=labels)#color = colors[j],

		#pl.text(0.9,1.95,r'$\Delta F = %.2f$' % j(0), color='r')
		#pl.text(0.9,1.85,r'$\Delta F = %.2f$' % j(1), color='g')
		#pl.text(0.9,1.75,r'$\Delta F = %.2f$' % j(2), color='b')
		#pl.text(0.9,1.65,r'$\Delta F = %.2f$' % j(3), color= 'k')


	pl.axis([0, 1.05, -0.1, 1.05])

#font = {'size': 22}
#matplotlib.rc('font',**font)

pl.text(0.9,0.97,r'$\Delta F = 0.0$', color='r')
pl.text(0.9,0.93,r'$\Delta F = 5.0$', color='g')
pl.text(0.9,0.89,r'$\Delta F = 10.0$', color='b')
pl.text(0.9,0.85,r'$\Delta F = 15.0$', color= 'k')

pl.legend(loc='lower right')
pl.title('Number Fraction of PEG400 in Pore for Varied $\Phi_{PEG400}$(Out), Fixed $\Phi_{PEG3500}$, $\Delta$F')
pl.xlabel('Number fraction of PEG400 $\Phi_{PEG400}$(Out)')
pl.ylabel('Number Fraction of $\Phi_{PEG400}$(In)')

pl.savefig('PhiS_In_Pore_Vary_Short.eps', dpi=600)


pl.show()



