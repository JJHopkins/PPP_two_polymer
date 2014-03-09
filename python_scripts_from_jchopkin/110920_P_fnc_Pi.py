#!/usr/bin/python

import numpy as np
from numpy import log,exp
import scipy.optimize as opt
import matplotlib.pyplot as pl
import matplotlib

Vbar = 1.0 # mL/g
alpha = 0.49
Nb = 35000. / 400 * 9 # Monomers / big molecule
Ns = 1000. / 400 * 9 # Monomers / small molecule
dfs = [0.0, 1.0, 2.0, 3.0]
phi_bs = [0.0, 0.10, 0.20, 0.30]

def f(p, phi, phib, df):
	""" Implicit definition of P(phi, dF) """
	return -p + exp(- df + (p - 1)*phi + Ns*(phi*(1 - p) + phib + (5./4)*alpha*p**(5./4)*phi**(5./4)*(p*phi - (9./5))- \
		phib/Nb - (5./4)*alpha*(phi + phib)**(5./4)*(phi + phib-(9./5))-\
		(1./2)*((1 - p*phi)**(2) - (1 - phi - phib)**(2)))) 

#Optimization algorithms, if opposite sign use brent, else use newton
def P(phi, phib, df):
	""" Numerically solve for partition coefficient as a
	    function of \phi_s """
	if f(0,phi,phib,df)*f(1,phi,phib,df) < 0:
		return opt.bisect(f, 0, 1, args=(phi,phib,df)) # Bisection method
	else:
		return opt.newton(f, 0.5, args=(phi,phib,df)) # Newton-Raphson

def osmotic_pressure(s, b):
	return s/Ns + b/Nb + (5./4)*alpha*(s + b)**(9./4)
	
# linestyles for different phi_big, colors for different deltaf's
lts = [':','-.','--','-']
colors = ['r','g','b','k','c','m','y']

pl.figure()
for i,phi_b in enumerate(phi_bs):
	phis = np.linspace(0, 1.0-phi_b, 100)
	ops= [osmotic_pressure(phi, phi_b) for phi in phis]
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
		#pl.semilogx(phis,ps,color=colors[j+1],linestyle = lts[i],label=labels)#color = colors[j],
                #pl.loglog(ops,ps,color=colors[j+2],linestyle = lts[i],label=labels)#color = colors[j],

		pps = [phi*P(phi,phi_b,df) for phi in phis]
		#pl.semilogx(ops,pps,color=colors[j],linestyle = lts[i],label=labels)#color = colors[j],
		pl.semilogx(ops,ps,color=colors[j],linestyle = lts[i],label=labels)#color = colors[j],
	#pl.axis([0, 1.0, -0.1, 1.0])
	pl.axis([0, 1.0, -0.1, 3.0])
#font = {'size': 22}
#matplotlib.rc('font',**font)
pl.text(0.02,0.68,r'$\Delta F = 0.0$', color='r')
pl.text(0.02,0.64,r'$\Delta F = 1.0$', color='g')
pl.text(0.02,0.59,r'$\Delta F = 2.0$', color='b')
pl.text(0.02,0.55,r'$\Delta F = 3.0$', color= 'k')

pl.legend(loc='upper left')
pl.title('Partition Coefficient vs Osmotic Pressure for Varied $\Phi_{PEG400}$(Out) and Fixed $\Phi_{PEG3500}$')
#pl.title('Number Fraction Of PEG400 In Pore for Varied $\Phi_{PEG400}$(Out) and Fixed $\Phi_{PEG3500}$')
pl.xlabel('Osmotic Pressure for Varied $\Phi_{PEG400}$(Out) and Fixed $\Phi_{PEG3500}$')
pl.ylabel('Partition Coefficient P as a function of $\Phi_{PEG400}$(Out), $\Phi_{PEG3500}$, $\Delta$F')
#pl.ylabel('Number Fraction $\Phi_{PEG400}$(In) as a function of $\Phi_{PEG400}$(Out), $\Phi_{PEG3500}$, $\Delta$F')

pl.savefig('110908_c_Part_Coeff_Mixed_VarySmall_Corrected.png')


pl.show()




