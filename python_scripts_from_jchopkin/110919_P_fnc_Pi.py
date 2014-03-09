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
Vbar = 1.0 # mL/g
Nb = 35000. / 400 * 9 # Monomers / big molecule
Ns = 1000. / 400 * 9 # Monomers / small molecule
dfs = [0.0, 1.0, 2.0, 3.0]
phi_bs = [0.0, 0.10, 0.20, 0.30]


# Define eqn that is to be solved numerically
#Variables: p=partition coefficent, phi=vol fract of short, phib=vol frac big
#df=Free energy cost by pore for entering pore, Ns and Nb are number of subunits
# for short and big polymers

# Define Osmotic Pressure equation as function of short and big fraction
def osmotic_pressure(s, b):
	return s/Ns + b/Nb + (5./4)*alpha*(s + b)**(9./4)

# Let short polymer volume fraction range from 0 to 60 percent 
phis = np.linspace(0, 0.60)
# Let big polymer volume fraction range from 0 to 30 percent
phib30 = np.linspace(0, 0.30)
phib20 = np.linspace(0, 0.20)
phib10 = np.linspace(0, 0.10)


def f(p, phi, phib, df):
	""" Implicit definition of P(phi, dF) """
	#return -log(p) - df + (p-1)*phi + \
	#    ( phi*(1-p) + \
	#      5./4*alpha*(phi*p)**(5./4)*(p*phi-9./5) - \
	#      (phib/Nb)-5./4*alpha*(phi+phib)**(5./4)*(phi+phib-9./5) - \
	#      1./2*((1-p*phi)**2 - (1-phi-phib)**2) ) * Ns

	return -p + exp(- df + (p - 1)*phi + Ns*(phi*(1 - p) + phib + (5./4)*alpha*p**(5./4)*phi**(5./4)*(p*phi - (9./5))- \
		phib/Nb - (5./4)*alpha*(phi + phib)**(5./4)*(phi + phib-(9./5))-\
		(1./2)*((1 - p*phi)**(2) - (1 - phi - phib)**(2)))) 

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

pl.figure(1)

# Plot 0 to 30% big polymer and 0 short polymer
pl.semilogy(phib10, osmotic_pressure(0.0, phib10),'b,', label=' 0% PEG1000')
pl.semilogy(phib20, osmotic_pressure(0.0, phib20),'g,', label=' 0% PEG1000')
pl.semilogy(phib30, osmotic_pressure(0.0, phib30),'r,', label=' 0% PEG1000')

# Plot using whole range of phis and  10,20,30% big polymer
pl.semilogy(phis+0.10, osmotic_pressure(phis, 0.10),'b-', label='10% PEG3500')
pl.semilogy(phis+0.20, osmotic_pressure(phis, 0.20),'g-', label='20% PEG3500')
pl.semilogy(phis+0.30, osmotic_pressure(phis, 0.30),'r-', label='30% PEG3500')

#phi = np.linspace(0,1)
#pl.semilogy(phi, osmotic_pressure(0, phi), label='hi')

pl.legend(loc='lower right')
pl.xlabel('total polymer concentration')
pl.ylabel('Osmotic Pressure')
pl.title('Osmotic Pressure of Mixed Polymers')

pl.show()


pl.figure(2)
for i,phi_b in enumerate(phi_bs):
	phis = np.linspace(0, 1.0-phi_b, 100)
#	pl.subplot(len(phi_bs),1,i+1)
#	for df in dfs:
	ops= [osmotic_pressure(phi, phi_b) for phi in phis]
	for j,df in enumerate(dfs):
		try: ps = [phi*P(phi,phi_b,df) for phi in phis]
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
		pl.semilogx(ops,ps,color=colors[j],linestyle = lts[i],label=labels)#color = colors[j],
	pl.axis([0, 2.0, -0.1, 1.2])

#font = {'size': 22}
#matplotlib.rc('font',**font)
pl.text(0.02,0.68,r'$\Delta F = 0.0$', color='r')
pl.text(0.02,0.64,r'$\Delta F = 1.0$', color='g')
pl.text(0.02,0.59,r'$\Delta F = 2.0$', color='b')
pl.text(0.02,0.55,r'$\Delta F = 3.0$', color= 'k')

pl.legend(loc='upper left')
#pl.title('Partition Coefficient vs Osmotic Pressure for Varied $\Phi_{PEG400}$(Out) and Fixed $\Phi_{PEG3500}$')
pl.title('Number Fraction Of PEG400 In Pore for Varied $\Phi_{PEG400}$(Out) and Fixed $\Phi_{PEG3500}$')
pl.xlabel('Osmotic Pressure for Varied $\Phi_{PEG400}$(Out) and Fixed $\Phi_{PEG3500}$')
#pl.ylabel('Partition Coefficient P as a function of $\Phi_{PEG400}$(Out), $\Phi_{PEG3500}$, $\Delta$F')
pl.ylabel('Number Fraction $\Phi_{PEG400}$(In) as a function of $\Phi_{PEG400}$(Out), $\Phi_{PEG3500}$, $\Delta$F')

pl.savefig('110908_c_Part_Coeff_Mixed_VarySmall_Corrected.png')


pl.show()



