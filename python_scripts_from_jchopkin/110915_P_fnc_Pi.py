#!/usr/bin/python

import numpy as np
from numpy import log,exp
import scipy.optimize as opt
import matplotlib.pyplot as pl
import matplotlib

# Things will go wrong a lot in the code below
# Ignore most of the warnings
#np.seterr(all='ignore')

alpha = 0.49 # for PEG
Ns = 9
Nb = 79
dfs = [0.0, 0.2, 2.0, 10.0]
phi_ss = [0.000001, 0.3, 0.50, 0.70]
#Pi =  phi/Ns + phib/Nb + (5./4)*alpha*(phi + phib)**(9./4)
# Define equation to be solved:
# Variables are: p=partition coefficent; df=free energy cost by pore of polymer entering pore;
# phis=vol fract of short polymers<=>num fract; phib= vol(numb) fract big polymers;
# Nb,Ns= subunits in big,short polymers; Chi=1/2

#def Pi(phi,phib):
#	return phi/Ns + phib/Nb + (5./4)*alpha*(phi + phib)**(9./4)
def f(p, phi, phib, df):
	""" Implicit definition of P(phi,phib, Pi, dF) """
	return -p + exp(- df + (p - 1)*phi + Ns*(phi*(1 - p) + phib + (5./4)*alpha*p**(5./4)*phi**(5./4)*(p*phi - (9./5))- \
			phib/Nb - (5./4)*alpha*(phi + phib)**(5./4)*(phi + phib-(9./5))-\
			(1./2)*((1 - p*phi)**(2) - (1 - phi - phib)**(2))))
	##### this is for p as a function of pi
	#return -p + exp(- df + (p - 1)*phi + Ns*(phi*(1 - p) + phib + (5./4)*alpha*p**(9./4)*phi**(9./4) - (9./4)*alpha*p**(5./4)*phi**(5./4) - \
	#	(1./2)*(1 - p*phi)**2 - Pi + phi/Ns + (9./4)*alpha*(phi + phib)**(5./4) + (1./2)*(1 - phi - phib)**2)

#def Pi(phi,phib):
#	return phi/Ns + phib/Nb + (5./4)*alpha*(phi + phib)**(9./4)

# Use brent optimization algorithm if f has opposite signs
def P(phi, phib, df):
	""" Numerically solve for partition coefficient as a
	function of \phib """
	#print 'P'
	if f(0,phi,phib,df)*f(1,phi,phib,df) < 0:
		#print 'brent'
		#return opt.brentq(f, 0, 1, args=(phi,df)) # Brent's method
		#return opt.brenth(f, 0, 1, args=(phi,df)) # Brent's method
		return opt.bisect(f, 0, 1, args=(phi,phib,df)) # Bisection method
		#x,r = opt.bisect(f, 0, 1, args=(phi,df), full_output=True) # Bisection method
		#print r.iterations
		#return x
	else:
		#print 'newton'
		return opt.newton(f, 0.5, args=(phi,phib,df)) # Newton-Raphson
	#print 'newtonpre'
       	#return opt.newton(f, 1.0, args=(phi,phib,df)) # Newton-Raphson
	#print 'newtonpost'

# Give different fixed phi differnt line types via i loop (more solid =>larger phi)
# Give different delta f's different colors via j loop
lts = [':','-.','--','-']
colors = ['r','g','b','k','c','m','y']

# Plot P as a function of phi big
#pl.figure()
#print 'here'
for i,phi_s in enumerate(phi_ss):
	pi = np.linspace(0, 1.0,100)
	phibb = np.linspace(0, 1.0-phi_s, 100)
	#def Pi(phi_s,phib):  
	#	return phi_s/Ns + phib/Nb + (5./4)*alpha*(phi_s + phib)**(9./4)
#	#Pis = [Pi(phi_s, phib) for phib in phibb]
	#def AntiPi(phi_s,phib):
	#	return  phib/Nb + (5./4)*alpha*(phi_s + phib)**(9./4)
#	AntiPis = [AntiPi(phi_s,phib) for phib in phibb]
#	pl.subplot(len(phi_ss),1,i+1) # use this if you want to see all phi separately
#	for df in dfs:
	for j,df in enumerate(dfs):
#		print 'j=%g'%j
#		print 'i=%g'%phi_s
		try: 
			#print 'try'
			#print Pi
			ps = [phi_s*P(phi_s,phib,df) - Pi(phi_s,phib) + AntiPi(phi_s,phib) + phi_s/Nb for phib in phibb]
			#print ps
	       	except: 
			#print 'except'
			continue
		#print 'try & exc r ok'
		#resids = [abs(f(p/phi_s,phi_s,phib,df)) for p,phib in zip(ps, phibb)]
		#max_resid = max(resids)
		#print 'Largest residual for df=%f: %e' % (df, max_resid)
		if j==3:
			labels=r'$\Phi_{PEG400} = %.2f$' % phi_s
		else:
			labels=None
			###########################################################
		pl.plot(Pi,ps, color=colors[j],linestyle = lts[i],label=labels)
	#	pl.plot(Pi,- Pi + AntiPi, color=colors[j],linestyle = lts[i],label=labels)
pl.axis([0.0, 1.3, -0.1, 1.0])


pl.text(0.04,0.88,r'$\Delta F = 0.0$', color= 'r')
pl.text(0.04,0.84,r'$\Delta F = 0.2$', color= 'g')
pl.text(0.04,0.80,r'$\Delta F = 2.0$', color= 'b')
pl.text(0.04,0.76,r'$\Delta F = 10.0$', color= 'k')

font = {'size':18}
matplotlib.rc('font',**font)

pl.legend(loc='center right')
pl.title('Number Fraction of PEG400 In Pore for Fixed  $\Phi_{PEG400}$ and Varied $\Phi_{PEG3500}$')
pl.xlabel('Number fraction of PEG3500 $\Phi_{PEG3500}$')
pl.ylabel('Number Fraction $\Phi_{PEG400}$(In) as a function of $\Phi_{PEG400}$(Out), $\Phi_{PEG3500}$, $\Delta$F)')

pl.savefig('110915_P_fnc_Pi_VarySmall.png')

pl.show()



