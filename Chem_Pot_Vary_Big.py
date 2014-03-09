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
Ns = 9.
Nb = 79.
dfs = [0.0, 0.5, 2.0]
#dfs = [0.0,-5.0, -10.0] #minimized content for poster
phi_bs = [0.05, 0.10, 0.15]
#phi_bs = [0.10, 0.30]# for poster
phi_water_min = 0.01

def f(p, phi, phib, df):
	""" Implicit definition of P(phi, dF) """
	return - p + exp( - df + Ns*(log((1 - p*phi)/(1 - phi - phib)) + \
		(p - 1)*phi - phib + (9./4)*alpha*((phi + phib)**(5./4) - (p*phi)**(5./4))))

def Mu_s_w(phi,phib): 
	return log(phi) + 1 - Ns*(log(1-phi-phib)+1-(1-phi-phib)) + Ns*(9./4)*alpha*(phi + phib)**(5./4) + \
		Ns*(log(1 - phi - phib) + 1 - (1 - phi - phib) - phi/Ns - phib/Nb + \
		(1./2)*(1 - (1 - phi - phib))**2 + (5./4)*alpha* \
		(phi + phib)**(9./4))

def term(phi,phib):
	return (9./4)*alpha*Ns*(phi + phib)**(5./4)

def Mu_s(phi,phib): 
	return phi + exp(1 - phi - (1-phi-phib)*Ns - phib*(Ns / Nb) +  \
	(1./2)*Ns*(1-phi-phib)**(2) - (5./4)*alpha*Ns*(phi + phib)**(9./4) + \
	(9./4)*alpha*Ns*(phi + phib)**(5./4))
	
#Optimization algorithms, if opposite sign use brent, else use newton
def P(phi, phib, df):
	""" Numerically solve for partition coefficient as a
	    function of \phi_s """
	if f(0,phi,phib,df)*f(1,phi,phib,df) < 0:
		return opt.bisect(f, 0, 1, args=(phi,phib,df), maxiter=500) # Bisection method
	else:
		return opt.newton(f, 1.0, args=(phi,phib,df), maxiter=5000) # Newton-Raphson

lts = [':','-.','--','-']
colors = ['k','r','b','g','c','m','y']

pl.figure()
for i,phi_b in enumerate(phi_bs):
	phis = np.linspace(0.01, 1.0-phi_b-phi_water_min, 100)
	x = np.linspace(0.01, 20, 100)
	#phis = np.linspace(0.0, 0.5, 100)
	#phis10 = np.linspace(0.11, 1.0-phi_water_min, 100)
	#phis20 = np.linspace(0.21, 1.0-phi_water_min, 100)
	for j,df in enumerate(dfs):
		try: ps = [P(phi,phi_b,df) for phi in phis]
	       	except: continue
		resids = [abs(f(p,phi,phi_b,df)) for p,phi in zip(ps, phis)]
		max_resid = max(resids)
		print 'Largest residual for df=%f: %e' % (df, max_resid)
		if j==0:
			labels=r'$\Phi_{PEG3500} = %.2f$' % phi_b
		else:
			labels=None
		#Mu_s_w_mix10 = [Mu_s_w(phi,10)/Mu_s_w(phi10,0.1) for phi,phi10 in zip(phis,phis10)]
		#Mu_s_w_mix20 = [Mu_s_w(phi,20) for phi in phis]
		Mu_s_w_mix = [Mu_s_w(phi,phi_b) for phi in phis]
		Mu_s_w_pure = [Mu_s_w(phi+phi_b,0) for phi in phis]
		#Mu_s_w_pure10 = [Mu_s_w(phi10,0) for phi10 in phis10]
		#Mu_s_w_pure20 = [Mu_s_w(phi20,0) for phi20 in phis20]
		#phi_pure = [phi + phi_b for phi in phis]
		
		#Mu_s_w_ratio = [Mu_s_w_mix/Mu_s_w_pure for phi in phis]
		#Mu_s_ratio = [Mu_s(phi,phi_b)/Mu_s(phi+phi_b,0) for phi in phis]
		
		#Mu_s_ratio = [Mu_s(phi,phi_b)/Mu_s(phi,0) for phi in phis]
		print 'Last phi tried for df=%f: phis=%.2f: phib=%.1f: total phi=%.2f' % (df, phi,phi_b,phi+phi_b)
		
		##pl.plot(phis+phi_b,Mu_s_w_ratio,color=colors[0],linestyle = lts[i],label=labels)
		pl.plot(Mu_s_w_pure,Mu_s_w_mix,color=colors[j],linestyle = lts[i],label=labels)
		pl.plot(x,x,color=colors[0],linestyle = lts[3])#,label=labels)
		
		####pl.plot(phis+phi_b,Mu_s_w_pure,color=colors[0])#,linestyle = lts[i],label=labels)
		####pl.plot(phis,Mu_s_w_mix,color=colors[1],linestyle = lts[i],label=labels)
		
		#pl.plot(phis,Mu_s_w_mix/Mu_s_w_pure,color=colors[2],linestyle = lts[i],label=labels)
		#pl.plot(phis20,Mu_s_w_pure20,color=colors[2],linestyle = lts[i],label=labels)
		#pl.plot(phis,Mu_s_w_mix,color=colors[3],linestyle = lts[i],label=labels)
		
		#pl.plot(phis+phi_b,Mu_s_w_ratio,color=colors[2],linestyle = lts[i],label=labels)
		
	#pl.axis([0, 1.0, 0.1, 2.4])
	pl.axis([0, 10.0, 0, 10.0])


pl.text(0.1,6.5,'1:1 fuduciary line', color='k')
pl.text(0.1,7.00,'$\mu_{mix}$/$\mu_{pure}$ for different $\phi_{bigs}$', color='b')

pl.legend(loc='upper left')
pl.title('$\mu_{PEG400}$(mixed) versus $\mu_{PEG400}$(pure) for Varied $\Phi_{PEG400}$ and Fixed $\Phi_{PEG3500}$')
pl.xlabel('$\mu_{PEG400}$(pure)')#otal Number fraction of PEG Monomers in Bath (big and small)')
pl.ylabel('$\mu_{mixed}$(mixed)')

pl.savefig('Chem_Pot_VaryShort', dpi=600)


pl.show()





