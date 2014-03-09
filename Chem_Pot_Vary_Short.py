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
phi_bs = [0.05, 0.10, 0.15]
phi_water_min = 0.01

def Mu_s_w(phi,phib): 
	return log(phi) + 1 - Ns*(log(1-phi-phib)+1-(1-phi-phib)) + Ns*(9./4)*alpha*(phi + phib)**(5./4) + \
		Ns*(log(1 - phi - phib) + 1 - (1 - phi - phib) - phi/Ns - phib/Nb + \
		(1./2)*(1 - (1 - phi - phib))**2 + (5./4)*alpha* \
		(phi + phib)**(9./4))

def Mu_s(phi,phib): 
	return phi + exp(1 - phi - (1-phi-phib)*Ns - phib*(Ns / Nb) +  \
	(1./2)*Ns*(1-phi-phib)**(2) - (5./4)*alpha*Ns*(phi + phib)**(9./4) + \
	(9./4)*alpha*Ns*(phi + phib)**(5./4))
	

phis = np.linspace(0.001, 0.7,200)
phib = np.linspace(0.1, 1.0,200)
phib30 = np.linspace(0.001, 0.30,30)
phib20 = np.linspace(0.001, 0.20,20)
phib10 = np.linspace(0.001, 0.10,10)

phis10diff = np.linspace(0.0, 0.9,200)
phib10diff = np.linspace(0.1, 1.0,200)


phis20diff = np.linspace(0.0, 0.8,200)
phib20diff = np.linspace(0.2, 1.0,200)


phis30diff = np.linspace(0.0, 0.7,200)
phib30diff = np.linspace(0.3, 1.0,200)

pl.figure(1)
#
#pl.plot(phib10, Mu_s_w(0.001, phib10),'b,', label=' 0% PEG400')
#pl.plot(phib20, Mu_s_w(0.001, phib20),'g,', label=' 0% PEG400')
#pl.plot(phib30, Mu_s_w(0.001, phib30),'r,', label=' 0% PEG400')
#
## Plot using whole range of phis and  10,20,30% big polymer
#pl.plot(phis+0.10, Mu_s_w(phis, 0.10),'b-', label='10% PEG3500')
#pl.plot(phis+0.20, Mu_s_w(phis, 0.20),'g-', label='20% PEG3500')
#pl.plot(phis+0.30, Mu_s_w(phis, 0.30),'r-', label='30% PEG3500')
#
#pl.plot(phib10, Mu_s(0.001, phib10),'b,', label=' 0% PEG400')
#pl.plot(phib20, Mu_s(0.001, phib20),'g,', label=' 0% PEG400')    
#pl.plot(phib30, Mu_s(0.001, phib30),'r,', label=' 0% PEG400')   
#                                                              
## Plot using whole range of phis and  10,20,30% big polymer   
#pl.plot(phis+0.10, Mu_s(phis, 0.10),'b-', label='10% PEG3500')
#pl.plot(phis+0.20, Mu_s(phis, 0.20),'g-', label='20% PEG3500') 
#pl.plot(phis+0.30, Mu_s(phis, 0.30),'r-', label='30% PEG3500') 


#pl.semilogx(phib10, Mu_s_w(0.001, phib10),'k,', label=' 0% PEG400')
#pl.semilogx(phib20, Mu_s_w(0.001, phib20),'g,', label=' 0% PEG400')
#pl.semilogx(phib30, Mu_s_w(0.001, phib30),'r,', label=' 0% PEG400')
#pl.semilogx(phis+0.10, Mu_s_w(phis, 0.10),'b-', label='10% PEG3500')
#pl.semilogx(phis+0.20, Mu_s_w(phis, 0.20),'g-', label='20% PEG3500')
#pl.semilogx(phis+0.30, Mu_s_w(phis, 0.30),'r-', label='30% PEG3500')
pl.semilogy(phib10, Mu_s(0.001, phib10),'b,', label=' 0% PEG400')
pl.semilogy(phib20, Mu_s(0.001, phib20),'g,', label=' 0% PEG400')    
pl.semilogy(phib30, Mu_s(0.001, phib30),'r,', label=' 0% PEG400')   
pl.semilogy(phib, Mu_s(0.001, phib),'k,', label=' 0% PEG400')   
pl.semilogy(phis+0.10, Mu_s(phis, 0.10),'b-', label='10% PEG3500')
pl.semilogy(phis+0.20, Mu_s(phis, 0.20),'g-', label='20% PEG3500') 
pl.semilogy(phis+0.30, Mu_s(phis, 0.30),'r-', label='30% PEG3500') 

#pl.semilogx(phib10, Mu_s_w(0.001, phib10),'k,', label=' 0% PEG400')
#pl.semilogx(phib20, Mu_s_w(0.001, phib20),'g,', label=' 0% PEG400')
#pl.semilogx(phib30, Mu_s_w(0.001, phib30),'r,', label=' 0% PEG400')
#pl.semilogx(phis+0.10, Mu_s_w(phis, 0.10),'b-', label='10% PEG3500')
#pl.semilogx(phis+0.20, Mu_s_w(phis, 0.20),'g-', label='20% PEG3500')
#pl.semilogx(phis+0.30, Mu_s_w(phis, 0.30),'r-', label='30% PEG3500')

#pl.semilogy(phib10, Mu_s(0.001, phib10),'b,', label=' 0% PEG400')
#pl.semilogy(phib20, Mu_s(0.001, phib20),'g,', label=' 0% PEG400')   
#pl.semilogy(phib30, Mu_s(0.001, phib30),'r,', label=' 0% PEG400')   
#pl.semilogy(phis+0.10, Mu_s(phis, 0.10),'b-', label='10% PEG3500')
#pl.semilogy(phis+0.20, Mu_s(phis, 0.20),'g-', label='20% PEG3500') 
#pl.semilogy(phis+0.30, Mu_s(phis, 0.30),'r-', label='30% PEG3500') 

pl.axis([0.09, 0.5, 0.04, 1.6])

pl.figure(2)

pl.plot(phib10diff, Mu_s(phis10diff, 0.10)-Mu_s(0.0, phib10diff),'b', label=' 10% PEG3500')
pl.plot(phib20diff, Mu_s(phis20diff, 0.20)-Mu_s(0.0, phib20diff),'g', label=' 20% PEG3500')
pl.plot(phib30diff, Mu_s(phis30diff, 0.30)-Mu_s(0.0, phib30diff),'r', label=' 30% PEG3500')
pl.axis([0.1, 0.5, 0.0, 0.2])
#pl.semilogy(phib20, Mu_s(0.001, phib20),'g,', label=' 0% PEG400')    
#pl.semilogy(phib30, Mu_s(0.001, phib30),'r,', label=' 0% PEG400')   
#pl.semilogy(phib, Mu_s(0.001, phib),'k,', label=' 0% PEG400')   

pl.show()

pl.legend(loc='upper left')
#pl.title('Chemical Potential of Small Polymers $\mu_{PEG400}$')
#pl.xlabel('Total Polymer Volume Fraction in Bath')#otal Number fraction of PEG Monomers in Bath (big and small)')
#pl.ylabel('$\mu_{PEG400}$')
pl.savefig('Chem_Pot_Vary_Short.eps', dpi=600)
pl.savefig('Chem_Pot_Vary_Short.pdf', dpi=600)

