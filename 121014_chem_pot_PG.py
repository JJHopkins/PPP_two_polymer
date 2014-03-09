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
Ns = 5
Nb = 77
#dfs = [0.0, 2.0, 4.0]
dfs = [-10.0, -5.0, 0.0,5.0,10.0]
phi_bs = [0.0, 0.15]
#phi_bs = [0.10, 0.30]# for poster
phi_water_min = 0.01

# Define eqn that is to be solved numerically
#Variables: p=partition coefficent, phi=vol fract of short, phib=vol frac big
#df=Free energy cost by pore for entering pore, Ns and Nb are number of subunits121014_chem_pot_PG.py
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
#
#def Mu_s(phi,phib): #this is with the chemical pot of water s.t. Ns*mu_w(in)=Ns*mu-w(out), check this 121015 
#	return log(phi) + 1 - Ns*(log(1-phi-phib)+1-(1-phi-phib)) + Ns*(9./4)*alpha*(phi + phib)**(5./4) + \
#		Ns*(log(1 - phi - phib) + 1 - (1 - phi - phib) - phi/Ns - phib/Nb + \
#		(1./2)*(1 - (1 - phi - phib))**2 + (5./4)*alpha* \
#		(phi + phib)**(9./4) - log(df))
def Mu_s(phi,phib): 
	return phi + exp(1 - phi - (1-phi-phib)*Ns - phib*(Ns / Nb) +  \
	(1./2)*Ns*(1-phi-phib)**(2) - (5./4)*alpha*Ns*(phi + phib)**(9./4) + \
	(9./4)*alpha*Ns*(phi + phib)**(5./4) - df)

#lts = [':','-.','--','-']
lts = ['--','-']
#markers = ['x', '^', 'd']
colors = ['b','g','k','y','r']

pl.figure()
for i,phi_b in enumerate(phi_bs):
	phis = np.linspace(0.0, 1.0-phi_b-phi_water_min, 100)
#	pl.subplot(len(phi_bs),1,i+1)
#	for df in dfs:
	for j,df in enumerate(dfs):
		try: ps = [P(phi,phi_b,df) for phi in phis]
	       	except: print 'hello';continue
		resids = [abs(f(p,phi,phi_b,df)) for p,phi in zip(ps, phis)]
		max_resid = max(resids)
		print 'Largest residual for df=%f: %e' % (df, max_resid)
		if j==2:
			labels=r'$\phi_{b} = %.2f$' % phi_b
		else:
			labels=None
		if i==0:
			labelss=r'$\Delta f = %.1f$' % df
		else:
			labelss=None
		#pps = [phi*P(phi,phi_b,df) for phi in phis]
	#	print 'Last phi tried for df=%f: phis=%.2f: phib=%.1f: total phi=%.2f' % (df, phi,phi_b,phi+phi_b)
	#	pl.plot(phis,ps,'k',marker= markers[j],linestyle = lts[i], linewidth= 0.5, markersize= 4, markevery= 7, markerfacecolor= 'w')#label=labels)#color = colors[j],
	#	pl.plot(-2.0,-2.0,'k',linestyle = lts[i], linewidth= 1.0,label=labels)#color = colors[j],
	#	pl.plot(-1.0,-1.0,'k',marker= markers[j],linestyle = 'none', markersize= 6,markerfacecolor= 'w',label=labelss)#color = colors[j],

		Mu_s_for_plot = log(Mu_s(phis,phi_b))

	#pl.plot(phis,Mu_s_for_plot) 
		pl.plot(Mu_s_for_plot,phis,color = colors[j],linestyle = lts[i], linewidth= 1.5, label=labels) 
		#pl.text(0.9,1.95,r'$\Delta F = %.2f$' % j(0), color='r')
		#pl.text(0.9,1.85,r'$\Delta F = %.2f$' % j(1), color='g')
		#pl.text(0.9,1.75,r'$\Delta F = %.2f$' % j(2), color='b')
		#pl.text(0.9,1.65,r'$\Delta F = %.2f$' % j(3), color= 'k')


#	pl.axis([0, 1.0, -0.20, 6.0])

#font = {'size': 22}
#matplotlib.rc('font',**font)

#pl.text(5.0,0.22,r'$\Delta f = 0.0$', color='b')
#pl.text(5.0,0.17,r'$\Delta f = 2.0$', color='g')
#pl.text(5.0,0.12,r'$\Delta f = 4.0$', color='r')

pl.text(2700,0.15,r'$\Delta F = -10.0$', color='b')
pl.text(2700,0.11,r'$\Delta F = -5.0$', color='g')
pl.text(2700,0.07,r'$\Delta F =\, \,0.0$', color='k')
pl.text(2700,0.03,r'$\Delta F = \,\, 5.0$', color='y')
pl.text(2700,-0.01,r'$\Delta F = \,\,10.0$', color= 'r')

#pl.text(0.85,3.60,r'$\Delta f = 0.0$', color='k')
#pl.text(0.85,3.30,r'$\Delta f = 1.0$', color='r')
#pl.text(0.85,3.00,r'$\Delta f = 2.5$', color='b')
#pl.text(0.85,2.70,r'$\Delta f = 5.0$', color='g')
#pl.text(0.85,2.40,r'$\Delta f = 10.0$', color='c')
#pl.text(0.85,2.10,r'$\Delta f = 15.0$', color= 'm')

#pl.legend(loc='lower left')
#pl.title('Partition Coefficient for Varied $\Phi_{PEG400}$(Out) and Fixed $\Phi_{PEG3500}$(Out)')
pl.title('$\Phi_{200}$(In) for Baths with 0% or 15% $\Phi_{3400}$ vs. $\mu_{s}( \phi_{s}  , \phi_{b} , \Delta f )$')
pl.xlabel(r'$\mu_{s}( \phi_{s}  , \phi_{b} , \Delta f )$', size = 'x-large')
#pl.xlabel(r'$\phi_{s} + \phi_{b}$', size = 'x-large')
pl.ylabel(r'$\phi_{s}$', size = 'x-large')
#pl.ylabel(r'$p  ( \phi_{s}  , \phi_{b} , \Delta f )$', size = 'x-large')
#pl.axis([-0.019, 0.9, -0.025, 0.9])

pl.savefig('121022_121014_chem_pot_PG_logMu_s.eps', dpi=600)


pl.show()




