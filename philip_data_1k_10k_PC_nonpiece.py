#!/usr/bin/python

import numpy as np
from numpy import log,exp
import scipy.optimize as opt
import matplotlib.pyplot as pl
import matplotlib
from scipy import interpolate

alpha = 0.49
Ns = 23
Nb = 230
dfs = [0.0,2.0, 4.0]
phi_bs = [0.05, 0.10, 0.15]
phi_water_min = 0.01
gs_PEG200=[0.835,0.724,0.609,0.499,0.402] #experimental conducatance of pure peg200
gs_PEG1k_add= [0.720,0.574,0.488,0.407,0.278] #experimental conductance of peg200 added to peg3400
percs=[3,6,9,12,15,] #percent of peg200 (either added-to or pure)

# Define eqn that is to be solved numerically
def f(p, phi, phib, df):
	""" Implicit definition of P(phi, dF) """
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

#Partition Coeff from experimental data
def g(a,b):
	return (0.939-a)/(0.939 - b) # g_max = 0.939

lts = [':','--','-']
markers = ['x', '^', 'd']
#colors = ['k','r','b','g','c','m','y']

pl.figure()
#Plot experimental data points and interpolated curve
for k,perc in enumerate(percs):
	partition = [g(a,b) for a,b in zip(gs_PEG1k_add,gs_PEG200)]
	print perc,partition[k]
	pl.plot(0.01*perc,partition[k],'b',marker= 's',linestyle = '--', linewidth= 2.5, markersize= 4, markevery= 1)
	
	#interpolate experimental data points
	data_perc =[.03,.06,.09,.12,.15] #percents of peg200 added
	data_conduct =[2.106,1.70,1.367,1.209,1.231] #conductance of peg200added
	z = interpolate.interp1d(data_perc,data_conduct) 
	data_perc_new = np.arange(0.03,.15,.001) #refined percent points
	data_conduct_new = z(data_perc_new) #new interpolated conductance points
	pl.plot(data_perc, data_conduct, 'o', data_perc_new, data_conduct_new, '-')

#Plot theory partition coeff curves for different fixed vals of deltaF and big polymer vol frac
for i,phi_b in enumerate(phi_bs):
	phis = np.linspace(0.0, 1.0-phi_b-phi_water_min, 100)
#	for df in dfs:
	for j,df in enumerate(dfs):
		try: ps = [P(phi,phi_b,df) for phi in phis]
	       	except: continue
		resids = [abs(f(p,phi,phi_b,df)) for p,phi in zip(ps, phis)]
		max_resid = max(resids)
		print 'Largest residual for df=%f: %e' % (df, max_resid)
		if j==0:
			labels=r'$\phi_{b} = %.2f$' % phi_b
		else:
			labels=None
		if i==2:
			labelss=r'$\Delta f = %.1f$' % df
		else:
			labelss=None
		#pps = [phi*P(phi,phi_b,df) for phi in phis]
		print 'Last phi tried for df=%f: phis=%.2f: phib=%.1f: total phi=%.2f' % (df, phi,phi_b,phi+phi_b)
		pl.plot(phis,ps,'k',marker= markers[j],linestyle = lts[i], linewidth= 0.5, markersize= 4, markevery= 10, markerfacecolor= 'w')#label=labels)#color = colors[j],
		#pl.plot(3,2.16,'b','-')
		#pl.plot(6,1.70,'b','-')
		pl.plot(-2.0,-2.0,'k',linestyle = lts[i], linewidth= 1.0,label=labels)#color = colors[j],
		pl.plot(-1.0,-1.0,'k',marker= markers[j],linestyle = 'none', markersize= 6,markerfacecolor= 'w',label=labelss)#color = colors[j],
		
#		pl.plot(perc,partition[k],'b',marker= 's',linestyle = '-', linewidth= 0.5, markersize= 4, markevery= 1)
		#pl.plot(exptphis,exptp,'r', markersize= 10,markerfacecolor= 'r')#color = colors[j],
#for k,perc in enumerate(percs):
#	partition = [g(a,b) for a,b in zip(gs_PEG1k_add,gs_PEG200)]
#	print perc,partition[k]
#	pl.plot(perc,partition[k],'b',marker= 's',linestyle = '-', linewidth= 0.5, markersize= 4, markevery= 1)

pl.text(0.002,1.35,'Exper. Data', color='b')
pl.text(0.002,1.24,'w/ 15% PEG10k', color='b')
pl.legend(loc='lower right')
pl.title('Partition Coefficient for Varied $\Phi_{PEG1,000}$ and Fixed $\Phi_{PEG10,000}$')
pl.xlabel(r'$\phi_{PEG1k}$', size = 'x-large')
pl.ylabel(r'$p  ( \phi_{PEG1k}  , \phi_{PEG10k} , \Delta f )$', size = 'x-large')
pl.axis([0.0, 0.50 ,0.0, 3.0])
pl.savefig('120618_PhilipData_1k_10k_PC_cspline1d.eps', dpi=600)


pl.show()







