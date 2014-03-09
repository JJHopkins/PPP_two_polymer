#!/usr/bin/python

import numpy as np
from numpy import log,exp
#from fitting import *
from scipy.optimize import curve_fit
import matplotlib.pyplot as pl
import matplotlib

alpha = 0.49
Ns = 23
Nb = 227
dfs = [0.0, 4.0, 9.0]
phi_bs = [0.05, 0.10, 0.15]
phi_water_min = 0.01
exptphis = [0.18,0.21,0.24,0.28,0.30]
exptp = [1.17,1.18,1.03,1.04,1.12]

#Fit function
def func(x,a,b,c):
	return (a*x**2 + b*x + c )
#Data
x = [0.0,2.0,5.0,10.0,15.0,20.0,25.0,30.0]
y = [0.939, 0.864,0.751,0.590,0.392,0.295,0.228,0.206]
yerr = [0.019,0.015,0.018,0.013,0.012, 0.011, 0.026,0.005]

p0 = [1.0,3.0,-1.0] #parameter ansantz
popt,punc,rc,d = curve_fit(func,x,y,p0,yerr) # rc is reduced chi, d is degrees of freedom
print 'optimal parameters: ', popt
print 'uncertainties of parameters: ', punc


#Optimization algorithms, if opposite sign use brent, else use newton
#def P(phi, phib, df):
#	""" Numerically solve for partition coefficient as a
#	    function of \phi_s """
#	if f(0,phi,phib,df)*f(1,phi,phib,df) < 0:
#		#print 'bisect'
#		#return opt.brentq(f, 0, 1, args=(phi,df)) # Brent's method
#		#return opt.brenth(f, 0, 1, args=(phi,df)) # Brent's method
#		return opt.bisect(f, 0, 1, args=(phi,phib,df), maxiter=500) # Bisection method
#		#x,r = opt.bisect(f, 0, 1, args=(phi,df), full_output=True) # Bisection method
#		#print r.iterations
#		#return x
#	else:
#		#print 'newton'
#		return opt.newton(f, 1.0, args=(phi,phib,df), maxiter=5000) # Newton-Raphson
#		#return opt.bisect(f, 0, 1, args=(phi,phib,df)) # Bisection method
#		# linestyles for different phi_big, colors for different deltaf's
#
#lts = [':','--','-']
##lts = ['--','-']#for poster
#markers = ['x', '^', 'd']
##colors = ['k','r','b','g','c','m','y']
#
#pl.figure()
#for i,phi_b in enumerate(phi_bs):
#	phis = np.linspace(0.0, 1.0-phi_b-phi_water_min, 100)
##	pl.subplot(len(phi_bs),1,i+1)
##	for df in dfs:
#	for j,df in enumerate(dfs):
#		try: ps = [P(phi,phi_b,df) for phi in phis]
#	       	except: continue
#		resids = [abs(f(p,phi,phi_b,df)) for p,phi in zip(ps, phis)]
#		max_resid = max(resids)
#		print 'Largest residual for df=%f: %e' % (df, max_resid)
#		if j==0:
#			labels=r'$\phi_{b} = %.2f$' % phi_b
#		else:
#			labels=None
#		if i==2:
#			labelss=r'$\Delta f = %.1f$' % df
#		else:
#			labelss=None
#		#pps = [phi*P(phi,phi_b,df) for phi in phis]
#		print 'Last phi tried for df=%f: phis=%.2f: phib=%.1f: total phi=%.2f' % (df, phi,phi_b,phi+phi_b)
#		pl.plot(phis,ps,'k',marker= markers[j],linestyle = lts[i], linewidth= 0.5, markersize= 4, markevery= 7, markerfacecolor= 'w')#label=labels)#color = colors[j],
#		pl.plot(-2.0,-2.0,'k',linestyle = lts[i], linewidth= 1.0,label=labels)#color = colors[j],
#		pl.plot(-1.0,-1.0,'k',marker= markers[j],linestyle = 'none', markersize= 6,markerfacecolor= 'w',label=labelss)#color = colors[j],
#
pl.plot(x,y,'r', markersize= 10,markerfacecolor= 'r')#color = colors[j],
		#pl.text(0.9,1.95,r'$\Delta F = %.2f$' % j(0), color='r')
		#pl.text(0.9,1.85,r'$\Delta F = %.2f$' % j(1), color='g')
		#pl.text(0.9,1.75,r'$\Delta F = %.2f$' % j(2), color='b')
		#pl.text(0.9,1.65,r'$\Delta F = %.2f$' % j(3), color= 'k')

pl.legend(loc='lower right')
#pl.title('Partition Coefficient for Varied $\Phi_{PEG400}$(Out) and Fixed $\Phi_{PEG3500}$(Out)')
pl.xlabel(r'$\phi_{s} + \phi_{b}$', size = 'x-large')
pl.ylabel(r'$p  ( \phi_{s}  , \phi_{b} , \Delta f )$', size = 'x-large')
pl.axis([0.0, 1.0, -0.05, 4])
pl.savefig('120608_conductance_200_3400.eps', dpi=600)


pl.show()





