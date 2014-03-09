#!/usr/bin/python

import numpy as np
from numpy import log,exp
import scipy.optimize as opt
import matplotlib.pyplot as pl
import matplotlib
from scipy import interpolate
from scipy.interpolate import UnivariateSpline
# Things will go wrong a lot in the code below
# Ignore most of the warnings
#np.seterr(all='ignore')

alpha = 0.49
Ns = 5
Nb = 77
dfs = [0.0,0.2,0.4,0.5,0.6]
#dfs = [-0.5,-0.3,-0.2,-0.1,0.0]
#dfs = [0.0,-5.0, -10.0] #minimized content for poster
phi_bs = [0.15, 0.15, 0.15]
#phi_bs = [0.05, 0.10, 0.15]
#phi_bs = [0.10, 0.30]# for poster
phi_water_min = 0.01
#exptphis = [0.03,0.06,0.09,0.12,0.15]
#exptp = [1.299,1.306,1.198,1.171,0.955]
# I duno where these came from, they're wrong, need to check which values I'm asking for, checked, ones below are ok...
#gs_PEG200=[0.835,0.724,0.609,0.499,0.402]
#gs_PEG200_add= [0.794,0.660,0.536,0.425,0.326]

percs=[3,6,9,12,15,]
#spline
gs_PEG200_s=      [0.819,0.699,0.594,0.504,0.429]
gs_PEG200_add_s = [0.794,0.659,0.536,0.425,0.325]
#interpolated
gs_PEG200=      [ 0.826,0.719,0.622,0.511,0.392]
gs_PEG200_add = [ 0.792,0.666,0.527,0.431,0.324]
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

def g(a,b):
	return (0.939-a)/(0.939 - b)
lts = [':','--','-']
#lts = ['--','-']#for poster
markers = ['x', '^', 'd','s','*']
#colors = ['k','r','b','g','c','m','y']

pl.figure()

for k,perc in enumerate(percs):
	partition = [g(a,b) for a,b in zip(gs_PEG200_add,gs_PEG200)]
	partition_s = [g(a,b) for a,b in zip(gs_PEG200_add_s,gs_PEG200_s)]
	
	print perc,partition[k]
	print perc,partition_s[k]
	
	pl.plot(0.01*perc,partition[k],'r',marker= 's',linestyle = '--', linewidth= 2.5, markersize= 4, markevery= 1)
	pl.plot(0.01*perc,partition_s[k],'b',marker= 's',linestyle = '--', linewidth= 2.5, markersize= 4, markevery= 1)
	
	y=[.03,.06,.09,.12,.15]
	z=[1.394,1.298,1.221,1.168,1.142]
	x = interpolate.UnivariateSpline(y,z,k=3)
	xx = interpolate.UnivariateSpline(y,z,k=2)
	yy=np.arange(0,.16,.01)
	zz=x(yy)
	zzz=xx(yy)
#	if k==0:
#		labels='3rd fit'
#	else:
#		labels=None
#	if k==1:
#		labelss='2nd fit' 
#	else:
#		labelss=None
	coeff = UnivariateSpline.get_coeffs(xx)
	#fit_params=polyfit(yy,zz), finish this, see chat with james on 6/19
	#evals = UnivariateSpline._call_(y,nu=0)
#	print coeff
	#print 'Derivs: %,% ' (x,deriv)
	#pl.plot(yy,zz,'b', label= labels)
	#pl.plot(yy,zz+1.4,'b', label= labels)
#	pl.plot(yy,zzz,'r',linestyle=':', linewidth=0.2)#label= labelss) 
#	pl.plot(yy,zzz+1.0,'r',linestyle=':', linewidth=0.2)#label= labelss) 
#	pl.plot(yy,zzz+0.8,'r',linestyle=':', linewidth=0.2)#label= labelss) 
#	pl.plot(yy,zzz+0.6,'r',linestyle=':', linewidth=0.2)#label= labelss) 
#	pl.plot(yy,zzz+0.4,'r',linestyle=':', linewidth=0.2)#label= labelss) 
#	pl.plot(yy,zzz+0.2,'r',linestyle=':', linewidth=0.2)#label= labelss) 

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
		pl.plot(-2.0,-2.0,'k',linestyle = lts[i], linewidth= 1.0,label=labels)#color = colors[j],
		pl.plot(-1.0,-1.0,'k',marker= markers[j],linestyle = 'none', markersize= 6,markerfacecolor= 'w',label=labelss)#color = colors[j],

		#pl.plot(exptphis,exptp,'r', markersize= 10,markerfacecolor= 'r')#color = colors[j],

#pl.text(0.03,1.30,'Data w/', color='r')
#pl.text(0.03,1.20,'15% PEG3400', color='r')
#pl.text(0.002,1.0,'Coeff:1.39,1.17,1.14', color='r')
pl.legend(loc='lower right')
pl.title('Partition Coefficient for Varied $\Phi_{PEG200}$ and Fixed $\Phi_{PEG3400}$')
pl.xlabel(r'$\phi_{200}$', size = 'x-large')
pl.ylabel(r'$p  ( \phi_{200}  , \phi_{3,400} , \Delta f )$', size = 'x-large')
pl.axis([0.0, 0.5 ,0.0,2.0])
pl.savefig('120619_philip_data_200_3400_PC_spline2_shifted.eps', dpi=750)#, dpi=600)


pl.show()





