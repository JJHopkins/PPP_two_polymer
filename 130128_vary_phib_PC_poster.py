#!/usr/bin/python
#testing changes for git

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
dfs = [0.0, 3.0, 5.0, 9.0] #2.5, 5.0, 10.0, 15.0]
#dfs = [0.0, -5.0, -10.0] #minimized content for poster
phi_ss = [0.10, 0.20]
#phi_ss = [0.10, 0.30] #minimized content for poster
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
		#return opt.brentq(f, 0, 1, args=(phi,phib,df)) # Brent's method
		#return opt.brenth(f, 0, 1, args=(phi,phib,df)) # Brent's method
		return opt.bisect(f, 0, 1, args=(phi,phib,df), maxiter=500) # Bisection method
		#x,r = opt.bisect(f, 0, 1, args=(phi,phib,df), full_output=True) # Bisection method
		#return opt.broyden2(f, 0, 1.0, args=(phi,phib,df)) # Bisection method
		#print r.iterations
		#return x
	else:
		#print 'newton'
		return opt.newton(f, 1.0, args=(phi,phib,df), maxiter=5000) # Newton-Raphson
		#return opt.bisect(f, 0, 1, args=(phi,phib,df)) # Bisection method
		# linestyles for different phi_big, colors for different deltaf's

#phi_s = 0.1
#phi_b = 0
#df = 0
#print 'hi', phi_b, P(phi_s, phi_b, df)
#exit()

lts = ['--','-']
#lts = ['--','-'] # minimized content for poster
colors = ['k','#2E9AFE','#04B404','r','g','c','m','y']
markers = ['x','^','d']
#markers = ['+','x','|']

pl.figure()
for i,phi_s in enumerate(phi_ss):
#	labels= r'$\phi_{s} (O) = %.2f$' % phi_s
	phibb = np.linspace(0.0, 1.0-phi_s-phi_water_min, 100)
#	pl.subplot(len(phi_bs),1,i+1)
#	for df in dfs:
	for j,df in enumerate(dfs):
		try: 
			ps = [P(phi_s,phib,df) for phib in phibb]
		except Exception as err: 
			print 'tried and failed for this df=%.3f:ps=this combo bad:phi_s=%.3f:phib=%.3f' % (df,phi_s,phib)
			print err
			continue
		#resids = [abs(f(p,phi_s,phib,df)) for p,phib in zip(ps, phibb)]

		p1 = [phib/phib for phib in phibb]#[P(phi_s,phib,df) for phib in phibb]
		resids = [abs(f(p,phi_s,phib,df)) for p,phib in zip(ps, phibb)]	
		#print resids
		#pl.plot(phibb,resids,color=colors[j],linestyle = lts[i])#color = colors[j],
		max_resid = max(resids)
		print 'Largest residual for df=%f: %e' % (df, max_resid)
		#labels= r'$\phi_{s} (O) = %.2f$' % phi_s
		if j==0:
			labels=r'$\phi_{s} = %.2f$' % phi_s
		#	marker = markers[3]
		else:
			labels=None
		if i==1:
			labelss=r'$\Delta f = %.1f$' % df
		else:
			labelss=None
		#pps = [phi_s*P(phi_s,phib,df) for phib in phibb]
		print 'Last phi tried for df=%f: phis=%.2f: phib=%.1f: total phi=%.2f' % (df, phi_s,phib,phi_s+phib)
		pl.plot(phibb,ps, color = colors[j], linewidth= 3.0,linestyle = lts[i])#,label=labels)#,color = colors[j],
		pl.plot(phibb,p1, 'k', linewidth = 0.5,linestyle = ':')#,label=labels)#,color = colors[j],
		#pl.plot(-2.0,-2.0, 'k',linewidth= 3.0,linestyle = lts[i],label=labels)#,color = colors[j],
		pl.plot(-1.0,-1.0, linestyle= '-',color = colors[j],linewidth= 3.0,label=labelss)#,color = colors[j],
	#pl.axis([0.0, 1.0, -0.2, 10])

#font = {'size': 22}
#matplotlib.rc('font',**font)
#####pl.text(0.02,12.4,r'$\Delta f = 0.0 %$',color='k')
#####pl.text(0.02,11.8,r'$\Delta f = 2.0$',color='k')
#####pl.text(0.02,11.2,r'$\Delta f = 5.0$',color='k') 
#pl.text(0.02,10.6,r'$\Delta f = 5.0$', color='g')
#pl.text(0.02,10.0,r'$\Delta f = 10.0$', color='c')
#pl.text(0.02,9.4,r'$\Delta F = 15.0$', color= 'm')

pl.legend(loc='lower right')
#pl.title('Partition Coefficient for Fixed $\Phi_{PEG400}$(Out) and Varied $\Phi_{PEG3500}$')
#pl.xlabel(r'$\phi_{s} + \phi_{b}$', size = 'x-large')
pl.ylabel(r'$p  ( \phi_{s}  , \phi_{b} , \Delta f )$', size = 'x-large')
pl.axis([0.0, 0.70, -0.1, 2.1])

pl.savefig('130128_vary_phib_PC_poster.eps', dpi=600)


pl.show()





