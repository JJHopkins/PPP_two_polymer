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
dfs = [0.0, 0.5, 2.0, 10.0]
phi_ss = [0.0, 0.01, 0.05, 0.1]

# Define equation to be solved:
# Variables are: p=partition coefficent; df=free energy cost by pore of polymer entering pore;
# phis=vol fract of short polymers<=>num fract; phib= vol(numb) fract big polymers;
# Nb,Ns= subunits in big,short polymers; Chi=1/2
def f(p, phi, phib, df):
	""" Implicit definition of P(phi,phib, dF) """
	#return -log(p) - df + (p-1)*phi + \
	#    ( phi*(1-p) + phib + \
	#      5./4*alpha*phi**(5./4)*p**(5./4)*(p*phi - (9./5)) - \
	#      phib/Nb - (5./4)*alpha*(phi + phib)**(5./4)*(phi + phib - (9./5)) - \
	#      1./2*((1 - p*phi)**2 - (1 - phi - phib)**2) ) * Ns 
	return -p + exp( - df + (p-1)*phi + \
             ( phi*(1-p) + phib + \
	      (5./4)*alpha*phi**(9./4)*p**(9./4)- \
	      (9./4)*alpha*phi**(5./4)*p**(5./4) - \
	      (1./2)*(1 - p*phi)**2 - phib/Nb - (5./4)*alpha*(phi + phib)**(9./4) + \
	      (9./4)*alpha*(phi+phib)**(5./4)+ (1./2)*(1-phi-phib)**2) * Ns) 




# Use brent optimization algorithm if f has opposite signs
# Use Newton-Raphson optimization algorith otherwise
def P(phi, phib, df):
	""" Numerically solve for partition coefficient as a
	    function of \phib """
	#print 'P'
	if f(0,phi,phib,df)*f(1,phi,phib,df) < 0:
		print 'brent'
		#return opt.brentq(f, 0, 1, args=(phi,df)) # Brent's method
		#return opt.brenth(f, 0, 1, args=(phi,df)) # Brent's method
		return opt.bisect(f, 0, 1, args=(phi,phib,df)) # Bisection method
		#x,r = opt.bisect(f, 0, 1, args=(phi,df), full_output=True) # Bisection method
		#print r.iterations
		#return x
	else:
		print 'newton'
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
print "here."
for i,phi_s in enumerate(phi_ss):
	print "i=",i
	phibb = np.linspace(0, 1.0-phi_s, 100)
#	pl.subplot(len(phi_ss),1,i+1) # use this if you want to see all phi separately
#	for df in dfs:
	for j,df in enumerate(dfs):
		print 'j=%g'%j
		try: 
			#print 'try'
			ps = [P(phi_s,phib,df) for phib in phibb]
	       	except: 
			#print 'except'
			continue
		#print 'try & exc r ok'
		resids = [abs(f(p,phi_s,phib,df)) for p,phib in zip(ps, phibb)]
		max_resid = max(resids)
		print 'Largest residual for df=%f: %e' % (df, max_resid)
		if j==3:
			labels=r'$\Phi_{PEG400} = %.2f$' % phi_s
		else:
			labels=None
		pl.plot(phibb,ps,color=colors[j],linestyle = lts[i],label=labels)
pl.axis([0, 1.3, -0.1, 25.0])

pl.text(0.04,15.0,r'$\Delta F = 0.0$', color= 'r')
pl.text(0.04,14.0,r'$\Delta F = 0.5$', color= 'g')
pl.text(0.04,13.0,r'$\Delta F = 2.0$', color= 'b')
pl.text(0.04,12.0,r'$\Delta F = 10.0$', color= 'k')

font = {'size':18}
matplotlib.rc('font',**font)

pl.legend(loc='center right')
pl.title('Partition Coefficent for Fixed  $\Phi_{PEG400}$ and Varied $\Phi_{PEG3500}$')
pl.xlabel('Number fraction of PEG3500 $\Phi_{PEG3500}$')
pl.ylabel('Partition Coefficent  P( $\Phi_{PEG400}$, $\Phi_{PEG3500}$, $\Delta$F)')

pl.savefig('Part_Coeff_Mixed_VaryBig_Corrected_DWB.png')

pl.show()


