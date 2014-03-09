#!/usr/bin/python
import numpy as np
from numpy import log,exp
import scipy.optimize as opt
import matplotlib.pyplot as pl
import matplotlib

alpha = 0.49
Ns = 2.27
Nm = 22.73#!
Nb = 227.27
dfs = [-1.0, 0.0, 0.5]
phi_ms = [0.30, 0.40]#!
phi_bs = [0.10, 0.20]
phi_water_min = 0.01

def f(p, phi, phim, phib, df):
	""" Implicit definition of P(phi, dF) """
	return - p + exp( - Ns*df + Ns*(log((1 - p*phi)/(1 - phi - phim - phib)) + \
		(p - 1)*phi - phim - phib + (9./4)*alpha*((phi + phim + phib)**(5./4) - (p*phi)**(5./4))))
def P(phi, phim, phib, df):
	""" Numerically solve for partition coefficient as a
	    function of \phi_s """
	if f(0,phi,phim,phib,df)*f(1,phi,phim,phib,df) < 0:
		return opt.bisect(f, 0, 1, args=(phi,phim,phib,df), maxiter=500) # Bisection method
	else:
		return opt.newton(f, 1.0, args=(phi,phim,phib,df), maxiter=5000) # Newton-Raphson

def fm(pm, phi, phim, phib, df):
	""" Implicit definition of P(phi, dF) """
	return - pm + exp( - Nm*df + Nm*(log((1 - pm*phim)/(1 - phi - phim - phib)) + \
		(pm - 1)*phim - phi - phib + (9./4)*alpha*((phi + phim + phib)**(5./4) - (pm*phim)**(5./4))))
def Pm(phi, phim, phib, df):
	""" Numerically solve for partition coefficient as a
	    function of \phi_s """
	if fm(0,phi,phim,phib,df)*fm(1,phi,phim,phib,df) < 0:
		return opt.bisect(fm, 0, 1, args=(phi,phim,phib,df), maxiter=500) # Bisection method
	else:
		return opt.newton(fm, 1.0, args=(phi,phim,phib,df), maxiter=5000) # Newton-Raphson


lts = [':','--','-']
markers = ['x', '^', 'd']
colors = ['b','g','r']

pl.figure()
for h,phi_m in enumerate(phi_ms):
	for i,phi_b in enumerate(phi_bs):
		phis = np.linspace(0.0, 1.0 - phi_m - phi_b - phi_water_min, 100)
		for j,df in enumerate(dfs):
			try: ps = [phi*P(phi,phi_m,phi_b,df) for phi in phis]
			#try: ps = [P(phi,phi_m,phi_b,df) for phi in phis]
			except: continue
			resids = [abs(f(p,phi,phi_m,phi_b,df)) for p,phi in zip(ps, phis)]
			max_resid = max(resids)
			print 'Largest residual for df=%f: %e' % (df, max_resid)
		#	if j==1:
		#		labels=r'$\Delta f = %.1f$' % df
		#	else:
		#		labels=None
		#	if i==2:
		#		labelss=r'$\phi_{b} = %.2f$' % phi_b
		#	else:
		#		labelss=None

			try: pms = [Pm(phi,phi_m,phi_b,df) for phi in phis]
			#try: pms = [phi_m*Pm(phi,phi_m,phi_b,df) for phi in phis]
			except: continue
			residsm = [abs(fm(pm,phi,phi_m,phi_b,df)) for pm,phi in zip(pms, phis)]
			max_residm = max(residsm)
			print 'Largest residual_m for df=%f: %e' % (df, max_residm)
		#	if j==1:
		#		labels=r'$\Delta f = %.1f$' % df
		#	else:
		#		labels=None
		#	if i==2:
		#		labelss=r'$\phi_{b} = %.2f$' % phi_b
		#	else:
		#		labelss=None

			print 'Last phi tried for df=%f: phis=%.2f: phib=%.1f: total phi=%.2f' % (df, phi,phi_b,phi+phi_m+phi_b)
			print 'Last phi tried_m for df=%f: phis=%.2f: phib=%.1f: total phi=%.2f' % (df, phi,phi_b,phi+phi_m+phi_b)
			
			#pl.plot(phis,ps,marker = markers[h],color = colors[j],linestyle = lts[i], linewidth= 0.5, markersize= 4, markevery= 7, markerfacecolor= 'w')
			#pl.plot(-2.0,-2.0,'k',linestyle = lts[i], linewidth= 1.0,label=labels)#color = colors[j],
			#pl.plot(-1.0,-1.0,'k',marker= markers[j],linestyle = 'none', markersize= 6,markerfacecolor= 'w',label=labelss)#color = colors[j],


			#pmms = [phi_m*pm for pm in pms]

			pl.plot(phis,pms,marker= markers[h],color=colors[j],linestyle = lts[i], linewidth= 1.0, markersize= 4, markevery= 7, markerfacecolor= 'w')
			#pl.plot(-2.0,-2.0,'k',linestyle = lts[i], linewidth= 1.0,label=labelss)#color = colors[j],
			#pl.plot(-1.0,-1.0,'k',marker= markers[h],linestyle = 'none', markersize= 6,markerfacecolor= 'w',label=labelss)#color = colors[j],

pl.plot(-2.0,-2.0,'k',marker = 'x',linewidth = 0.0, label=r'$\Phi_{m} = 0.30$')#color = colors[j],
pl.plot(-2.0,-2.0,'k',marker = '^',linewidth = 0.0, label=r'$\Phi_{m} = 0.40$')#color = colors[j],

pl.plot(-2.0,-2.0,'k',linestyle = ':', label=r'$\Phi_{b} = 0.10$')
pl.plot(-2.0,-2.0,'k',linestyle = '--',label=r'$\Phi_{b} = 0.20$')

pl.plot(-2.0,-2.0,color = 'b',label=r'$\Delta f = -1.0$')#color = colors[j],
pl.plot(-2.0,-2.0,color = 'g',label=r'$\Delta f = 0.0$')#color = colors[j],
pl.plot(-2.0,-2.0,color = 'r',label=r'$\Delta f = 0.5$')#color = colors[j],
pl.legend(loc='lower right')

#pl.title('Partition Coefficient for Varied $\Phi_{PEG400}$(Out) and Fixed $\Phi_{PEG3500}$(Out)')
#pl.xlabel(r'$\phi_{s} + \phi_{b}$', size = 'x-large')
#pl.ylabel(r'$p  ( \phi_{s}  , \phi_{b} , \Delta f )$', size = 'x-large')
#pl.axis([0, 0.96, -0.05, 1.0])


#number fraction for smalls
#pl.title('$\Phi_{PEG100}$(in) for varied $\Phi_{PEG100}$, fixed $\Phi_{PEG1k}$(out), $\Phi_{PEG10k}$(out)')
#pl.xlabel(r'$\Phi_{s}^{out}$', size = 'x-large')
#pl.ylabel(r'$\Phi_{s}^{i}$($\Phi_{s}^{o}$,$\Phi_{m}^{o}$,$\Phi_{b}^{o}$,$\Delta f$)', size = 'x-large')
#pl.axis([0, 0.6, 0.0, 1.0])
#pl.savefig('121129_PC_3_poly_NFs_phis.eps', dpi=600)

#number fraction for medium
pl.title('$\Phi_{PEG1k}$(in) for varied $\Phi_{PEG100}$(out only), fixed $\Phi_{PEG1k}$(out), $\Phi_{PEG10k}$(out)')
pl.xlabel(r'$\Phi_{s}^{out}$', size = 'x-large')
pl.ylabel(r'$\Phi_{m}^{i}$($\Phi_{s}^{o}$,$\Phi_{m}^{o}$,$\Phi_{b}^{o}$,$\Delta f$)', size = 'x-large')
pl.axis([0, 1.0, 0.0, 1.0])
pl.savefig('121128_PC_3_poly_NFm_phim.eps', dpi=600)

#partition coeff for mediums
#pl.title(r'$P_{PEG1k}$ for varied $\Phi_{PEG100}^{out\ only}$ fixed $\Phi_{PEG1k}^{out}$ $\Phi_{PEG10k}^{out}$')
#pl.xlabel(r'$\Phi_{s}^{out}$', size = 'x-large')
#pl.ylabel(r'P$_{m}$($\Phi_{s}^{o}$,$\Phi_{b}^{o}$,$\Delta f$)', size = 'x-large')
#pl.axis([0.0, 0.61, 0.0, 8.0])
#pl.savefig('121129_PC_3_poly_Pm_phis_out_only_zoom.eps', dpi=600)

#partition coeff for smalls
#pl.title('$P_{PEG100}$ for varied $\Phi_{PEG100}^{out}$ fixed $\Phi_{PEG1k}^{out}$ $\Phi_{PEG10k}^{out}$')
#pl.xlabel(r'$\Phi_{s}^{out}$', size = 'x-large')
#pl.ylabel(r'P$_{s}$($\Phi_{m}^{o}$,$\Phi_{b}^{o}$,$\Delta f$)', size = 'x-large')
#pl.axis([0.0, 0.61, 0.0, 8.0])
#pl.savefig('121128_PC_3_poly_Ps_in_vs_phis_out.eps', dpi=600)

pl.show()





