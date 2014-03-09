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
dfs = [0.75, 1.0, 1.25]
phi_ms = [0.15, 0.20]#!
phi_bs = [0.05, 0.10]
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

def fm(pm, phi, phim, phib, phi_in, df):
	""" Implicit definition of P(phi, dF) """
	return - pm + exp( - Nm*df + Nm*(log((1 - pm*phim)/(1 - phi - phim - phib)) + \
		(pm - 1)*phim - phi - phib + (9./4)*alpha*((phi + phim + phib)**(5./4) - (phi_in + pm*phim)**(5./4))))
def Pm(phi, phim, phib, phi_in, df):
	""" Numerically solve for partition coefficient as a
	    function of \phi_s """
	if fm(0,phi,phim,phib,phi_in,df)*fm(1,phi,phim,phib,phi_in,df) < 0:
		#print 'bisect'
		return opt.bisect(fm, 0, 1, args=(phi,phim,phib,phi_in,df), maxiter=5000) # Bisection method
	else:
		#print 'newton'
		return opt.newton(fm, 10.0, args=(phi,phim,phib,phi_in,df), maxiter=5000) # Newton-Raphson


lts = [':','--','-']
markers = ['x', '^', 'd']
colors = ['b','g','r']

pl.figure()
#for h,phi_m in enumerate(phi_ms):
#	for i,phi_b in enumerate(phi_bs):
#		phis = np.linspace(0.0, 1.0 - phi_m - phi_b - phi_water_min, 100)
#		for j,df in enumerate(dfs):
#			#try: ps = [phi*P(phi,phi_m,phi_b,df) for phi in phis]
#			try: ps = [P(phi,phi_m,phi_b,df) for phi in phis]
#			except: continue
#			resids = [abs(f(p,phi,phi_m,phi_b,df)) for p,phi in zip(ps, phis)]
#			max_resid = max(resids)
#			print 'Largest residual for df=%f: %e' % (df, max_resid)
#		#	if j==1:
#		#		labels=r'$\Delta f = %.1f$' % df
#		#	else:
#		#		labels=None
#		#	if i==2:
#		#		labelss=r'$\phi_{b} = %.2f$' % phi_b
#		#	else:
#		#		labelss=None

for h,phi_m in enumerate(phi_ms):
	#print (phi_m)
	for i,phi_b in enumerate(phi_bs):
		#print (phi_b)
		phis = np.linspace(0.0, 0.7 - phi_m - phi_b - phi_water_min, 100)
		#phi_ins = np.linspace(0.0, 0.3,100)# -phi_m - phi_water_min, 100)
		for j,df in enumerate(dfs):
			#print (df)
			try: ps = [P(phi,phi_m,phi_b,df) for phi in phis]
			except: continue
			resids = [abs(f(p,phi,phi_m,phi_b,df)) for p,phi in zip(ps, phis)]
			max_resid = max(resids)
			#print 'Largest residual_s for df=%f: %e' % (df, max_resid)

			phis_in = [phi*p for phi,p in zip(phis,ps)]

			#print (phis)
			#print (phis_in)
			try: pms = [Pm(phi,phi_m,phi_b,phi_in,df) for phi,phi_in in zip(phis,phis_in)]
			
			#try: pms = [phi_m*Pm(phi,phi_m,phi_b,df) for phi in phis]
				#print 'yes'
		#	except Exception as err: 
		#		print 'tried and failed for this df=%.3f:ps=this combo bad:phi_s=%.3f:phib=%.3f' % (df,phi_b,phi_m)
		#		print err
		#	continue
			except: continue
		#	print 'no'
			#residsm = [abs(fm(pm,phi,phi_m,phi_b,phi_in,df)) for pm,phin_in,phi in zip(pms, phi_ins, phis)]
		#	pmms = [pms*phi_ms]
			#max_residm = max(residsm)
			#print 'Largest residual_m for df=%f: %e' % (df, max_residm)
		#	if j==1:
		#		labels=r'$\Delta f = %.1f$' % df
		#	else:
		#		labels=None
		#	if i==2:
		#		labelss=r'$\phi_{b} = %.2f$' % phi_b
		#	else:
		#		labelss=None
			#print (phi_m)	
			#print (pms)
			#print (phis)
			#print (phis_in)
			print 'Last phi tried_m for phi_m=%.2f: phi_b=%.2f: df=%.2f: phis=%.2f: phi_in=%.2f:' % (phi_m,phi_b,df,phi,phi_in)
			#print 'Last phi tried_m for df=%f: phis=%.2f: phib=%.1f: total phi=%.2f' % (df, phi,phi_b,phi+phi_m+phi_b)
			
			#pl.plot(phis,phis_in,marker = markers[h],color = colors[j],linestyle = lts[i], linewidth= 0.7, markersize= 4, markevery= 7, markerfacecolor= 'w')
			#pl.plot(-2.0,-2.0,'k',linestyle = lts[i], linewidth= 1.0,label=labels)#color = colors[j],
			#pl.plot(-1.0,-1.0,'k',marker= markers[j],linestyle = 'none', markersize= 6,markerfacecolor= 'w',label=labelss)#color = colors[j],

			pl.plot(phis,pms,marker= markers[h],color=colors[j],linestyle = lts[i], linewidth= 0.7, markersize= 4, markevery= 7, markerfacecolor= 'w')
			#pl.plot(-2.0,-2.0,'k',linestyle = lts[i], linewidth= 1.0,label=labelss)#color = colors[j],
			#pl.plot(-1.0,-1.0,'k',marker= markers[h],linestyle = 'none', markersize= 6,markerfacecolor= 'w',label=labelss)#color = colors[j],

pl.plot(-0.5,-0.5,'k',marker = 'x',linewidth = 0.0, label=r'$\Phi_{m} = 0.15$')#color = colors[j],
pl.plot(-0.5,-0.5,'k',marker = '^',linewidth = 0.0, label=r'$\Phi_{m} = 0.20$')#color = colors[j],

pl.plot(-0.5,-0.5,'k',linestyle = ':', label=r'$\Phi_{b} = 0.05$')
pl.plot(-0.5,-0.5,'k',linestyle = '--',label=r'$\Phi_{b} = 0.10$')

pl.plot(-0.5,-0.5,color = 'b',label=r'$\Delta f = N*0.75$')#color = colors[j],
pl.plot(-0.5,-0.5,color = 'g',label=r'$\Delta f = N*1.00$')#color = colors[j],
pl.plot(-0.5,-0.5,color = 'r',label=r'$\Delta f = N*1.26$')#color = colors[j],
pl.legend(loc='upper left')
#pl.title('Partition Coefficient for Varied $\Phi_{PEG400}$(Out) and Fixed $\Phi_{PEG3500}$(Out)')
#pl.xlabel(r'$\phi_{s} + \phi_{b}$', size = 'x-large')
#pl.ylabel(r'$p  ( \phi_{s}  , \phi_{b} , \Delta f )$', size = 'x-large')
#pl.axis([0, 0.96, -0.05, 1.0])


#pl.title('$\Phi_{PEG1k}$(in) for varied $\Phi_{PEG100}$(out only), fixed $\Phi_{PEG1k}$(out), $\Phi_{PEG10k}$(out)')
#pl.xlabel(r'$\Phi_{s}^{out}$', size = 'x-large')
#pl.ylabel(r'$\Phi_{m}^{i}$($\Phi_{s}^{o}$,$\Phi_{m}^{o}$,$\Phi_{b}^{o}$,$\Delta f$)', size = 'x-large')
#pl.axis([0, 1.0, 0.0, 1.0])
#pl.savefig('121128_PC_3_poly_phim_in_vs_phis_out_only.eps', dpi=600)


pl.title(r'$P_{PEG1k}$ for varied $\Phi_{PEG100}^{out}$ fixed $\Phi_{PEG1k}^{out}$ $\Phi_{PEG10k}^{out}$')
pl.xlabel(r'$\Phi_{s}^{out}$', size = 'x-large')
pl.ylabel(r'P$_{m}$($\Phi_{s}^{i,o}$,$\Phi_{b}^{o}$,$\Delta f$)', size = 'x-large')
#
#pl.set_xticks([0.01,0.1,0.5,1])
#pl.axis.set_xticklabels(['0.01','0.1','0.5','1'])
#pl.axis.grid()
pl.axis([0.0, 0.5, 0.0, 0.005])
#pl.minorticks_on()
pl.savefig('121129_PC_3_poly_Pm_in_with_true_phis_at_2300hr.eps', dpi=600)

#pl.title('$P_{PEG100}$ for varied $\Phi_{PEG100}^{out}$ fixed $\Phi_{PEG1k}^{out}$ $\Phi_{PEG10k}^{out}$')
#pl.xlabel(r'$\Phi_{s}^{out}$', size = 'x-large')
#pl.ylabel(r'P$_{s}$($\Phi_{m}^{o}$,$\Phi_{b}^{o}$,$\Delta f$)', size = 'x-large')
#pl.axis([0.0, 0.61, 0.0, 8.0])
#pl.savefig('121128_PC_3_poly_Ps_in_vs_phis_out.eps', dpi=600)

pl.show()





