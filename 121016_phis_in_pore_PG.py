#!/usr/bin/python
import numpy as np
from numpy import log,exp
import scipy.optimize as opt
import matplotlib.pyplot as pl
import matplotlib

#alpha = 0.49
alphas = [0.10,0.90]
Ns = 5
Nb = 77
dfs = [-6.0,-3.0, 0.0]
phi_bs = [0.0, 0.15]
phi_water_min = 0.01

def f(p, phi, phib, df):
	""" Implicit definition of P(phi, dF) """
	return - p + exp( - df + Ns*(log((1 - p*phi)/(1 - phi - phib)) + \
		(p - 1)*phi - phib + (9./4)*alpha*((phi + phib)**(5./4) - (p*phi)**(5./4))))
#Optimization algorithms, if opposite sign use brent, else use newton
def P(phi, phib, df):
	""" Numerically solve for partition coefficient as a
	    function of \phi_s """
	if f(0,phi,phib,df)*f(1,phi,phib,df) < 0:
		return opt.bisect(f, 0, 1, args=(phi,phib,df), maxiter=500) # Bisection method
	else:
		return opt.newton(f, 1.0, args=(phi,phib,df), maxiter=5000) # Newton-Raphson

#lts = [':','-.','--','-',':']
lts = ['--','-'] #['--',':','-.','-']
colors = ['b','g','k']#'r','y','k','c']#'c','m','y','k']
markers = ['*','+','^','d', '|', 'x']

pl.figure()
for k, alpha in enumerate(alphas):
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
			if j==2:
				#labels = r'$\Phi_{3400} = %.2f$' % phi_b
				labels= r'$\alpha = %.2f$' % alpha
			else:
				labels=None
			
			pps = [phi*P(phi,phi_b,df) for phi in phis]
			non_pps_in = [1.0 - phi*P(phi,phi_b,df) for phi in phis]
			non_pps_out = [1.0 - phi - phi_b for phi in phis]
			phis_shift = [0.01 + phi for phi in phis]
			print 'Last phi tried for df=%f: phis=%.2f: phib=%.1f: total phi=%.2f' % (df, phi,phi_b,phi+phi_b)
			
			#pl.plot(phis+phi_b,pps,color=colors[j],linestyle = lts[i], linewidth = 1.5, label=labels)#color = colors[j],
			#pl.plot(phis+phi_b,non_pps,color=colors[j],linestyle = lts[i], linewidth = 1.0, label=labels)#color = colors[j],
			pl.plot(phis,non_pps_in,marker = markers[k], markevery = 7, color=colors[j], linestyle = lts[i], linewidth = 1.0, label=labels)#color = colors[j],linestyle = lts[j],marker= markers[j],marker = markers[k], markevery = 10,
		#	pl.plot(phis,pps,color=colors[j], linestyle = lts[i], linewidth = 1.5, label=labels)#color = colors[j],linestyle = lts[j],marker= markers[j]
			
			#PG's data percentage and conductance values for pure PEG200
			p200 = [0.00,0.02,0.05,0.10,0.15,0.25]
			#p200 = [0.00,0.01,0.025,0.05,0.075,0.125]
			c200 = [0.939,0.864,0.751,0.590,0.392,0.295]
			p200yerr = [0.019,0.015,0.018,0.013,0.012,0.011]
			# PG's data percentage and conductance values for PEG200 (in PEG 3400)
			p200_add = [0.03,0.06,0.09,0.12,0.15]
			#p200_add = [0.015,0.03,0.045,0.06,0.075]
			c200_add = [0.792,0.666,0.527,0.431,0.324]
			p200_addyerr = [0.009,0.015,0.021,0.014,0.004]
			pl.errorbar(p200,c200, yerr = p200yerr, fmt='o', color = 'm')
			pl.errorbar(p200_add,c200_add,yerr = p200_addyerr, fmt='s', color = 'c')
			# Markers for PG's data legend
			pl.plot(0.81,0.76, marker='o', color = 'm')
			pl.plot(0.81,0.72, marker='s', color = 'c')
		pl.grid(True)


#pl.text(0.82,0.93,r'$\Delta F = -5.0$', color='b')
#pl.text(0.82,0.89,r'$\Delta F = -4.0$', color='g')
#pl.text(0.82,0.85,r'$\Delta F = -3.0$', color='r')
#pl.text(0.82,0.81,r'$\Delta F = -2.0$', color='y')
#pl.text(0.82,0.77,r'$\Delta F = \,\,0.0$', color= 'k')

pl.text(0.82,0.93,r'$\Delta F = -5.0$', color='b')
pl.text(0.82,0.89,r'$\Delta F = -3.0$', color='g')
pl.text(0.82,0.85,r'$\Delta F = \,\,0.0$', color='k')
#pl.text(0.82,0.81,r'$\Delta F = \,\, 3.0$', color='y')
#pl.text(0.82,0.77,r'$\Delta F = \,\,5.0$', color= 'r')
#legend for PG's data
pl.text(0.82,0.73,r'$\Phi_{3400} = \,\,0.00$', color= 'k')
pl.text(0.82,0.69,r'$\Phi_{3400} = \,\,0.15$', color= 'k')

pl.legend(loc='lower right')

#pl.title(r'Number Fraction of PEG200 in Pore for Varied $\Phi_{PEG200}$(Out), Fixed $\Phi_{PEG3400}$, $\Delta$F')
#pl.xlabel(r'Total Number fraction of polymer in bath $\Phi_{200}$ + $\Phi_{3400}$')
#pl.ylabel(r'Number Fraction of $\Phi_{400}$(In)')
#pl.savefig('121016_phis_in_pore_PG.eps', dpi=600)

#pl.title('Number Fraction of Water in Pore $\Phi_{water}(in)$=$\Phi_{Total}(in)$-$\Phi_{200}(in)$ For Fixed $\Phi_{3400}$')
#pl.xlabel('Total polymer number fraction in bath ($\Phi_{200}$ + $\Phi_{3400}$)')
#pl.ylabel('$\Phi_{water}$(In)')
#pl.savefig('121016_water_in_pore_fixed_phib_PG.eps', dpi=600)


#pl.title('Volume fraction of small polymer in vs. small polymer out')
#pl.xlabel('$\Phi_{200}(out)$')
#pl.ylabel('$\Phi_{200}(in)$')
#pl.savefig('121025_121016_phis_in_vs_phis_out.eps', dpi=600)

pl.title(r'Volume fraction of water in pore $\Phi_{w}(i)$ vs. $\Phi_{s}(o)$,Nb=77,Ns=5,$\alpha$=0.1,0.9')# for fixed $\Phi_{3400}$')
pl.xlabel('Small Polymer volume fraction in bath $\Phi_{s}(out)$')
pl.ylabel('$\Phi_{water}(in)$')
pl.savefig('121026_121016_phiw_in_vs_phis_out_Nb_77_Ns_5_alpha_10_90.eps', dpi=600)

pl.show()




