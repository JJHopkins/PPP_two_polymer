#!/usr/bin/python

import numpy as np
from numpy import log,exp
import scipy.optimize as opt
import matplotlib.pyplot as pl
import matplotlib
from matplotlib.backends.backend_pdf import PdfPages
#pp = PdfPages('121220_CP_3_poly_mixed_in_pore_multi_same_df.pdf')

alpha = 0.49
Ns = 2.27
Nm = 22.73
Nb = 227.27
dfs = [7.0, 8.0]
chi = 1.88
phis_df_cont= [0.05,0.10]
phi_ms = [0.05, 0.10]
phi_bs = [0.05, 0.10]

phi_water_min = 0.01

def f(p, phi, phim, phib, df):
	""" Implicit definition of P(phi, dF) """
	return - p + exp( - df + Ns*(log((1 - p*phi)/(1 - phi - phim - phib)) + \
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
	return - pm + exp( - df + Nm*(log((1 - pm*phim)/(1 - phi - phim - phib)) + \
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


#def Mu_s_w(phi,phib): 
#	return log(phi) + 1 - Ns*(log(1-phi-phib)+1-(1-phi-phib)) + Ns*(9./4)*alpha*(phi + phib)**(5./4) + \
#		Ns*(log(1 - phi - phib) + 1 - (1 - phi - phib) - phi/Ns - phib/Nb + \
#		(1./2)*(1 - (1 - phi - phib))**2 + (5./4)*alpha* \
#		(phi + phib)**(9./4))

def Mu_s_out(phism,phimd,phibg): 
	return log(phism) + 1 - phism - (1-phism-phimd-phibg)*Ns - phimd*(Ns / Nm) - phibg*(Ns/Nb) - \
	(chi -(1./2))*Ns + (1./2)*Ns*(1-phism-phimd-phibg)**(2) - (5./4)*alpha*Ns*(phism + phimd+phibg)**(9./4) + \
	(9./4)*alpha*Ns*(phism + phimd + phibg)**(5./4) 

def Mu_RP_s_out(phism,phibg): 
	return log(phism) + 1 - phism - (1-phism-phibg)*Ns -  phibg*(Ns/Nb) - \
	(chi -(1./2))*Ns + (1./2)*Ns*(1-phism-phibg)**(2) - (5./4)*alpha*Ns*(phism +phibg)**(9./4) + \
	(9./4)*alpha*Ns*(phism + phibg)**(5./4) 


#def Mu_mat_s_out(phism,phimd,phibg): 
#	return chi*Ns*(1 + phimd + phibg) +Ns*(phism - 1)*(phism + phimd + phibg) + \
#	(9./4)*alpha*Ns*(-phism + 1)*(phism + phimd + phibg)**(5./4) - (1./2)*Ns*(phism + phimd +phibg)**2 +\
#	alpha*Ns*(phism + phimd + phibg)**(9./4) + log(phism) + 1
#
def Mu_mat_s_out(phism,phimd,phibg): 
	return chi*Ns*(phism + phimd + phibg) - (1./2)*Ns*(phism + phimd +phibg)**2 + alpha*Ns*(phism + phimd + phibg)**(9./4) + chi*Ns - \
	Ns*(phism + phimd +phibg) + (9./4)*alpha*Ns*(phism + phimd + phibg)**(5./4)  + log(phism) + 1


lts = [':','--','-']
markers = ['x', '*']
markerss = ['x','*']
colors = ['b','g','r']
markerfacecolors = ['b','g']

# Things will go wrong a lot in the code below
# Ignore most of the warnings
#np.seterr(all='ignore')

for h,phi_m in enumerate(phi_ms):
	for i,phi_b in enumerate(phi_bs):
		for j,df in enumerate(dfs):
			phis = np.linspace(0.0, 0.10, 100)
			print 'df discrete= %.2f'%(df)
			try: ps = [P(phi,phi_m,phi_b,df) for phi in phis]
			except: continue
			resids = [abs(f(p,phi,phi_m,phi_b,df)) for p,phi in zip(ps, phis)]
			max_resid = max(resids)
			#print 'Largest residual_s for df=%f: %e' % (df, max_resid)

			phis_ins = [phi*p for phi,p in zip(phis,ps)]
			try: pms = [Pm(phi,phi_m,phi_b,phi_ins,df) for phi,phi_ins in zip(phis,phis_ins)]
			except: continue
			phims_in = [phi_m*pm for pm in pms]
			
			#Hand wavey 3 poly:
			Mus_md_s_out = [Mu_s_out(phi, phi_m,phi_b) for phi in phis]
			Mus_nomd_s_out = [Mu_s_out(phi, 0,phi_b) for phi in phis]
			
			#RP's 2 poly:
			Mus_RP_nomd_s_out = [Mu_RP_s_out(phi,phi_b) for phi in phis]
			
			#Mathematica derivative for 3 poly:
			Mus_mat_md_s_out = [Mu_mat_s_out(phi, phi_m,phi_b) for phi in phis]
			Mus_mat_nomd_s_out = [Mu_mat_s_out(phi, 0,phi_b) for phi in phis]

			pl.figure()
			pl.plot(phis, Mus_mat_md_s_out,marker = 'x',color = 'b',linestyle = lts[i],linewidth= 0.7,markersize=4,markevery=10,markerfacecolor='w')
			pl.plot(phis, Mus_mat_nomd_s_out,marker='^',color = 'g',linestyle = lts[i],linewidth= 0.7,markersize=4,markevery=4,markerfacecolor='k')
			pl.plot(phis,Mus_RP_nomd_s_out,marker='d',color = 'r',linestyle = lts[i],linewidth= 0.7,markersize=4,markevery=7,markerfacecolor='m')
			#pl.plot(phis, Mus_m_s_in_diff,marker=markers[h],color=colors[j],linestyle =lts[i],linewidth= 0.7,markersize=4,markevery=7,markerfacecolor='w')
			pl.title(r'Fig.11 Small(black), Medium(color) Chemical Potentials in Pore')
			pl.ylabel(r'$\mu_{s,m}^{i}$', size = 'x-large')
			pl.xlabel(r'$\Phi_{s}^{o}$', size = 'x-large')
                        pl.savefig('mu_rp')
pl.show()






