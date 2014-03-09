#!/usr/bin/python

import numpy as np
from numpy import log,exp
import scipy.optimize as opt
import matplotlib.pyplot as pl
import matplotlib
from matplotlib.backends.backend_pdf import PdfPages
pp = PdfPages('121221_CP_ref_3_poly_mixed_in_pore_multi_same_df.pdf')

alpha = 0.49
Ns = 2.27
Nm = 22.73
Nb = 227.27
dfs = [7.0, 8.0]
phis_df_cont= [0.05,0.10]
phi_ms = [0.05, 0.10]
phi_bs = [0.05, 0.10]
phi_water_min = 0.01

#print 'CONSTANTS AND VALUES'
#print 'alpha = 0.49'
#print 'Ns = 2.27'
#print 'Nm = 22.73'
#print 'Nb = 227.27'
#print 'dfs = [2.5, 5.0, 7.5]'
#print 'phis_df_cont= [0.10,0.20]'
#print 'phi_ms = [0.15, 0.20]'
#print 'phi_bs = [0.05, 0.10]'
#print 'phi_water_min = 0.01'
#print 'dfs_cont = np.linspace(2.5,7.5 , 100)'
#print ' '

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
def Mu_s_in(phism,phimd): 
	return log(phism) + 1 - phism - (1-phism-phimd)*Ns - phimd*(Ns / Nm) +  \
	(1./2)*Ns*(1-phism-phimd)**(2) - (5./4)*alpha*Ns*(phism + phimd)**(9./4) + \
	(9./4)*alpha*Ns*(phism + phimd)**(5./4) 

def Mu_m_in(phism,phimd): 
	return log(phimd) + 1 - phimd - (1-phism-phimd)*Nm - phism*(Nm / Ns) +  \
	(1./2)*Nm*(1-phism-phimd)**(2) - (5./4)*alpha*Nm*(phism + phimd)**(9./4) + \
	(9./4)*alpha*Nm*(phism + phimd)**(5./4) 

def Mu_w(s, m, b):
	return log(1 - s - m - b) - (1 - s - m - b) + 1 - s/Ns - m/Nm - b/Nb + \
		(1./2)*(1 - (1 - s - m - b))**2 - (5./4)*alpha* \
		(1 - (1 - s - m - b))**(9./4)
lts = [':','--','-']
markers = ['x', '*']
markerss = ['x','*']
markersss = ['d','o']
colors = ['b','g','r']
markerfacecolors = ['b','g']
# Things will go wrong a lot in the code below
# Ignore most of the warnings
#np.seterr(all='ignore')

print ' '
print 'CONTINUOUS DFs and DISCRETE phi_s'
for h,phi_m in enumerate(phi_ms):
	print ' '
	print 'phi_m cont= %.2f'%(phi_m)
	for i,phi_b in enumerate(phi_bs):
		print 'phi_b cont= %.2f'%(phi_b)
		dfs_cont = np.linspace(7.0, 8.0 , 100)
		for k,phi_df_cont in enumerate(phis_df_cont):
			print 'phis_df_cont= %.2f'%(phi_df_cont)
			try: ps_df_cont = [P(phi_df_cont,phi_m,phi_b,df_cont) for df_cont in dfs_cont]
			except: continue
			phis_in_df_cont = [phi_df_cont*p_df_cont for p_df_cont in ps_df_cont]
			try: pms_df_cont=[Pm(phi_df_cont,phi_m,phi_b,phi_in_df_cont,df_cont) for phi_in_df_cont,df_cont in zip(phis_in_df_cont,dfs_cont)]
			except: continue
			phim_in_df_cont = [phi_m*pm_df_cont for pm_df_cont in pms_df_cont]
			
			pl.figure(1)
			pl.plot(dfs_cont,ps_df_cont,marker=markerss[h],color = colors[k],linestyle = lts[i],linewidth= 0.7,markersize=4,markevery=7)
			pl.title(r'Fig.1  $P_{s}$ for varied $\Delta f$, fixed $\Phi_{s}^{out}$, $\Phi_{m}^{out\ only}$, $\Phi_{b}^{out\ only}$',size = 'x-large')
			pl.xlabel(r'$\Delta f$', size = 'x-large')
			pl.ylabel(r'P$_{s}$($\Phi_{s}^{i,o}$,$\Phi_{m}^{o}$,$\Phi_{b}^{o}$,$\Delta f$)')
                        pl.savefig('1')
			
			pl.figure(2)
			pl.plot(dfs_cont,phis_in_df_cont,marker=markerss[h],color = colors[k],linestyle = lts[i],linewidth= 0.7,markersize=4,markevery=7)
			pl.title(r'Fig.2 $\Phi_{s}^{in}$ for varied $\Delta f$,fixed $\Phi_{s}^{out}$,$\Phi_{m}^{out\ only}$,$\Phi_{b}^{out\ only}$',size = 'x-large')
			pl.xlabel(r'$\Delta f$', size = 'x-large')
			pl.ylabel(r'$\Phi_{s}^{in}$($\Phi_{s}^{o}$,$\Phi_{m}^{o}$,$\Phi_{b}^{o}$,$\Delta f$)')
                        pl.savefig('2')
			
			pl.figure(3)
			pl.plot(dfs_cont,pms_df_cont,marker=markerss[h],color = colors[k],linestyle = lts[i],linewidth= 0.7,markersize=4,markevery=7)
			pl.title(r'Fig.3 $P_{m}$ for varied $\Delta f$, fixed $\Phi_{s}^{out}$, $\Phi_{m}^{in,\ out}$, $\Phi_{b}^{out\ only}$', size = 'x-large')
			pl.xlabel(r'$\Delta f$', size = 'x-large')
			pl.ylabel(r'P$_{m}$($\Phi_{s}^{i,o}$,$\Phi_{m}^{i,o}$,$\Phi_{b}^{o}$,$\Delta f$)', size = 'x-large')
                        pl.savefig('3')
			
			pl.figure(4)
			pl.plot(dfs_cont,phim_in_df_cont,marker=markerss[h],color=colors[k],linestyle=lts[i],linewidth=0.7,markersize=4,markevery=7)
			pl.title(r'Fig.4 $\Phi_{m}^{in}$ for varied $\Delta f$,fixed $\Phi_{s}^{out}$,$\Phi_{m}^{out}$, $\Phi_{b}^{out\ only}$', size = 'x-large')
			pl.xlabel(r'$\Delta f$', size = 'x-large')
			pl.ylabel(r'$\Phi_{m}^{i}$($\Phi_{s}^{i,o}$,$\Phi_{m}^{o}$,$\Phi_{b}^{o}$,$\Delta f$)', size = 'x-large')
                        pl.savefig('4')
pp.savefig(1)
pp.savefig(2)
pp.savefig(3)
pp.savefig(4)

print ' '
print 'DISCRETE DFs and CONTINUOUS phi_s'
for h,phi_m in enumerate(phi_ms):
	print ' '
	print 'phi_m discrete= %.2f'%(phi_m)
	for i,phi_b in enumerate(phi_bs):
		print 'phi_b discrete= %.2f'%(phi_b)
		for j,df in enumerate(dfs):
			phis = np.linspace(0.0, 0.10, 100)
			print 'df discrete= %.2f'%(df)
			try: ps = [P(phi,phi_m,phi_b,df) for phi in phis]
			except: continue
			resids = [abs(f(p,phi,phi_m,phi_b,df)) for p,phi in zip(ps, phis)]
			max_resid = max(resids)
			phis_ins = [phi*p for phi,p in zip(phis,ps)]
			try: pms = [Pm(phi,phi_m,phi_b,phi_ins,df) for phi,phi_ins in zip(phis,phis_ins)]
			except: continue
			phims_in = [phi_m*pm for pm in pms]
			Mus_m_in = [Mu_m_in(phi_in, phim_in) for phi_in,phim_in in zip(phis_ins,phims_in)]
			Mus_in_s_only = [Mu_s_in(phi_in, 0) for phi_in in phis_ins]
			Mus_m_s_in_diff = [Mum_in-Mus_in_sonly for Mum_in,Mus_in_sonly in zip(Mus_m_in,Mus_in_s_only)]
			NmMus_w = [Nm*Mu_w(phi,phi_m,phi_b) for phi in phis]
			NsMus_w = [Ns*Mu_w(phi,phi_m,phi_b) for phi in phis]
			Mus_m_ref = [Mu_m_in(phi_in, phim_in)-Nm*Mu_w(phi_in,phim_in,0) for phi_in,phim_in in zip(phis_ins,phims_in)]
			Mus_s_ref = [Mu_s_in(phi_in, 0)-Ns*Mu_w(phi_in,0,0) for phi_in in phis_ins]

			pl.figure(5)	
			pl.plot(phis,ps,marker = markers[h],color = colors[j],linestyle = lts[i], linewidth= 0.7, markersize= 4, markevery= 7, markerfacecolor= 'w')
			pl.title(r'Fig.5 $P_{s}$ for varied $\Phi_{s}^{out}$ fixed $\Phi_{m}^{out\ only}$ $\Phi_{b}^{out\ only}$', size = 'x-large')
			pl.xlabel(r'$\Phi_{s}^{out}$', size = 'x-large')
			pl.ylabel(r'P$_{s}$($\Phi_{s}^{i,o}$,$\Phi_{m}^{o}$,$\Phi_{b}^{o}$,$\Delta f$)', size = 'x-large')
			pl.axis([0.0, 0.101, 0.0, 0.002])
                        pl.savefig('5')
			
			pl.figure(6)	
			pl.plot(phis,phis_ins,marker=markers[h],color=colors[j],linestyle = lts[i],linewidth= 0.7, markersize= 4, markevery= 7, markerfacecolor= 'w')
			pl.title(r'Fig.6 $\Phi_{s}^{in}$ for varied $\Phi_{s}^{out}$ fixed $\Phi_{m}^{out\ only}$ $\Phi_{b}^{out\ only}$', size = 'x-large')
			pl.xlabel(r'$\Phi_{s}^{out}$', size = 'x-large')
			pl.ylabel(r'$\Phi_{s}^{in}$($\Phi_{s}^{o}$,$\Phi_{m}^{o}$,$\Phi_{b}^{o}$,$\Delta f$)')
                        pl.savefig('6')

			pl.figure(7)
			pl.plot(phis,pms,marker= markers[h],color=colors[j],linestyle = lts[i], linewidth= 0.7, markersize= 4, markevery= 7, markerfacecolor= 'w')
			pl.title(r'Fig.7 $P_{m}$ for varied $\Phi_{s}^{out}$ fixed $\Phi_{m}^{out}$ $\Phi_{b}^{out\ only}$', size = 'x-large')
			pl.xlabel(r'$\Phi_{s}^{out}$', size = 'x-large')
			pl.ylabel(r'P$_{m}$($\Phi_{s}^{i,o}$,$\Phi_{m}^{i,o}$,$\Phi_{b}^{o}$,$\Delta f$)', size = 'x-large')
                        pl.savefig('7')
		
			pl.figure(8)
			pl.plot(phis,phims_in,marker=markers[h],color = colors[j],linestyle = lts[i],linewidth= 0.7,markersize=4,markevery=7,markerfacecolor='w')
			pl.title(r'Fig.8 $\Phi_{m}^{in}$ for varied $\Phi_{s}^{out}$ fixed $\Phi_{m}^{out}$ $\Phi_{b}^{out\ only}$', size = 'x-large')
			pl.xlabel(r'$\Phi_{s}^{o}$', size = 'x-large')
			pl.ylabel(r'$\Phi_{m}^{i}$($\Phi_{s}^{i,o}$,$\Phi_{m}^{o}$,$\Phi_{b}^{o}$,$\Delta f$)', size = 'x-large')
                        pl.savefig('8')
			
			pl.figure(9)
			pl.plot(phis, Mus_in_s_only,marker=markers[h],color = colors[j],linestyle = lts[i],linewidth= 0.7,markersize=4,markevery=7,markerfacecolor='w')
			pl.title(r'Fig.9 Chemical Potential of small polymers in pore (no PEG1k in)', size = 'x-large')
			pl.ylabel(r'$\mu_{s}^{i}(\Phi_{s}^{i,o},\Phi_{m}^{o},\Phi_{b}^{o},\Delta f$)', size = 'x-large')
			pl.xlabel(r'$\Phi_{s}^{o}$', size = 'x-large')
                        pl.savefig('9')

			pl.figure(10)
			pl.plot(phis,Mus_m_in,marker=markers[h],color = colors[j],linestyle = lts[i],linewidth= 0.7,markersize=4,markevery=7,markerfacecolor='w')
			pl.title(r'Fig.10 Medium Polymer Chemical Potential in pore with small')
			pl.ylabel(r'$\mu_{m}^{i}$($\Phi_{s}^{i,o}$,$\Phi_{m}^{i,o}$,$\Phi_{b}^{o}$,$\Delta f$)', size = 'x-large')
			pl.xlabel(r'$\Phi_{s}^{o}$', size = 'x-large')
		        pl.savefig('10')
		
			pl.figure(11)
			pl.plot(phis, Mus_in_s_only,marker=markers[h],color = 'k',linestyle = lts[i],linewidth= 0.7,markersize=4,markevery=10,markerfacecolor='w')
			pl.plot(phis,Mus_m_in,marker=markers[h],color = colors[j],linestyle = lts[i],linewidth= 0.7,markersize=4,markevery=7,markerfacecolor='w')
			pl.title(r'Fig.11 Small(black), Medium(color) Chemical Potentials in Pore')
			pl.ylabel(r'$\mu_{s,m}^{i}$', size = 'x-large')
			pl.xlabel(r'$\Phi_{s}^{o}$', size = 'x-large')
                        pl.savefig('11')

			pl.figure(12)
			pl.plot(phis, NsMus_w,marker=markers[h],color = 'k',linestyle = lts[i],linewidth= 0.7,markersize=4,markevery=10,markerfacecolor='w')
			pl.plot(phis, NmMus_w,marker=markers[h],color =  colors[j],linestyle = lts[i],linewidth= 0.7,markersize=4,markevery=10,markerfacecolor='w')
			pl.title(r'Fig.12 Number of displaced water molecules* Chemical Potential of water')
			pl.ylabel(r'Ns,m$\ mu_{w}$', size = 'x-large')
			pl.xlabel(r'$\Phi_{s}^{o}$', size = 'x-large')
                        pl.savefig('12')

			pl.figure(13)
			pl.plot(phis, NsMus_w,marker=markersss[h],color = 'k',linestyle = lts[i],linewidth= 0.7,markersize=4,markevery=5,markerfacecolor='k')
			pl.plot(phis, NmMus_w,marker=markersss[h],color = colors[j],linestyle=lts[i],linewidth= 0.7,markersize=4,markevery=5,markerfacecolor=colors[j])
			pl.plot(phis,Mus_s_ref,marker=markersss[h],color='k',linestyle = lts[i],linewidth= 0.7,markersize=4,markevery=10,markerfacecolor='w')
			pl.plot(phis,Mus_m_ref,marker=markersss[h],color=colors[j],linestyle = lts[i],linewidth= 0.7,markersize=4,markevery=10,markerfacecolor='w')
			pl.title(r'Fig.13 Replacement Chemical potential')
			pl.ylabel(r'$\mu_{ref}$= ($\mu_{s,b}$-Ns,m$\ mu_{w}$)/kT', size = 'x-large')
			pl.xlabel(r'$\Phi_{s}^{o}$', size = 'x-large')
                        pl.savefig('13')
#	pl.figure()
		#	pl.semilogy(0.005,0.005, color = 'b',linewidth = 0.7, label=r'$\phi_{s} = 0.10$')#color = colors[j],
		#	pl.semilogy(0.005,0.005, color = 'g',linewidth = 0.7, label=r'$\phi_{s} = 0.20$')#color = colors[j],
		#	pl.semilogy(0.005,0.005,'k',marker = 'x',linewidth = 0.0, label=r'$\phi_{m} = 0.15$')#color = colors[j],
		#	pl.semilogy(0.005,0.005,'k',marker = '^',linewidth = 0.0, label=r'$\phi_{m} = 0.20$')#color = colors[j],
		#	pl.semilogy(0.005,0.005,'k',linestyle = ':', label=r'$\phi_{b} = 0.05$')
		#	pl.semilogy(0.005,0.005,'k',linestyle = '--',label=r'$\phi_{b} = 0.10$')
		#	pl.figure(16)
		#	pl.semilogy(phis,Mus_in_diff,marker=markers[h],color=colors[j],linestyle = lts[i],linewidth= 0.7,markersize=4,markevery=7,markerfacecolor='w')
		#	pl.title(r'chemical potential of small polymers in pore with mediums minus chemical potential of small polymers only allowed in pore$')
pl.figure(14)
pl.title(r'Fixed $\Phi_{s,m,b}^{out}$ vs. Range $\Delta f_{s}$=$\Delta f_{m}$: Legend for plots 1-4')
pl.plot(-0.005,-0.005, color = 'b',linewidth = 1.2, label=r'$\phi_{s} = 0.05$')#color = colors[j],
pl.plot(-0.005,-0.005, color = 'g',linewidth = 1.2, label=r'$\phi_{s} = 0.10$')#color = colors[j],
pl.plot(-0.005,-0.005,'k',marker = 'x',linewidth = 0.0, label=r'$\phi_{m} = 0.05$')#color = colors[j],
pl.plot(-0.005,-0.005,'k',marker = '*',linewidth = 0.0, label=r'$\phi_{m} = 0.10$')#color = colors[j],
pl.plot(-0.005,-0.005,'k',linestyle = ':', label=r'$\phi_{b} = 0.05$')
pl.plot(-0.005,-0.005,'k',linestyle = '--',label=r'$\phi_{b} = 0.10$')
pl.legend(loc='center')
pl.axis([0.0, 5.0, 0.0, 5.0])
pl.savefig('14')

pl.figure(15)
pl.title(r'Fixed df_same, phi_m,b(out) vs. range phi_s(out): Legend for plots 5-11')
#pl.title(r'Fixed $\Delta f_{same}$, $\Phi_{m,b}^{out}$ vs. Range $\Phis_{s}^{out}$: Legend for plots 7-15')
pl.plot(0.005,0.005, color = 'b',linewidth = 1.2, label=r'df = 7.0')#color = colors[j],
pl.plot(0.005,0.005, color = 'g',linewidth = 1.2, label=r'df = 8.0')#color = colors[j],
#pl.plot(0.005,0.005, color = 'r',linewidth = 1.2, label=r'df = 8.0')#color = colors[j],
pl.plot(0.005,0.005,'k',marker = 'x',linewidth = 0.0, label=r'$\phi_{m} = 0.05$')#color = colors[j],
pl.plot(0.005,0.005,'k',marker = '*',linewidth = 0.0, label=r'$\phi_{m} = 0.10$')#color = colors[j],
#pl.plot(0.005,0.005,'k',marker = 'd',linewidth = 0.0, label=r'$\phi_{m} = 0.10$')#color = colors[j],
pl.plot(0.005,0.005,'k',linestyle = ':', label=r'$\phi_{b} = 0.05$')
pl.plot(0.005,0.005,'k',linestyle = '--',label=r'$\phi_{b} = 0.10$')
pl.legend(loc='center')
pl.axis([0.0, 5.0, 0.0, 5.0])
pl.savefig('15')
#pl.legend(loc='center')

pp.savefig(5)
pp.savefig(6)
pp.savefig(7)
pp.savefig(8)
pp.savefig(9)
pp.savefig(10)
pp.savefig(11)
pp.savefig(12)
pp.savefig(13)
pp.savefig(14)
pp.savefig(15)
pp.close()
pl.show()





