#!/usr/bin/python

import numpy as np
from numpy import log,exp
import scipy.optimize as opt
import matplotlib.pyplot as pl
import matplotlib
from matplotlib.backends.backend_pdf import PdfPages
pp = PdfPages('121219_CP_3_poly_mixed_in_pore_multi.pdf')

alpha = 0.49
Ns = 2.27
Nm = 22.73
Nb = 227.27
dfs = [0.75, 1.0, 1.25]
phis_df_cont= [0.10,0.20]
phi_ms = [0.15, 0.20]
phi_bs = [0.05, 0.10]
phi_water_min = 0.01

print 'CONSTANTS AND VALUES'
print 'alpha = 0.49'
print 'Ns = 2.27'
print 'Nm = 22.73'
print 'Nb = 227.27'
print 'dfs = [0.75, 1.0, 1.25]'
print 'phis_df_cont= [0.10,0.20]'
print 'phi_ms = [0.15, 0.20]'
print 'phi_bs = [0.05, 0.10]'
print 'phi_water_min = 0.01'
print 'dfs_cont = np.linspace(0.75, 1.25 , 100)'
print ' '

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


#def Mu_s_w(phi,phib): 
#	return log(phi) + 1 - Ns*(log(1-phi-phib)+1-(1-phi-phib)) + Ns*(9./4)*alpha*(phi + phib)**(5./4) + \
#		Ns*(log(1 - phi - phib) + 1 - (1 - phi - phib) - phi/Ns - phib/Nb + \
#		(1./2)*(1 - (1 - phi - phib))**2 + (5./4)*alpha* \
#		(phi + phib)**(9./4))

def Mu_s_in(phism,phimd): 
	return log(phism) + 1 - phism - (1-phism-phimd)*Ns - phimd*(Ns / Nm) +  \
	(1./2)*Ns*(1-phism-phimd)**(2) - (5./4)*alpha*Ns*(phism + phimd)**(9./4) + \
	(9./4)*alpha*Ns*(phism + phimd)**(5./4) 

def Mu_m_in(phism,phimd): 
	return log(phimd) + 1 - phimd - (1-phism-phimd)*Nm - phism*(Nm / Ns) +  \
	(1./2)*Nm*(1-phism-phimd)**(2) - (5./4)*alpha*Nm*(phism + phimd)**(9./4) + \
	(9./4)*alpha*Nm*(phism + phimd)**(5./4) 

lts = [':','--','-']
markers = ['x', '*', 'd']
markerss = ['x','*']
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
		dfs_cont = np.linspace(0.75, 1.25 , 100)
		#phis = np.linspace(0.0, 0.7 - phi_m - phi_b - phi_water_min, 100)
		for k,phi_df_cont in enumerate(phis_df_cont):
			#phi_ins = np.linspace(0.0, 0.3,100)# -phi_m - phi_water_min, 100)
			#dfs_cont = np.linspace(0.75, 1.25 , 100)
			print 'phis_df_cont= %.2f'%(phi_df_cont)
			#print (dfs_cont)
			try: ps_df_cont = [P(phi_df_cont,phi_m,phi_b,df_cont) for df_cont in dfs_cont]
			#print 'phi_s for cont df for ps is %.3f' % (phis_s_df_cont)
			except: continue
			#print 'phi_s for cont df for ps is %.2f' % (phi_s_df_cont)
			phis_in_df_cont = [phi_df_cont*p_df_cont for p_df_cont in ps_df_cont]
			try: pms_df_cont=[Pm(phi_df_cont,phi_m,phi_b,phi_in_df_cont,df_cont) for phi_in_df_cont,df_cont in zip(phis_in_df_cont,dfs_cont)]
			#print 'phi_s for cont df for pm is %f' % (phi_s_df_cont)
			#print 'pm for cont df is %' % (pms_df_cont)
			except: continue
			phim_in_df_cont = [phi_m*pm_df_cont for pm_df_cont in pms_df_cont]
			
#			pl.figure()
#			pl.plot(-0.005,-0.005, color = 'b',linewidth = 0.7, label=r'$\phi_{s} = 0.10$')#color = colors[j],
#			pl.plot(-0.005,-0.005, color = 'g',linewidth = 0.7, label=r'$\phi_{s} = 0.20$')#color = colors[j],
#			pl.plot(-0.005,-0.005,'k',marker = 'x',linewidth = 0.0, label=r'$\phi_{m} = 0.15$')#color = colors[j],
#			pl.plot(-0.005,-0.005,'k',marker = '^',linewidth = 0.0, label=r'$\phi_{m} = 0.20$')#color = colors[j],
#			pl.plot(-0.005,-0.005,'k',linestyle = ':', label=r'$\phi_{b} = 0.05$')
#			pl.plot(-0.005,-0.005,'k',linestyle = '--',label=r'$\phi_{b} = 0.10$')
			
			pl.figure(1)
			pl.plot(dfs_cont,ps_df_cont,marker=markerss[h],color = colors[k],linestyle = lts[i],linewidth= 0.7,markersize=4,markevery=7)
			pl.title(r'FOR DIFFERENT dfs:Fig.1$P_{s}$ for varied $\Delta f$, fixed $\Phi_{s}^{o}$, $\Phi_{m}^{o\ only}$, $\Phi_{b}^{o\ only}$')
			pl.xlabel(r'$\Delta f$', size = 'x-large')
			pl.ylabel(r'P$_{s}$($\Phi_{s}^{i,o}$,$\Phi_{m}^{o}$,$\Phi_{b}^{o}$,$\Delta f$)', size = 'x-large')
			pl.axis([0.7, 1.3, 0.1, 0.7])
			pl.legend(loc='upper right')
			
#			pl.figure()
#			pl.plot(-0.005,-0.005, color = 'b',linewidth = 0.7, label=r'$\phi_{s} = 0.10$')#color = colors[j],
#			pl.plot(-0.005,-0.005, color = 'g',linewidth = 0.7, label=r'$\phi_{s} = 0.20$')#color = colors[j],
#			pl.plot(-0.005,-0.005,'k',marker = 'x',linewidth = 0.0, label=r'$\phi_{m} = 0.15$')#color = colors[j],
#			pl.plot(-0.005,-0.005,'k',marker = '^',linewidth = 0.0, label=r'$\phi_{m} = 0.20$')#color = colors[j],
#			pl.plot(-0.005,-0.005,'k',linestyle = ':', label=r'$\phi_{b} = 0.05$')
#			pl.plot(-0.005,-0.005,'k',linestyle = '--',label=r'$\phi_{b} = 0.10$')
			
			pl.figure(2)
			pl.plot(dfs_cont,phis_in_df_cont,marker=markerss[h],color = colors[k],linestyle = lts[i],linewidth= 0.7,markersize=4,markevery=7)
			pl.title(r'Fig.2 $\Phi_{PEG100}^{in}$ for varied $\Delta f$,fixed $\Phi_{PEG100}^{out}$,$\Phi_{PEG1k}^{out\ only}$,$\Phi_{PEG10k}^{out\ only}$')
			pl.xlabel(r'$\Delta f$', size = 'x-large')
			pl.ylabel(r'$\Phi_{s}^{in}$($\Phi_{s}^{o}$,$\Phi_{m}^{o}$,$\Phi_{b}^{o}$,$\Delta f$)', size = 'x-large')
			pl.axis([0.745, 1.255, 0.01, 0.13])
			pl.legend(loc='upper right')
			
#			pl.figure()
#			pl.plot(-0.005,-0.005, color = 'b',linewidth = 0.7, label=r'$\phi_{s} = 0.10$')#color = colors[j],
#			pl.plot(-0.005,-0.005, color = 'g',linewidth = 0.7, label=r'$\phi_{s} = 0.20$')#color = colors[j],
#			pl.plot(-0.005,-0.005,'k',marker = 'x',linewidth = 0.0, label=r'$\phi_{m} = 0.15$')#color = colors[j],
#			pl.plot(-0.005,-0.005,'k',marker = '^',linewidth = 0.0, label=r'$\phi_{m} = 0.20$')#color = colors[j],
#			pl.plot(-0.005,-0.005,'k',linestyle = ':', label=r'$\phi_{b} = 0.05$')
#			pl.plot(-0.005,-0.005,'k',linestyle = '--',label=r'$\phi_{b} = 0.10$')
			
			pl.figure(3)
			pl.plot(dfs_cont,pms_df_cont,marker=markers[h],color = colors[k],linestyle = lts[i],linewidth= 0.7,markersize=4,markevery=7)
			pl.title(r'Fig.3 $P_{PEG1k}$ for varied $\Delta f$, fixed $\Phi_{PEG100}^{out}$, $\Phi_{PEG1k}^{in,\ out}$, $\Phi_{PEG10k}^{out\ only}$')
			pl.xlabel(r'$\Delta f$', size = 'x-large')
			pl.ylabel(r'P$_{m}$($\Phi_{s}^{i,o}$,$\Phi_{m}^{i,o}$,$\Phi_{b}^{o}$,$\Delta f$)', size = 'x-large')
			pl.axis([0.745, 1.255, 0.0, 0.016])
			pl.legend(loc='upper right')

			#pl.figure()
			#pl.semilogy(0.005,0.005, color = 'b',linewidth = 0.7, label=r'$\phi_{s} = 0.10$')#color = colors[j],
			#pl.semilogy(0.005,0.005, color = 'g',linewidth = 0.7, label=r'$\phi_{s} = 0.20$')#color = colors[j],
			#pl.semilogy(0.005,0.005,'k',marker = 'x',linewidth = 0.0, label=r'$\phi_{m} = 0.15$')#color = colors[j],
			#pl.semilogy(0.005,0.005,'k',marker = '^',linewidth = 0.0, label=r'$\phi_{m} = 0.20$')#color = colors[j],
			#pl.semilogy(0.005,0.005,'k',linestyle = ':', label=r'$\phi_{b} = 0.05$')
#			pl.semilogy(0.005,0.005,'k',linestyle = '--',label=r'$\phi_{b} = 0.10$')
			
			pl.figure(4)
			pl.semilogy(dfs_cont,pms_df_cont,marker=markers[h],color=colors[k],linestyle=lts[i],linewidth=0.7,markersize=4,markevery=7)
			pl.title(r'Fig.4 Zoom Log Fig.3 $P_{PEG1k}$ for varied $\Delta f$,fixed $\Phi_{s}^{out}$,$\Phi_{m}^{in,\ out}$,$\Phi_{b}^{out\ only}$')
			pl.xlabel(r'$\Delta f$', size = 'x-large')
			pl.ylabel(r'P$_{m}$($\Phi_{s}^{i,o}$,$\Phi_{m}^{i,o}$,$\Phi_{b}^{o}$,$\Delta f$)', size = 'x-large')
			pl.axis([0.745, 1.255, 0.0, 0.02])
			pl.legend(loc='upper right')

#			pl.figure()
#			pl.plot(-0.005,-0.005, color = 'b',linewidth = 0.7, label=r'$\phi_{s} = 0.10$')#color = colors[j],
#			pl.plot(-0.005,-0.005, color = 'g',linewidth = 0.7, label=r'$\phi_{s} = 0.20$')#color = colors[j],
#			pl.plot(-0.005,-0.005,'k',marker = 'x',linewidth = 0.0, label=r'$\phi_{m} = 0.15$')#color = colors[j],
#			pl.plot(-0.005,-0.005,'k',marker = '^',linewidth = 0.0, label=r'$\phi_{m} = 0.20$')#color = colors[j],
#			pl.plot(-0.005,-0.005,'k',linestyle = ':', label=r'$\phi_{b} = 0.05$')
#			pl.plot(-0.005,-0.005,'k',linestyle = '--',label=r'$\phi_{b} = 0.10$')
			
			pl.figure(5)
			pl.plot(dfs_cont,phim_in_df_cont,marker=markers[h],color=colors[k],linestyle=lts[i],linewidth=0.7,markersize=4,markevery=7)
			pl.title(r'Fig.5 $\Phi_{PEG1k}^{in}$ for varied $\Delta f$,fixed $\Phi_{PEG100}^{out}$,$\Phi_{PEG1k}^{out}$, $\Phi_{PEG10k}^{out\ only}$')
			pl.xlabel(r'$\Delta f$', size = 'x-large')
			pl.ylabel(r'$\Phi_{m}^{i}$($\Phi_{s}^{i,o}$,$\Phi_{m}^{o}$,$\Phi_{b}^{o}$,$\Delta f$)', size = 'x-large')
			pl.axis([0.745, 1.255, 0.0, 0.0032])
			pl.legend(loc='upper right')

			#pl.figure()
			#pl.semilogy(0.005,0.005, color = 'b',linewidth = 0.7, label=r'$\phi_{s} = 0.10$')#color = colors[j],
			#pl.semilogy(0.005,0.005, color = 'g',linewidth = 0.7, label=r'$\phi_{s} = 0.20$')#color = colors[j],
			#pl.semilogy(0.005,0.005,'k',marker = 'x',linewidth = 0.0, label=r'$\phi_{m} = 0.15$')#color = colors[j],
			#pl.semilogy(0.005,0.005,'k',marker = '^',linewidth = 0.0, label=r'$\phi_{m} = 0.20$')#color = colors[j],
			#pl.semilogy(0.005,0.005,'k',linestyle = ':', label=r'$\phi_{b} = 0.05$')
			#pl.semilogy(0.005,0.005,'k',linestyle = '--',label=r'$\phi_{b} = 0.10$')
			
			pl.figure(6)
			pl.semilogy(dfs_cont,phim_in_df_cont,marker=markers[h],color=colors[k],linestyle =lts[i],linewidth=0.7,markersize=4,markevery=7)
			pl.title(r'Fig.6 Zoom Log Fig.5 $\Phi_{PEG1k}^{in}$ for varied $\Delta f$, fixed $\Phi_{s}^{out}$, $\Phi_{m}^{out}$ $\Phi_{b}^{out\ only}$')
			pl.xlabel(r'$\Delta f$', size = 'x-large')
			pl.ylabel(r'$\Phi_{m}^{i}$($\Phi_{s}^{i,o}$,$\Phi_{m}^{o}$,$\Phi_{b}^{o}$,$\Delta f$)', size = 'x-large')
			pl.axis([0.745, 1.255, 0.0, 0.0035])
			pl.legend(loc='upper right')
pp.savefig(1)
pp.savefig(2)
pp.savefig(3)
pp.savefig(4)
pp.savefig(5)
pp.savefig(6)

print ' '
print 'DISCRETE DFs and CONTINUOUS phi_s'
for h,phi_m in enumerate(phi_ms):
	print ' '
	print 'phi_m discrete= %.2f'%(phi_m)
	for i,phi_b in enumerate(phi_bs):
		print 'phi_b discrete= %.2f'%(phi_b)
		for j,df in enumerate(dfs):
			phis = np.linspace(0.0, 0.7 - phi_m - phi_b - phi_water_min, 100)
			print 'df discrete= %.2f'%(df)
			try: ps = [P(phi,phi_m,phi_b,df) for phi in phis]
			except: continue
			resids = [abs(f(p,phi,phi_m,phi_b,df)) for p,phi in zip(ps, phis)]
			max_resid = max(resids)
			#print 'Largest residual_s for df=%f: %e' % (df, max_resid)

			phis_ins = [phi*p for phi,p in zip(phis,ps)]

			#print (phis)
			#print (phis_in)
			try: pms = [Pm(phi,phi_m,phi_b,phi_ins,df) for phi,phi_ins in zip(phis,phis_ins)]
			
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
			#print 'Last phi tried_m for phi_m=%.2f: phi_b=%.2f: df=%.2f: phis=%.2f: phi_in=%.2f:' % (phi_m,phi_b,df,phi,phi_in)
			#print 'Last phi tried_m for df=%f: phis=%.2f: phib=%.1f: total phi=%.2f' % (df, phi,phi_b,phi+phi_m+phi_b)


			phims_in = [phi_m*pm for pm in pms]
			
			#Mus_s_in = [Mu_s_in(phis*p, phi_m*pm)for pm in pms]
#			Mus_s_in = [Mu_s_in(phi_in, phim_in) for phi_in,phim_in in zip(phis_ins,phims_in)]
			Mus_m_in = [Mu_m_in(phi_in, phim_in) for phi_in,phim_in in zip(phis_ins,phims_in)]
			Mus_in_s_only = [Mu_s_in(phi_in, 0) for phi_in in phis_ins]
#			Mus_in_diff = [Mu_s_in-Mu_in_s_only for Mu_s_in,Mu_in_s_only in zip(Mus_s_in,Mus_in_s_only)]
			Mus_m_s_in_diff = [Mum_in-Mus_in_sonly for Mum_in,Mus_in_sonly in zip(Mus_m_in,Mus_in_s_only)]

	#		pl.figure()
	#		pl.plot(-0.005,-0.005,'k',marker = 'x',linewidth = 0.0, label=r'$\Phi_{s} = 0.10$')#color = colors[j],
	#		pl.plot(-0.005,-0.005,'k',marker = '^',linewidth = 0.0, label=r'$\Phi_{s} = 0.20$')#color = colors[j],
	#		pl.plot(-0.005,-0.005,'k',marker = 'x',linewidth = 0.0, label=r'$\Phi_{m} = 0.15$')#color = colors[j],
	#		pl.plot(-0.005,-0.005,'k',marker = '^',linewidth = 0.0, label=r'$\Phi_{m} = 0.20$')#color = colors[j],
	#		pl.plot(-0.005,-0.005,'k',linestyle = ':', label=r'$\Phi_{b} = 0.05$')
	#		pl.plot(-0.005,-0.005,'k',linestyle = '--',label=r'$\Phi_{b} = 0.10$')
	#		pl.plot(-0.005,-0.005,color = 'b',label=r'$\Delta f = N*0.75$')#color = colors[j],
	#		pl.plot(-0.005,-0.005,color = 'g',label=r'$\Delta f = N*1.00$')#color = colors[j],
	#		pl.plot(-0.005,-0.005,color = 'r',label=r'$\Delta f = N*1.25$')#color = colors[j],
			pl.figure(7)	
			pl.plot(phis,ps,marker = markers[h],color = colors[j],linestyle = lts[i], linewidth= 0.7, markersize= 4, markevery= 7, markerfacecolor= 'w')
			pl.title(r'$P_{PEG100}$ for varied $\Phi_{PEG100}^{out}$ fixed $\Phi_{PEG1k}^{out\ only}$ $\Phi_{PEG10k}^{out\ only}$')
			pl.xlabel(r'$\Phi_{s}^{out}$', size = 'x-large')
			pl.ylabel(r'P$_{s}$($\Phi_{s}^{i,o}$,$\Phi_{m}^{o}$,$\Phi_{b}^{o}$,$\Delta f$)', size = 'x-large')
			#pl.axis([-0.05, 0.5, 0.0, 15])
#			pl.legend(loc='upper right')

	#		pl.figure()
	#		pl.plot(-0.005,-0.005,'k',marker = 'x',linewidth = 0.0, label=r'$\Phi_{s} = 0.10$')#color = colors[j],
	#		pl.plot(-0.005,-0.005,'k',marker = '^',linewidth = 0.0, label=r'$\Phi_{s} = 0.20$')#color = colors[j],
	#		pl.plot(-0.005,-0.005,'k',marker = 'x',linewidth = 0.0, label=r'$\Phi_{m} = 0.15$')#color = colors[j],
	#		pl.plot(-0.005,-0.005,'k',marker = '^',linewidth = 0.0, label=r'$\Phi_{m} = 0.20$')#color = colors[j],
	#		pl.plot(-0.005,-0.005,'k',linestyle = ':', label=r'$\Phi_{b} = 0.05$')
	#		pl.plot(-0.005,-0.005,'k',linestyle = '--',label=r'$\Phi_{b} = 0.10$')
	#		pl.plot(-0.005,-0.005,color = 'b',label=r'$\Delta f = N*0.75$')#color = colors[j],
	#		pl.plot(-0.005,-0.005,color = 'g',label=r'$\Delta f = N*1.00$')#color = colors[j],
	#		pl.plot(-0.005,-0.005,color = 'r',label=r'$\Delta f = N*1.25$')#color = colors[j],
			pl.figure(8)	
			pl.plot(phis,phis_ins,marker=markers[h],color=colors[j],linestyle = lts[i],linewidth= 0.7, markersize= 4, markevery= 7, markerfacecolor= 'w')
			pl.title(r'$\Phi_{PEG100}^{in}$ for varied $\Phi_{PEG100}^{out}$ fixed $\Phi_{PEG1k}^{out\ only}$ $\Phi_{PEG10k}^{out\ only}$')
			pl.xlabel(r'$\Phi_{s}^{out}$', size = 'x-large')
			pl.ylabel(r'$\Phi_{PEG100}^{in}$($\Phi_{s}^{o}$,$\Phi_{m}^{o}$,$\Phi_{b}^{o}$,$\Delta f$)', size = 'x-large')
			pl.axis([0, 0.5, 0.0, 0.4])
#			pl.legend(loc='upper right')

	#		pl.figure()
	#		pl.plot(-0.005,-0.005,'k',marker = 'x',linewidth = 0.0, label=r'$\Phi_{s} = 0.10$')#color = colors[j],
	#		pl.plot(-0.005,-0.005,'k',marker = '^',linewidth = 0.0, label=r'$\Phi_{s} = 0.20$')#color = colors[j],
	#		pl.plot(-0.005,-0.005,'k',marker = 'x',linewidth = 0.0, label=r'$\Phi_{m} = 0.15$')#color = colors[j],
	#		pl.plot(-0.005,-0.005,'k',marker = '^',linewidth = 0.0, label=r'$\Phi_{m} = 0.20$')#color = colors[j],
	#		pl.plot(-0.005,-0.005,'k',linestyle = ':', label=r'$\Phi_{b} = 0.05$')
	#		pl.plot(-0.005,-0.005,'k',linestyle = '--',label=r'$\Phi_{b} = 0.10$')
	#		pl.plot(-0.005,-0.005,color = 'b',label=r'$\Delta f = N*0.75$')#color = colors[j],
	#		pl.plot(-0.005,-0.005,color = 'g',label=r'$\Delta f = N*1.00$')#color = colors[j],
	#		pl.plot(-0.005,-0.005,color = 'r',label=r'$\Delta f = N*1.25$')#color = colors[j],
			pl.figure(9)
			pl.plot(phis,pms,marker= markers[h],color=colors[j],linestyle = lts[i], linewidth= 0.7, markersize= 4, markevery= 7, markerfacecolor= 'w')
			pl.title(r'$P_{PEG1k}$ for varied $\Phi_{PEG100}^{out}$ fixed $\Phi_{PEG1k}^{out}$ $\Phi_{PEG10k}^{out\ only}$')
			pl.xlabel(r'$\Phi_{s}^{out}$', size = 'x-large')
			pl.ylabel(r'P$_{m}$($\Phi_{s}^{i,o}$,$\Phi_{m}^{i,o}$,$\Phi_{b}^{o}$,$\Delta f$)', size = 'x-large')
			pl.axis([0, 0.6, 0.0, 0.6])
#			pl.legend(loc='upper right')

		#	pl.figure()
		#	pl.semilogy(0.005,0.005, color = 'b',linewidth = 0.7, label=r'$\phi_{s} = 0.10$')#color = colors[j],
		#	pl.semilogy(0.005,0.005, color = 'g',linewidth = 0.7, label=r'$\phi_{s} = 0.20$')#color = colors[j],
		#	pl.semilogy(0.005,0.005,'k',marker = 'x',linewidth = 0.0, label=r'$\phi_{m} = 0.15$')#color = colors[j],
		#	pl.semilogy(0.005,0.005,'k',marker = '^',linewidth = 0.0, label=r'$\phi_{m} = 0.20$')#color = colors[j],
		#	pl.semilogy(0.005,0.005,'k',linestyle = ':', label=r'$\phi_{b} = 0.05$')
		#	pl.semilogy(0.005,0.005,'k',linestyle = '--',label=r'$\phi_{b} = 0.10$')
			pl.figure(10)
			pl.plot(phis,pms,marker= markers[h],color=colors[j],linestyle = lts[i], linewidth= 0.7, markersize= 4, markevery= 7, markerfacecolor= 'w')
			pl.title(r'Zoom $P_{PEG1k}$ for varied $\Phi_{PEG100}^{out}$ fixed $\Phi_{PEG1k}^{out}$ $\Phi_{PEG10k}^{out\ only}$')
			pl.xlabel(r'$\Phi_{s}^{out}$', size = 'x-large')
			pl.ylabel(r'P$_{m}$($\Phi_{s}^{i,o}$,$\Phi_{m}^{i,o}$,$\Phi_{b}^{o}$,$\Delta f$)', size = 'x-large')
			pl.axis([0.0, 0.4, 0.0, 0.0004])
#			pl.legend(loc='upper right')

		#	pl.figure()
		#	pl.plot(-0.005,-0.005,'k',marker = 'x',linewidth = 0.0, label=r'$\Phi_{s} = 0.10$')#color = colors[j],
		#	pl.plot(-0.005,-0.005,'k',marker = '^',linewidth = 0.0, label=r'$\Phi_{s} = 0.20$')#color = colors[j],
		#	pl.plot(-0.005,-0.005,'k',marker = 'x',linewidth = 0.0, label=r'$\Phi_{m} = 0.15$')#color = colors[j],
		#	pl.plot(-0.005,-0.005,'k',marker = '^',linewidth = 0.0, label=r'$\Phi_{m} = 0.20$')#color = colors[j],
		#	pl.plot(-0.005,-0.005,'k',linestyle = ':', label=r'$\Phi_{b} = 0.05$')
		#	pl.plot(-0.005,-0.005,'k',linestyle = '--',label=r'$\Phi_{b} = 0.10$')
		#	pl.plot(-0.005,-0.005,color = 'b',label=r'$\Delta f = N*0.75$')#color = colors[j],
		#	pl.plot(-0.005,-0.005,color = 'g',label=r'$\Delta f = N*1.00$')#color = colors[j],
		#	pl.plot(-0.005,-0.005,color = 'r',label=r'$\Delta f = N*1.25$')#color = colors[j],
			pl.figure(11)
			pl.plot(phis,phims_in,marker=markers[h],color = colors[j],linestyle = lts[i],linewidth= 0.7,markersize=4,markevery=7,markerfacecolor='w')
			pl.title(r'$\Phi_{PEG1k}^{in}$ for varied $\Phi_{PEG100}^{out}$ fixed $\Phi_{PEG1k}^{out}$ $\Phi_{PEG10k}^{out\ only}$')
			pl.xlabel(r'$\Delta f$', size = 'x-large')
			pl.ylabel(r'$\Phi_{m}^{i}$($\Phi_{s}^{i,o}$,$\Phi_{m}^{o}$,$\Phi_{b}^{o}$,$\Delta f$)', size = 'x-large')
			pl.axis([0, 0.5, 0.0, 0.1])
#			pl.legend(loc='upper right')
			
		#	pl.figure()
		#	pl.semilogy(0.005,0.005, color = 'b',linewidth = 0.7, label=r'$\phi_{s} = 0.10$')#color = colors[j],
		#	pl.semilogy(0.005,0.005, color = 'g',linewidth = 0.7, label=r'$\phi_{s} = 0.20$')#color = colors[j],
		#	pl.semilogy(0.005,0.005,'k',marker = 'x',linewidth = 0.0, label=r'$\phi_{m} = 0.15$')#color = colors[j],
		#	pl.semilogy(0.005,0.005,'k',marker = '^',linewidth = 0.0, label=r'$\phi_{m} = 0.20$')#color = colors[j],
		#	pl.semilogy(0.005,0.005,'k',linestyle = ':', label=r'$\phi_{b} = 0.05$')
		#	pl.semilogy(0.005,0.005,'k',linestyle = '--',label=r'$\phi_{b} = 0.10$')
			pl.figure(12)
			pl.plot(phis,phims_in,marker=markers[h],color = colors[j],linestyle = lts[i],linewidth= 0.7,markersize=4,markevery=7,markerfacecolor='w')
			pl.title(r'Zoomed $\Phi_{PEG1k}^{in}$ for varied $\Phi_{PEG100}^{out}$ fixed $\Phi_{PEG1k}^{out}$ $\Phi_{PEG10k}^{out\ only}$')
			pl.xlabel(r'$\Delta f$', size = 'x-large')
			pl.ylabel(r'$\Phi_{m}^{i}$($\Phi_{s}^{i,o}$,$\Phi_{m}^{o}$,$\Phi_{b}^{o}$,$\Delta f$)', size = 'x-large')
			pl.axis([0.0, 0.4, 0.0, 0.0004])
#			pl.legend(loc='upper right')


			pl.figure(13)
#			pl.plot(phis, Mus_s_in,marker=markers[h],color = colors[j],linestyle = lts[i],linewidth= 0.7,markersize=4,markevery=7,markerfacecolor='w')
			pl.plot(phis, Mus_in_s_only,marker=markers[h],color = 'k',linestyle = lts[i],linewidth= 0.7,markersize=4,markevery=7,markerfacecolor='w')
			pl.title(r'Chemical Potential of small polymers in pore with and without medium polymers allowed in pore$')
		#	pl.legend(loc='upper right')
		
		#	pl.figure()
		#	pl.plot(-0.005,-0.005,'k',marker = 'x',linewidth = 0.0, label=r'$\Phi_{s} = 0.10$')#color = colors[j],
		#	pl.plot(-0.005,-0.005,'k',marker = '^',linewidth = 0.0, label=r'$\Phi_{s} = 0.20$')#color = colors[j],
		#	pl.plot(-0.005,-0.005,'k',marker = 'x',linewidth = 0.0, label=r'$\Phi_{m} = 0.15$')#color = colors[j],
		#	pl.plot(-0.005,-0.005,'k',marker = '^',linewidth = 0.0, label=r'$\Phi_{m} = 0.20$')#color = colors[j],
		#	pl.plot(-0.005,-0.005,'k',linestyle = ':', label=r'$\Phi_{b} = 0.05$')
		#	pl.plot(-0.005,-0.005,'k',linestyle = '--',label=r'$\Phi_{b} = 0.10$')
		#	pl.plot(-0.005,-0.005,color = 'b',label=r'$\Delta f = N*0.75$')#color = colors[j],
		#	pl.plot(-0.005,-0.005,color = 'g',label=r'$\Delta f = N*1.00$')#color = colors[j],
		#	pl.plot(-0.005,-0.005,color = 'r',label=r'$\Delta f = N*1.25$')#color = colors[j],
			pl.figure(14)
			pl.plot(phis, Mus_m_in,marker=markers[h],color = colors[j],linestyle = lts[i],linewidth= 0.7,markersize=4,markevery=7,markerfacecolor='w')
			pl.title(r'Chemical Potential of medium polymers in pore with small polymers allowed in pore$')
		#	pl.legend(loc='upper right')

			pl.figure(15)
			pl.plot(phis, Mus_m_s_in_diff,marker=markers[h],color = colors[j],linestyle = lts[i],linewidth= 0.7,markersize=4,markevery=7,markerfacecolor='w')
			pl.title(r'chemical potential of medium polymers in pore minus chemical potential of small polymers in pore$')
#
#
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
pl.figure(16)

pl.title(r'Legend for plots 1-6')
pl.plot(-0.005,-0.005, color = 'b',linewidth = 0.7, label=r'$\phi_{s} = 0.10$')#color = colors[j],
pl.plot(-0.005,-0.005, color = 'g',linewidth = 0.7, label=r'$\phi_{s} = 0.20$')#color = colors[j],
pl.plot(-0.005,-0.005,'k',marker = 'x',linewidth = 0.0, label=r'$\phi_{m} = 0.15$')#color = colors[j],
pl.plot(-0.005,-0.005,'k',marker = '*',linewidth = 0.0, label=r'$\phi_{m} = 0.20$')#color = colors[j],
pl.plot(-0.005,-0.005,'k',linestyle = ':', label=r'$\phi_{b} = 0.05$')
pl.plot(-0.005,-0.005,'k',linestyle = '--',label=r'$\phi_{b} = 0.10$')
pl.legend(loc='center')

pl.figure(17)

pl.title(r'Legend for plots 7-15')
pl.plot(0.005,0.005, color = 'b',linewidth = 0.7, label=r'df = N*0.75')#color = colors[j],
pl.plot(0.005,0.005, color = 'g',linewidth = 0.7, label=r'df = N*1.00')#color = colors[j],
pl.plot(0.005,0.005, color = 'r',linewidth = 0.7, label=r'df = N*1.25')#color = colors[j],
pl.plot(0.005,0.005,'k',marker = 'x',linewidth = 0.0, label=r'$\phi_{m} = 0.15$')#color = colors[j],
pl.plot(0.005,0.005,'k',marker = '*',linewidth = 0.0, label=r'$\phi_{m} = 0.20$')#color = colors[j],
pl.plot(0.005,0.005,'k',linestyle = ':', label=r'$\phi_{b} = 0.05$')
pl.plot(0.005,0.005,'k',linestyle = '--',label=r'$\phi_{b} = 0.10$')
pl.legend(loc='center')

#	pl.figure(16)
		#	pl.semilogy(phis,Mus_in_diff,marker=markers[h],color=colors[j],linestyle = lts[i],linewidth= 0.7,markersize=4,markevery=7,markerfacecolor='w')
		#	pl.title(r'chemical potential of small polymers in pore with mediums minus chemical potential of small polymers only allowed in pore$')
#pl.axis([0.0, 0.4, 0.0, 0.0004])
#pl.figure(17)
#
#pl.plot(-0.005,-0.005,'k',marker = 'x',linewidth = 0.0, label=r'$\Phi_{s} = 0.10$')#color = colors[j],
#pl.plot(-0.005,-0.005,'k',marker = '^',linewidth = 0.0, label=r'$\Phi_{s} = 0.20$')#color = colors[j],
#
#pl.plot(-0.005,-0.005,'k',marker = 'x',linewidth = 0.0, label=r'$\Phi_{m} = 0.15$')#color = colors[j],
#pl.plot(-0.005,-0.005,'k',marker = '^',linewidth = 0.0, label=r'$\Phi_{m} = 0.20$')#color = colors[j],
#
#pl.plot(-0.005,-0.005,'k',linestyle = ':', label=r'$\Phi_{b} = 0.05$')
#pl.plot(-0.005,-0.005,'k',linestyle = '--',label=r'$\Phi_{b} = 0.10$')
#
#pl.plot(-0.005,-0.005,color = 'b',label=r'$\Delta f = N*0.75$')#color = colors[j],
#pl.plot(-0.005,-0.005,color = 'g',label=r'$\Delta f = N*1.00$')#color = colors[j],
#pl.plot(-0.005,-0.005,color = 'r',label=r'$\Delta f = N*1.25$')#color = colors[j],
#pl.legend(loc='center')
#
#pl.title(r'Legend', size='x-large')
#			#pl.xlabel(r'$\Delta f$', size = 'x-large')
#			#pl.ylabel(r'$\Phi_{m}^{i}$($\Phi_{s}^{i,o}$,$\Phi_{m}^{o}$,$\Phi_{b}^{o}$,$\Delta f$)', size = 'x-large')
#			pl.axis([0.0, 1, 0, 1])
#

#pl.legend(loc='upper right')
#pl.show()
pp.savefig(7)
pp.savefig(8)
pp.savefig(9)
pp.savefig(10)
pp.savefig(11)
pp.savefig(12)
#pp.savefig(13)
pp.savefig(14)
pp.savefig(15)
pp.savefig(16)
pp.savefig(17)
pp.close()



