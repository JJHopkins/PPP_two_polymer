#!/usr/bin/python

from matplotlib import pyplot as pl
import numpy as np
from pylab import * 
from numpy import log,exp

# Define Constants
Vbar = 1.0 # mL/g
Nb = 10000. / 400 * 9 # Monomers / big molecule
Nm = 1000. / 400 * 9 # Monomers / med molecule
Ns = 100. / 400 * 9 # Monomers / small molecule
alpha = 0.49 # alpha_tilda= alpha*(1/2 - chi)**(3/4)

# Define Osmotic Pressure equation as function of short and big fraction
def osmotic_pressure(s, m, b):
	#return s/Ns + b/Nb + (5./4)*alpha*(s + b)**(9./4)
	return - log(1 - s - m - b) + (1 - s - m - b) - 1 + s/Ns + m/Nm + b/Nb - \
		(1./2)*(1 - (1 - s - m - b))**2 + (5./4)*alpha* \
		(1 - (1 - s - m - b))**(9./4)

# Let short polymer volume fraction range from 0 to 60 percent 
phis = np.linspace(0, 0.60,60)
phim = np.linspace(0, 0.60,60)
# Let big polymer volume fraction range from 0 to 30 percent

phis05 = np.linspace(0, 0.15,15)
phis10 = np.linspace(0, 0.10,10)
phis15 = np.linspace(0, 0.05,5)

phim05 = np.linspace(0, 0.15,15)
phim10 = np.linspace(0, 0.10,10)
phim15 = np.linspace(0, 0.05,5)

phib05 = np.linspace(0, 0.15,15)
phib10 = np.linspace(0, 0.10,10)
phib15 = np.linspace(0, 0.05,5)

phib10 = np.linspace(0, 0.10,10)
phib20 = np.linspace(0, 0.20,20)
phib30 = np.linspace(0, 0.30,30)
#pl.figure(1)
#pl.axis([-0.01,1.0,-0.0,0.1])
# Plot 0 to 30% big polymer and 0 short polymer
#pl.plot(log(phib10), log(osmotic_pressure(0.0, phib10)),'b,', label=' 0% PEG400')
#pl.plot(log(phib20), log(osmotic_pressure(0.0, phib20)),'g,', label=' 0% PEG400')
#pl.plot(log(phib30), log(osmotic_pressure(0.0, phib30)),'r,', label=' 0% PEG400')

# Plot using whole range of phis and  10,20,30% big polymer
#pl.plot(log(phis+0.10), log(osmotic_pressure(phis, 0.10)),'b-', label='10% PEG3500')
#pl.plot(log(phis+0.20), log(osmotic_pressure(phis, 0.20)),'g-', label='20% PEG3500')
#pl.plot(log(phis+0.30), log(osmotic_pressure(phis, 0.30)),'r-', label='30% PEG3500')
#pl.axis([-3.0,-0.5,-6.0,-1.0])

#pl.figure(1)
#pl.semilogy(phib10, osmotic_pressure(0.0, phib10),'b,', label=' 0% PEG400')
#pl.semilogy(phis+0.10, osmotic_pressure(phis, 0.10),'b-', label='10% PEG3500')
#pl.legend(loc='lower right')
#pl.xlabel('Total Polymer Number Fraction')
#pl.ylabel('Osmotic Pressure   (dimensionless)')
#pl.title('Osmotic Pressure of Mixed Polymers')
#pl.axis([0.05,0.8,0.0005,0.8])
#pl.savefig('111207_PiMixed_10.eps', dpi=600)
#pl.show()
#
#pl.figure(2)
#pl.semilogy(phib20, osmotic_pressure(0.0, phib20),'g,', label=' 0% PEG400')
#pl.semilogy(phis+0.20, osmotic_pressure(phis, 0.20),'g-', label='20% PEG3500')
#pl.legend(loc='lower right')
#pl.xlabel('Total Polymer Number Fraction')
#pl.ylabel('Osmotic Pressure   (dimensionless)')
#pl.title('Osmotic Pressure of Mixed Polymers')
#pl.axis([0.05,0.8,0.0005,0.8])
#pl.savefig('111207_PiMixed_20.eps', dpi=600)
#pl.show()
#
#pl.figure(3)
#pl.semilogy(phib30, osmotic_pressure(0.0, phib30),'r,', label=' 0% PEG400')
#pl.semilogy(phis+0.30, osmotic_pressure(phis, 0.30),'r-', label='30% PEG3500')
#pl.legend(loc='lower right')
#pl.xlabel('Total Polymer Number Fraction')
#pl.ylabel('Osmotic Pressure   (dimensionless)')
#pl.title('Osmotic Pressure of Mixed Polymers')
#pl.axis([0.05,0.8,0.0005,0.8])
#pl.savefig('111207_PiMixed_30.eps', dpi=600)
#pl.show()


pl.figure(4)# Plot using whole range of phis and  10,20,30% big polymer
ax=axes()
pl.semilogy(phis+0.10, osmotic_pressure(phis, 0.05, 0.05),'b-',markerfacecolor='w',markevery=4,linewidth=0.5, markersize=3, label='$\phi_m=0.05,\phi_b=0.05$')
pl.semilogy(phis+0.20, osmotic_pressure(phis, 0.10, 0.10),'b-',markerfacecolor='w',markevery=4,linewidth=0.5, markersize=4, label='$\phi_m=0.10,\phi_b=0.10$')
pl.semilogy(phis+0.30, osmotic_pressure(phis, 0.15, 0.15),'b-',markerfacecolor='w',markevery=4,linewidth=0.5, markersize=4, label='$\phi_m=0.15,\phi_b=0.15$')

pl.semilogy(phis+0.10, osmotic_pressure(phis, 0.0, 0.10),'g-',markerfacecolor='w',markevery=4,linewidth=0.5, markersize=3, label='$\phi_m=0.0, \phi_b= 0.10$')
pl.semilogy(phis+0.20, osmotic_pressure(phis, 0.0, 0.20),'g-',markerfacecolor='w',markevery=4,linewidth=0.5, markersize=4, label='$\phi_m=0.0, \phi_b= 0.20$')
pl.semilogy(phis+0.30, osmotic_pressure(phis, 0.0, 0.30),'g-',markerfacecolor='w',markevery=4,linewidth=0.5, markersize=4, label='$\phi_m=0.0, \phi_b= 0.30$')

pl.semilogy(phib05 + phim05, osmotic_pressure(0.0, phim05, phib05),'b--',markerfacecolor='w',linewidth=0.5, markersize=2)#, label=' $\Phi_s$0% PEG400')
pl.semilogy(phib10 + phim10, osmotic_pressure(0.0, phim10, phib10),'b--',markerfacecolor='w',linewidth=0.5, markersize=2)#, label=' 0% PEG400')
pl.semilogy(phib15 + phim15, osmotic_pressure(0.0, phim15, phib15),'b--',markerfacecolor='w',linewidth=0.5, markersize=2)#, label=' 0% PEG400')

pl.semilogy(phib10, osmotic_pressure(0.0, 0.0, phib10),'g--',markerfacecolor='w',linewidth=0.5, markersize=2)#, label=' $\Phi_s$0% PEG400')
pl.semilogy(phib20, osmotic_pressure(0.0, 0.0, phib20),'g--',markerfacecolor='w',linewidth=0.5, markersize=2)#, label=' 0% PEG400')
pl.semilogy(phib30, osmotic_pressure(0.0, 0.0, phib30),'g--',markerfacecolor='w',linewidth=0.5, markersize=2)#, label=' 0% PEG400')
#idx=range(0, len(x),5)
#plot(x[idx],y[idx],'ro')

pl.legend(loc='lower right')
pl.xlabel(r'$\phi_{s} + \phi_{m} + \phi_{b}$', size= 'large')# size= 'x-large')
pl.ylabel(r'$\~ {\Pi}$', size= 'x-large')
pl.title('Osmotic Pressure of Three PEG Sizes in Water,$\phi_{s}$ Varied ')
ax.set_xticks([0.01,0.1,0.5,1])
ax.set_xticklabels(['0.01','0.1','0.5','1'])
ax.grid()
pl.axis([0.09,0.57,0.003,0.35])
pl.minorticks_on()
pl.savefig('121121_OP_3_poly_semilog_phis_varied.eps', dpi=600)

pl.figure(5)# Plot using whole range of phis and  10,20,30% big polymer
ax=axes()
pl.semilogy(phim+0.10, osmotic_pressure(0.05, phim, 0.05),'b-',markerfacecolor='w',markevery=4,linewidth=0.5, markersize=3, label='$\phi_s=\phi_b=0.05$')
pl.semilogy(phim+0.20, osmotic_pressure(0.10, phim, 0.10),'b-',markerfacecolor='w',markevery=4,linewidth=0.5, markersize=4, label='$\phi_s=\phi_b=0.10$')
pl.semilogy(phim+0.30, osmotic_pressure(0.15, phim, 0.15),'b-',markerfacecolor='w',markevery=4,linewidth=0.5, markersize=4, label='$\phi_s=\phi_b=0.15$')
                                             
pl.semilogy(phim+0.10, osmotic_pressure(0.0, phim, 0.10),'g-',markerfacecolor='w',markevery=4,linewidth=0.5, markersize=3, label='$\phi_s=0,\phi_b= 0.10$')
pl.semilogy(phim+0.20, osmotic_pressure(0.0, phim, 0.20),'g-',markerfacecolor='w',markevery=4,linewidth=0.5, markersize=4, label='$\phi_s=0,\phi_b= 0.20$')
pl.semilogy(phim+0.30, osmotic_pressure(0.0, phim, 0.30),'g-',markerfacecolor='w',markevery=4,linewidth=0.5, markersize=4, label='$\phi_s=0,\phi_b= 0.30$')

pl.semilogy(phib05 + phis05, osmotic_pressure(phis05, 0.0, phib05),'b--',markerfacecolor='w',linewidth=0.5, markersize=2)#, label='$\phi_m=0,\phi_s=\phi_b=0.05$')#, label=' $\Phi_s$0% PEG400')
pl.semilogy(phib10 + phis10, osmotic_pressure(phis10, 0.0, phib10),'b--',markerfacecolor='w',linewidth=0.5, markersize=2)#, label='$\phi_m=0,\phi_s=\phi_b=0.10$')#, label=' 0% PEG400')
pl.semilogy(phib15 + phis15, osmotic_pressure(phis15, 0.0, phib15),'b--',markerfacecolor='w',linewidth=0.5, markersize=2)#, label='$\phi_m=0,\phi_s=\phi_b=0.15$')#, label=' 0% PEG400')

pl.semilogy(phib10, osmotic_pressure(0.0,0.0, phib10),'g--',markerfacecolor='w',linewidth=0.5, markersize=2)#, label='$\phi_s=\phi_m=0,\phi_b=0.10$')#, label=' $\Phi_s$0% PEG400')
pl.semilogy(phib20, osmotic_pressure(0.0,0.0, phib20),'g--',markerfacecolor='w',linewidth=0.5, markersize=2)#, label='$\phi_s=\phi_m=0,\phi_b=0.20$')#, label=' 0% PEG400')
pl.semilogy(phib30, osmotic_pressure(0.0,0.0, phib30),'g--',markerfacecolor='w',linewidth=0.5, markersize=2)#, label='$\phi_s=\phi_m=0,\phi_b=0.30$')#, label=' 0% PEG400')
#idx=range(0, len(x),5)
#plot(x[idx],y[idx],'ro')

pl.legend(loc='lower right')
pl.xlabel(r'$\phi_{s} + \phi_{m} + \phi_{b}$', size= 'large')# size= 'x-large')
pl.ylabel(r'$\~ {\Pi}$', size= 'x-large')
pl.title(r'Osmotic Pressure of Three PEG Sizes in Water,$\phi_{m}$ Varied ')
ax.set_xticks([0.01,0.1,0.5,1])
ax.set_xticklabels(['0.01','0.1','0.5','1'])
ax.grid()
pl.axis([0.09,0.57,0.003,0.35])
pl.minorticks_on()
pl.savefig('121121_OP_3_poly_semilog_phim_varied.eps', dpi=600)



pl.figure(6)# Plot using whole range of phis and  10,20,30% big polymer
ax=axes()
pl.loglog(phis+0.10, osmotic_pressure(phis, 0.05, 0.05),'b-',markerfacecolor='w',markevery=4,linewidth=0.5, markersize=3, label='$\phi_m=0.05,\phi_b=0.05$')
pl.loglog(phis+0.20, osmotic_pressure(phis, 0.10, 0.10),'b-',markerfacecolor='w',markevery=4,linewidth=0.5, markersize=4, label='$\phi_m=0.10,\phi_b=0.10$')
pl.loglog(phis+0.30, osmotic_pressure(phis, 0.15, 0.15),'b-',markerfacecolor='w',markevery=4,linewidth=0.5, markersize=4, label='$\phi_m=0.15,\phi_b=0.15$')
         
pl.loglog(phis+0.10, osmotic_pressure(phis, 0.0, 0.10),'g-',markerfacecolor='w',markevery=4,linewidth=0.5, markersize=3, label='$\phi_m=0.0, \phi_b= 0.10$')
pl.loglog(phis+0.20, osmotic_pressure(phis, 0.0, 0.20),'g-',markerfacecolor='w',markevery=4,linewidth=0.5, markersize=4, label='$\phi_m=0.0, \phi_b= 0.20$')
pl.loglog(phis+0.30, osmotic_pressure(phis, 0.0, 0.30),'g-',markerfacecolor='w',markevery=4,linewidth=0.5, markersize=4, label='$\phi_m=0.0, \phi_b= 0.30$')
         
pl.loglog(phib05 + phim05, osmotic_pressure(0.0, phim05, phib05),'b--',markerfacecolor='w',linewidth=0.5, markersize=2)#, label=' $\Phi_s$0% PEG400')
pl.loglog(phib10 + phim10, osmotic_pressure(0.0, phim10, phib10),'b--',markerfacecolor='w',linewidth=0.5, markersize=2)#, label=' 0% PEG400')
pl.loglog(phib15 + phim15, osmotic_pressure(0.0, phim15, phib15),'b--',markerfacecolor='w',linewidth=0.5, markersize=2)#, label=' 0% PEG400')
         
pl.loglog(phib10, osmotic_pressure(0.0, 0.0, phib10),'g--',markerfacecolor='w',linewidth=0.5, markersize=2)#, label=' $\Phi_s$0% PEG400')
pl.loglog(phib20, osmotic_pressure(0.0, 0.0, phib20),'g--',markerfacecolor='w',linewidth=0.5, markersize=2)#, label=' 0% PEG400')
pl.loglog(phib30, osmotic_pressure(0.0, 0.0, phib30),'g--',markerfacecolor='w',linewidth=0.5, markersize=2)#, label=' 0% PEG400')
#idx=range(0, len(x),5)
#plot(x[idx],y[idx],'ro')

pl.legend(loc='lower right')
pl.xlabel(r'$\phi_{s} + \phi_{m} + \phi_{b}$', size= 'large')# size= 'x-large')
pl.ylabel(r'$\~ {\Pi}$', size= 'x-large')
pl.title('Osmotic Pressure of Three PEG Sizes in Water,$\phi_{s}$ Varied ')
ax.set_xticks([0.01,0.1,0.5,1])
ax.set_xticklabels(['0.01','0.1','0.5','1'])
ax.grid()
pl.axis([0.09,0.57,0.003,0.35])
pl.minorticks_on()
pl.savefig('121121_OP_3_poly_loglog_phis_varied.eps', dpi=600)

pl.figure(7)# Plot using whole range of phis and  10,20,30% big polymer
ax=axes()
pl.loglog(phim+0.10, osmotic_pressure(0.05, phim, 0.05),'b-',markerfacecolor='w',markevery=4,linewidth=0.5, markersize=3, label='$\phi_s=\phi_b=0.05$')
pl.loglog(phim+0.20, osmotic_pressure(0.10, phim, 0.10),'b-',markerfacecolor='w',markevery=4,linewidth=0.5, markersize=4, label='$\phi_s=\phi_b=0.10$')
pl.loglog(phim+0.30, osmotic_pressure(0.15, phim, 0.15),'b-',markerfacecolor='w',markevery=4,linewidth=0.5, markersize=4, label='$\phi_s=\phi_b=0.15$')
                                 
pl.loglog(phim+0.10, osmotic_pressure(0.0, phim, 0.10),'g-',markerfacecolor='w',markevery=4,linewidth=0.5, markersize=3, label='$\phi_s=0,\phi_b= 0.10$')
pl.loglog(phim+0.20, osmotic_pressure(0.0, phim, 0.20),'g-',markerfacecolor='w',markevery=4,linewidth=0.5, markersize=4, label='$\phi_s=0,\phi_b= 0.20$')
pl.loglog(phim+0.30, osmotic_pressure(0.0, phim, 0.30),'g-',markerfacecolor='w',markevery=4,linewidth=0.5, markersize=4, label='$\phi_s=0,\phi_b= 0.30$')
   
pl.loglog(phib05 + phis05, osmotic_pressure(phis05, 0.0, phib05),'b--',markerfacecolor='w',linewidth=0.5, markersize=2)#, label='$\phi_m=0,\phi_s=\phi_b=0.05$')
pl.loglog(phib10 + phis10, osmotic_pressure(phis10, 0.0, phib10),'b--',markerfacecolor='w',linewidth=0.5, markersize=2)#, label='$\phi_m=0,\phi_s=\phi_b=0.10$')
pl.loglog(phib15 + phis15, osmotic_pressure(phis15, 0.0, phib15),'b--',markerfacecolor='w',linewidth=0.5, markersize=2)#, label='$\phi_m=0,\phi_s=\phi_b=0.15$')
  
pl.loglog(phib10, osmotic_pressure(0.0,0.0, phib10),'g--',markerfacecolor='w',linewidth=0.5, markersize=2)#, label='$\phi_s=\phi_m=0,\phi_b=0.10$')
pl.loglog(phib20, osmotic_pressure(0.0,0.0, phib20),'g--',markerfacecolor='w',linewidth=0.5, markersize=2)#, label='$\phi_s=\phi_m=0,\phi_b=0.20$')
pl.loglog(phib30, osmotic_pressure(0.0,0.0, phib30),'g--',markerfacecolor='w',linewidth=0.5, markersize=2)#, label='$\phi_s=\phi_m=0,\phi_b=0.30$')
#idx=range(0, len(x),5)
#plot(x[idx],y[idx],'ro')

pl.legend(loc='lower right')
pl.xlabel(r'$\phi_{s} + \phi_{m} + \phi_{b}$', size= 'large')# size= 'x-large')
pl.ylabel(r'$\~ {\Pi}$', size= 'x-large')
pl.title(r'Osmotic Pressure of Three PEG Sizes in Water,$\phi_{m}$ Varied ')
ax.set_xticks([0.01,0.1,0.5,1])
ax.set_xticklabels(['0.01','0.1','0.5','1'])
ax.grid()
pl.axis([0.09,0.57,0.003,0.35])
pl.minorticks_on()
pl.savefig('121121_OP_3_poly_loglog_phim_varied.eps', dpi=600)
#pylab.rcParams.update

#f=pl.figure(5)
## Plot using whole range of phis and  10,20,30% big polymer
#ax=axes()
#ax.loglog(phib10, osmotic_pressure(0.0, phib10),'k--',markerfacecolor='w',linewidth=0.5, markersize=2)#,'b,', label=' 0% PEG400')
#ax.loglog(phib20, osmotic_pressure(0.0, phib20),'k--',markerfacecolor='w',linewidth=0.5, markersize=2)#,'g,', label=' 0% PEG400')
#pl.loglog(phib30, osmotic_pressure(0.0, phib30),'k--',markerfacecolor='w',linewidth=0.5, markersize=2)#,'r,', label=' 0% PEG400')
#ax.loglog(phis+0.10, osmotic_pressure(phis, 0.10),'k-s',markerfacecolor='w',markevery=4,linewidth=0.5, markersize=3, label='$\phi_b= 0.10$')#'b-', label='10% PEG3500')
#ax.loglog(phis+0.20, osmotic_pressure(phis, 0.20),'k-^',markerfacecolor='w',markevery=4,linewidth=0.5, markersize=4, label='$\phi_b= 0.20$')#'g-', label='20% PEG3500')
#ax.loglog(phis+0.30, osmotic_pressure(phis, 0.30),'k-d',markerfacecolor='w',markevery=4,linewidth=0.5, markersize=4, label='$\phi_b= 0.30$')#'r-', label='30% PEG3500')
#
#
#
##phi = np.linspace(0,1)
##pl.semilogy(phi, osmotic_pressure(0, phi), label='hi')
#
##pl.legend(loc='lower right')
#ax.set_xlabel(r'$\phi_{s} + \phi_{b}$', size= 'x-large')
##pl.xlabel('Total Polymer Number Fraction')
#ax.set_ylabel(r'$\~ {\Pi}$', size= 'x-large')
##pl.ylabel('Osmotic Pressure   (dimensionless)')
##pl.title('Osmotic Pressure of Mixed Polymers')
#ax.set_xticks([0.01,0.1,0.5,1])
#ax.set_xticklabels(['0.01','0.1','0.5','1'])
#ax.grid()
#
##ax.grid(True)
#pl.axis([0.09,0.57,0.003,0.35])
#
## --- Set xtick labels:
##ax.set_xticklabels(['0.01','','','','0.1'])
#
##ax.set_xlim(0.09,0.57)
##ax.set_ylim(
#pl.savefig('121121_OP_3_poly.eps', dpi=600)
pl.show()



