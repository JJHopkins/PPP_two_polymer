#!/usr/bin/python

from matplotlib import pyplot as pl
import numpy as np

# Define Constants
Vbar = 1.0 # mL/g
Nb = 35000. / 400 * 9 # Monomers / big molecule
Ns = 1000. / 400 * 9 # Monomers / small molecule
alpha = 0.49

# Define Osmotic Pressure equation as function of short and big fraction
def osmotic_pressure(s, b,k):
#	return s/Ns + b/Nb + (5./4)*alpha*(s + b)**(9./4)
	return (s/Ns + b/Nb + (5./4)*alpha*(s + b)**(9./4))-(0./Ns + k/Nb + (5./4)*alpha*(0. + k)**(9./4))

# Let short polymer volume fraction range from 0 to 60 percent 
phis = np.linspace(0, 0.90)
# Let big polymer volume fraction range from 0 to 30 percent
phib30 = np.linspace(0, 0.30,num=500)
phib20 = np.linspace(0, 0.20,num=500)
phib10 = np.linspace(0, 0.10,num=500)
#phi = np.linspace(0,1,num=500)
#pl.semilogy(phi, osmotic_pressure(0, phi), label='hi')

pl.figure()

# Plot 0 to 30% big polymer and 0 short polymer
#pl.semilogy(phib10, osmotic_pressure(0.0, phib10),'b,', label=' 0% PEG1000')
#pl.semilogy(phib20, osmotic_pressure(0.0, phib20),'g,', label=' 0% PEG1000')
#pl.semilogy(phib30, osmotic_pressure(0.0, phib30),'r,', label=' 0% PEG1000')


# Plot using whole range of phis and  10,20,30% big polymer
#pl.semilogy(phis+0.10, osmotic_pressure(phis, 0.10),'b-', label='10% PEG3500')
#pl.semilogy(phis+0.20, osmotic_pressure(phis, 0.20),'g-', label='20% PEG3500')
#pl.semilogy(phis+0.30, osmotic_pressure(phis, 0.30),'r-', label='30% PEG3500')

pl.semilogy(phis+0.10, osmotic_pressure(phis, 0.10, 0.10),'b-', label='10% PEG3500')
pl.semilogy(phis+0.20, osmotic_pressure(phis, 0.20, 0.20),'g-', label='20% PEG3500')
pl.semilogy(phis+0.30, osmotic_pressure(phis, 0.30, 0.30),'r-', label='30% PEG3500')


pl.legend(loc='lower right')
pl.xlabel('total polymer concentration')
pl.ylabel('Osmotic Pressure')
pl.title('Osmotic Pressure of Mixed Polymers')
pl.axis([0, 1.0, 10**(-3), 10**(0)])

pl.savefig('Osmotic_Pressure.png')
pl.show()

