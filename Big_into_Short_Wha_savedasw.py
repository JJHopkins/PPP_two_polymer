#!/usr/bin/python

from matplotlib import pyplot as pl
import numpy as np
from numpy import log,exp

# Define Constants
Vbar = 1.0 # mL/g
Nb = 35000. / 400 * 9 # Monomers / big molecule
Ns = 400. / 400 * 9 # Monomers / small molecule
alpha_tilda = 0.49 # alpha_tilda= alpha*(1/2 - chi)**(3/4)
phi_ss = np.linspace(0.05, 0.30 10)
phibb = np.linspace(0.20, 0.05,10)

def osmotic_pressure(s, b):
	return - log(1 - s - b) + (1 - s - b) - 1 + s/Ns + b/Nb - \
		(1./2)*(1 - (1 - s - b))**2 + (5./4)*alpha_tilda* \
		(1 - (1 - s - b))**(9./4)
#def total(x,y):
#	return x + y

# Plot 0 to 30% big polymer and 0 short polymer
#pl.plot(log(phib10), log(osmotic_pressure(0.0, phib10)),'b,', label=' 0% PEG400')
#pl.plot(log(phib20), log(osmotic_pressure(0.0, phib20)),'g,', label=' 0% PEG400')
#pl.plot(log(phib30), log(osmotic_pressure(0.0, phib30)),'r,', label=' 0% PEG400')

# Plot using whole range of phis and  10,20,30% big polymer
#pl.plot(phi_ss,osmotic_pressure(phi_ss, 0.0))
#pl.plot(log(phis+0.20), log(osmotic_pressure(phis, 0.20)),'g-', label='20% PEG3500')
#pl.plot(log(phis+0.30), log(osmotic_pressure(phis, 0.30)),'r-', label='30% PEG3500')
#pl.axis([-3.0,-0.5,-6.0,-1.0])

pl.figure()
pl.plot(phi_ss, log(osmotic_pressure(phi_ss,phibb)))
pl.plot(phibb, log(osmotic_pressure(0.0,phibb)))
pl.plot(phi_ss, log(osmotic_pressure(phi_ss, 0.0)))
#pl.plot(phibb,total)
pl.xlabel('Polymer Number Fraction')
pl.ylabel('log(Osmotic Pressure)   (dimensionless)')
pl.title('Osmotic Pressure')

#pl.axis([-0.01,0.8,0.0001,0.8])
pl.savefig('Big_into_Short_OS.eps', dpi=600)
pl.show()


