#!/usr/bin/python

from matplotlib import pyplot as pl
import numpy as np
from numpy import log,exp

# Define Constants
Vbar = 1.0 # mL/g
Nb = 35000. / 400 * 9 # Monomers / big molecule
Ns = 400. / 400 * 9 # Monomers / small molecule
alpha_tilda = 0.49 # alpha_tilda= alpha*(1/2 - chi)**(3/4)
phis = np.linspace(0.01, 0.55, 10)
phib = np.linspace(0.15, 0.01,10)

phi_ss = np.linspace(0.01, 0.45, 10)
phibb = np.linspace(0.01, 0.15,10)

def osmotic_pressure(s, b):
	return - log(1 - s - b) + (1 - s - b) - 1 + s/Ns + b/Nb - \
		(1./2)*(1 - (1 - s - b))**2 + (5./4)*alpha_tilda* \
		(1 - (1 - s - b))**(9./4)
#def total(x,y):
#	return x + y

pl.figure()
pl.plot(phis+phib, (osmotic_pressure(phis,phib)), 'b')
#pl.plot(phi_ss+phibb, (osmotic_pressure(phi_ss, 0.0)), 'r')
#pl.plot(phi_ss+phibb, (osmotic_pressure(0.0,phibb)), 'g')
#pl.plot(phi_ss+phibb, (osmotic_pressure(phi_ss,phibb)), 'c')
#pl.plot(phi_ss+phibb, (osmotic_pressure(phi_ss, 0.0)) + (osmotic_pressure(0.0,phibb)), 'k')
#pl.plot(phibb,total)

pl.xlabel('Polymer Number Fraction')
pl.ylabel('(Osmotic Pressure)   (dimensionless)')
pl.title('Osmotic Pressure')

#pl.axis([-0.01,0.8,0.0001,0.8])
pl.savefig('Big_into_Short_OS.eps', dpi=600)
pl.show()


