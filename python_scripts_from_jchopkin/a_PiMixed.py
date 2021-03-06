#!/usr/bin/python

from matplotlib import pyplot as pl
import numpy as np

# Define Constants
Vbar = 1.0 # mL/g
Nb = 35000. / 400 * 9
Ns = 1000. / 400 * 9
alpha = 0.49

# Define Osmotic Pressure equation as function of short and big fraction
def osmotic_pressure(phis, phib):
	return phis/Ns + phib/Nb + 5./4*alpha*(phis + phib)**(9./4)
# Let short polymer volume fraction range from 0 to 60 percent 
phis = np.linspace(0, 60)
# Let big polymer volume fraction range from 0 to 30 percent
phib = np.linspace(0, 30)
#phib20 = np.linspace(0, 20)
#phib10 = np.linspace(0, 10)


# Plot using whole range of phis and  10,20,10% big polymer
#pl.loglog(phis, osmotic_pressure(phis, 0), label='00%')
pl.figure(1)
pl.subplots_adjust(wspace=0.001)
pl. subplot(122)
pl.semilogy(phis, osmotic_pressure(phis, 10), label='10% PEG3500')
pl.semilogy(phis, osmotic_pressure(phis, 20), label='20% PEG3500')
pl.semilogy(phis, osmotic_pressure(phis, 30), label='30% PEG3500')

pl.legend(loc='lower right')
pl.xlabel('Numb % added PEG1000')
#pl.ylabel('Osmotic Pressure')
pl.title('of Mixed Polymers')
pl.axis([0, 70, 10**(-1), 10**(5)])

# Plot 0 to 30% big polymer and 0 short polymer
pl.subplot(121)
pl.semilogy(phib, osmotic_pressure(0, phib),'r--', label=' 0% PEG1000')
#pl.semilogy(phib20, osmotic_pressure(0, phib20),'g--', label=' 0% PEG1000'
#pl.semilogy(phib10, osmotic_pressure(0, phib10),'b--', label=' 0% PEG1000'
pl.legend(loc='lower right')
pl.xlabel('Numb % PEG3500')
pl.ylabel('ln(Osmotic Pressure)')
pl.title('Osmotic Pressure ')
pl.axis([0, 30, 10**(-1), 10**(5)])

pl.show()

