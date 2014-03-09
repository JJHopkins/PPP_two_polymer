#!/usr/bin/python

from matplotlib import pyplot as pl
import numpy as np
from math import e

# Define Constants

#Nb = 35000. / 400 * 9 # Monomers / big molecule
Ns = 1000. / 400 * 9 # Monomers / small molecule
k = 8.617*10**(-5) # (eV)/(K) 1.380*10**(-23) J/K
T = 300. # Kelvin
alpha = 0.49
chi = 0.50

##Define mu_inside and mu_outside as their own fncs of phi_s_outside
#def mso(so):
#    return (1. - so)+ (5./4)*alpha*(so)**(9./4) - (9./4)*(alpha*(so)**(5./4) - chi*(1 - so)**2)
#def msi(so):
#    return (1. - (1. - so)) + (5./4)*alpha*(1 - so)**(9./4) - (9./4)*(alpha*(1. - so)**(5./4) - chi*(1. - (1. - so))**2)
#def chem_pot(so,f):
#	return e**((f) - so + (1 - so) - (mso(so) + msi(so))*Ns)   



# define partition_coefficent= exp[-(delta mu_s-delta f)/kT]   ---->  exp[ (1/kT)*(delta f - (phi_s_out-phi_s_in + mu_bar_s_in-mu_bar_s_out)*Ns)] -----> exp[ (1/kT)*(delta f - (so-(1-so) + msi(so)-mso(so))*Ns)]
# JCH don't forget to put (1./(k*T))* back in, gave "overflow encountered in power"error
def part_coeff(so,f):
    return e**(((f) - so + (1. - so) - (((1. - (1. - so)) + (5./4)*alpha*(1. - so)**(9./4) - (9./4)*(alpha*(1. - so)**(5./4) - chi*(1. - (1. - so))**2)) - ((1. - so)+ (5./4)*alpha*(so)**(9./4) - (9./4)*(alpha*(so)**(5./4) - chi*(1. - so)**2)))*Ns))


pl.figure()

phiso = np.linspace(0.00, 0.99, num=500)
pl.plot(phiso, (part_coeff(phiso, -0.90)),'b--', label=' delta f = -0.90')
pl.plot(phiso, (part_coeff(phiso, -0.30)),'g--', label=' delta f = -0.30')
pl.plot(phiso, (part_coeff(phiso, 0.0)),'r--', label=' delta f = 0.0')
#pl.plot(phiso, (part_coeff(phiso, (1 - phiso)/(0.01 + phiso))),'r--', label=' delta f = p')
pl.plot(phiso, (part_coeff(phiso, phiso)),'r--', label=' delta f = phiso')
pl.plot(phiso, (part_coeff(phiso, 0.30)),'b--,', label=' delta f = 0.30')
pl.plot(phiso, (part_coeff(phiso, 0.90)),'g--,', label=' delta f = 0.90')
#pl.plot(phiso, (part_coeffXS(1. - phiso, 0.90)),'g--,', label=' delta f = 0.90')


pl.legend(loc='lower right')
pl.xlabel('Number Fraction Polymers Outside Pore')
pl.ylabel('Partition Coefficent ($/phi/$in/$/phi/$out)')
pl.title('Partition Coefficent vs Number Fraction Outside Pore')
pl.axis([0., 1.0, 0., 1.25*10**(26)])

#pl.savefig('Partition_Coefficent.png')
pl.show()

