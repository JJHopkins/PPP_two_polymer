#!/usr/bin/python

import numpy as np
from numpy import log,exp
import scipy.optimize as opt
import matplotlib.pyplot as pl
import matplotlib
alpha = 0.49
Ns = 9
Nb = 79
df = 0
phib = 0.10
phi = .10

def f(p, phi, phib, df):
	""" Implicit definition of P(phi, dF) """
	return - log(p) + exp(- df + Ns*(log((1 - p*phi)/(1 - phi - phib)) + \
		(p - 1)*phi - phib + (9./4)*alpha*((phi + phib)**(5./4) - (p*phi)**(5./4))))

#Optimization algorithms, if opposite sign use brent, else use newton
def P(phi, phib, df):
	""" Numerically solve for partition coefficient as a
	    function of \phi_s """
	if f(0,phi,phib,df)*f(1,phi,phib,df) < 0:
		print 'bisect'
		return opt.bisect(f, 0, 1, args=(phi,phib,df)) # Bisection method
	else:
		return opt.newton(f, 1.0, args=(phi,phib,df)) # Newton-Raphson
	print P

