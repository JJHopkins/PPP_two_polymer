#!/usr/bin/python

# Import commands
from pylab import *
from scipy.signal import cspline1d, cspline1d_eval

# Coarse values
x_coarse = [0,2,5,10,15,25]
y_coarse = [0.939,0.864,0.751,0.590,0.392,0.295]

# Fine values to interpolate to:
#x_fine = linspace(0,2*pi,1e6)
x_fine = linspace(0,26,1)

# Linear interpolation to full range
y_lininterp=interp(x_fine,x_coarse,y_coarse)

# --- Cubic spline interpolation to full range
# fine coefficients:
cspline_coeffs = cspline1d(y_coarse)
y_cspline = cspline1d_eval(cspline_coeffs,x_fine,
                           dx=x_coarse[1]-x_coarse[0],
                           x0=x_coarse[0])


figure()
plot(x_coarse,y_coarse,'ko',label='y=sin(x)')
plot(x_fine,y_lininterp,'r',label='linear interpolation')
plot(x_fine,y_cspline,'k',label='cubic spline fit')
legend()
show()
