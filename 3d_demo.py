import pylab as p
from mpl_toolkits.mplot3d import Axes3D 
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as pl
from numpy import *


xi = arange(-5, 5, 0.25)
yi = arange(-5, 5, 0.25)
#xi = linspace(min(X), max(X))
#yi = linspace(min(Y), max(Y))
#zi = sin(-xi*yi)
zi = arange(-5, 5, 0.25)

fig = pl.figure()
ax = Axes3D(fig)
xim,yim = meshgrid(xi,yi)
ax.plot_surface(xim,yim,zi)

pl.show()
#
#ax = fig.gca(projection='3d')
#X = arange(-5, 5, 0.25)
#Y = arange(-5, 5, 0.25)
#X, Y = meshgrid(X, Y)
#R = sqrt(X**2 + Y**2)
#Z = sin(R)
#surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm,
#        linewidth=0, antialiased=False)
##ax.set_zlim(-1.01, 1.01)
#
##ax.zaxis.set_major_locator(LinearLocator(10))
##ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
#
#fig.colorbar(surf, shrink=0.5, aspect=5)
#
#plt.show()
#
