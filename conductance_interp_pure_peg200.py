import numpy as np
from numpy import exp
import matplotlib.pyplot as pl
from scipy import interpolate
from pylab import *
from scipy.signal import cspline1d

p200 = [0,2,5,10,15,25]
c200 = [0.939,0.864,0.751,0.590,0.392,0.295]
f = interpolate.UnivariateSpline(p200,c200,k=2)
#f3 = interpolate.UnivariateSpline(p200,c200,k=1)
f3 = interpolate.interp1d(p200,c200)
#fcs = cspline1d(p200,c200)
p200yerr = [0.019,0.015,0.018,0.013,0.012,0.011]
errorbar(p200,c200,p200yerr)

p200_add = [3,6,9,12,15]
c200_add = [0.792,0.666,0.527,0.431,0.324]
g = interpolate.UnivariateSpline(p200_add,c200_add,k=2)
#g3 = interpolate.UnivariateSpline(p200_add,c200_add,k=1)
g3 = interpolate.interp1d(p200_add,c200_add)
p200_addyerr = [0.009,0.015,0.021,0.014,0.004]
errorbar(p200_add,c200_add,p200_addyerr)

#p1k = [0,15,30]
#c1k = [0.939,0.413,0.224]
#h = interpolate.UnivariateSpline(p1k,c1k)
#h3 = interpolate.interp1d(p1k,c1k)

p1k_add = [3,6,9,12,15]
c1k_add = [0.723,0.563,0.504,0.396,0.281]
j = interpolate.UnivariateSpline(p1k_add,c1k_add,k=2)
#j3 = interpolate.UnivariateSpline(p1k_add,c1k_add,k=1)
j3 = interpolate.interp1d(p1k_add,c1k_add)
p1k_addyerr = [0.043,0.045,0.092,0.066,0.026]
errorbar(p1k_add,c1k_add,p1k_addyerr)

p3400 = [0,2,5,10,15,20,25]
c3400 = [0.93,0.930,0.954,0.993,0.936,0.843,0.526]
k = interpolate.UnivariateSpline(p3400,c3400,k=2)
#k3 = interpolate.UnivariateSpline(p3400,c3400,k=1)
k3 = interpolate.interp1d(p3400,c3400)
p3400_yerr = [0.019,0.028,0.020,0.020,0.011,0.086,0.190]
errorbar(p3400,c3400,p3400_yerr)

p10k = [0,7.5,15,20,25,30]
c10k = [0.939,0.980,1.011,1.071,0.980,0.522]
m = interpolate.UnivariateSpline(p10k,c10k,k=2)
#m3 = interpolate.UnivariateSpline(p10k,c10k,k=1)
m3 = interpolate.interp1d(p10k,c10k)
p10k_yerr = [0.019,0.008,0.021,0.046,0.046,0.106]
errorbar(p10k,c10k,p10k_yerr)

#f = interpolate.BarycentricInterpolator(x,y)
#f = interpolate.interp1d(x,y)

perc200 = np.arange(3,25,1)
cond200 = f(perc200)
cond200_3 = f3(perc200)

perc200_add = np.arange(3,15,1)
cond200_add = g(perc200_add)
cond200_3_add = g3(perc200_add)

#perc1k = np.arange(3,30,1)
#cond1k = h(perc1k)
#cond1k_3 = h3(perc1k)

perc1k_add = np.arange(3,15,1)
cond1k_add = j(perc1k_add)
cond1k_3_add = j3(perc1k_add)


perc3400 = np.arange(3,25,1)
cond3400 = k(perc3400)
cond3400_3 = k3(perc3400)

perc10k = np.arange(3,30,1)
cond10k = m(perc10k)
cond10k_3 = m3(perc10k)
	
perc_expt = [3,6,9,12,15]

print 'Conductance for PEG200 percentages:       3=%.3f,%.3f   6=%.3f,%.3f   9=%.3f,%.3f   12=%.3f,%.3f   15=%.3f,%.3f' % (f(3),f3(3),f(6),f3(6),f(9),f3(9),f(12),f3(12),f(15),f3(15))

print 'Conductance for PEG200_added percentages: 3=%.3f,%.3f   6=%.3f,%.3f   9=%.3f,%.3f   12=%.3f,%.3f   15=%.3f,%.3f' % (g(3),g3(3),g(6),g3(6),g(9),g3(9),g(12),g3(12),g(15),g3(15))



#print 'Conductance for PEG1k  percentages: 3=%.3f: 6=%.3f: 9=%.3f: 12=%.3f: 15=%.3f' % (h(3),h(6

print 'Conductance for PEG1k_added percentages:  3=%.3f,%.3f   6=%.3f,%.3f   9=%.3f,%.3f   12=%.3f,%.3f   15=%.3f,%.3f' % (j(3),j3(3),j(6),j3(6),j(9),j3(9),j(12),j3(12),j(15),j3(15))



print 'Conductance for PEG3400 percentages:      3=%.3f,%.3f   6=%.3f,%.3f   9=%.3f,%.3f   12=%.3f,%.3f   15=%.3f,%.3f' % (k(3),k3(3),k(6),k3(6),k(9),k3(9),k(12),k3(12),k(15),k3(15))



print 'Conductance for PEG10k percentages:       3=%.3f,%.3f   6=%.3f,%.3f   9=%.3f,%.3f   12=%.3f,%.3f   15=%.3f,%.3f' % (m(3),m3(3),m(6),m3(6),m(9),m3(9),m(12),m3(12),m(15),m3(15))



pl.plot(perc_expt,f(perc_expt),marker='o', markerfacecolor= 'k',linestyle = 'none', markersize= 6,color = 'k')
pl.plot(perc_expt,g(perc_expt),marker='D', markerfacecolor= 'k',linestyle = 'none', markersize= 6,color = 'k')
#pl.plot(perc_expt,h(perc_expt),marker='+', markerfacecolor= 'k',linestyle = 'none', markersize= 6,color ='k')
pl.plot(perc_expt,j(perc_expt),marker='s', markerfacecolor= 'k',linestyle = 'none', markersize= 6,color = 'k')
pl.plot(perc_expt,k(perc_expt),marker='^', markerfacecolor= 'k',linestyle = 'none', markersize= 6,color = 'k')
pl.plot(perc_expt,m(perc_expt),marker='x', markerfacecolor= 'k',linestyle = 'none', markersize= 6,color = 'k')


#pl.plot(perc_expt,f3(perc_expt),marker='o', markerfacecolor= 'w',linestyle = ':', markersize= 6,color = 'k')
pl.plot(perc_expt,g3(perc_expt),marker='D', markerfacecolor= 'w',linestyle = ':', markersize= 6,color = 'k')
#pl.plot(perc_expt,3h(perc_expt),marker='+', markerfacecolor= w',linestyle = ':', markersize= 6,color ='k')
pl.plot(perc_expt,j3(perc_expt),marker='s', markerfacecolor= 'w',linestyle = ':', markersize= 6,color = 'k')
pl.plot(perc_expt,k3(perc_expt),marker='^', markerfacecolor= 'w',linestyle = ':', markersize= 6,color = 'k')
pl.plot(perc_expt,m3(perc_expt),marker='x', markerfacecolor= 'w',linestyle = ':', markersize= 6,color = 'k')

#pl.plot(6,f(6),marker='^',linestyle = 'none', markersize= 6,color = 'r', label=r'$6\ percent:\ %.3f\ nS$' % f(6))
#pl.plot(9,f(9),marker='^',linestyle = 'none', markersize= 6,color = 'b', label=r'$9\ percent:\ %.3f\ nS$' % f(9))
#pl.plot(12,f(12),marker='^',linestyle = 'none', markersize= 6,color = 'c',label=r'$12\ percent:\ %.3f\ nS$' % f(12))
#pl.plot(15,f(15),marker='^',linestyle = 'none', markersize= 6,color = 'm', label=r'$15\ percent:\ %.3f\ nS$' % f(15))
#

#pl.plot(p200,c200,'o',linestyle=':')
pl.plot(perc200,cond200,'-.',linewidth=2.5, color='b')
#pl.plot(p200_add,c200_add,'D',
pl.plot(perc200_add,cond200_add,'-.',linewidth=2.5,color='g')
#pl.plot(p1k,c1k,'+',perc1k,cond1k,'-',color='g')
#pl.plot(p1k_add,c1k_add,'s',
pl.plot(perc1k_add,cond1k_add,'-.',linewidth=2.5,color='r')
#pl.plot(p3400,c3400,'^',
pl.plot(perc3400,cond3400,'-.',linewidth=2.5,color='c')
#pl.plot(p10k,c10k,'x',
pl.plot(perc10k,cond10k,'-.',linewidth=2.5,color='m')

#pl.plot(perc200,cond200_3,':', color='b')
#pl.plot(p200_add,c200_add,'D',perc200_add,cond200_3_add,':',color='r')
##pl.plot(p1k,c1k,'+',perc1k,cond1k,'-',color='g')
#pl.plot(p1k_add,c1k_add,'s',perc1k_add,cond1k_3_add,':',color='c')
#pl.plot(perc3400,cond3400_3,linestyle = ':',color='y')
#pl.plot(p10k,c10k,'x',perc10k,cond10k_3,':',color='m')

pl.plot(-2,-2,'-',color='b', label='PEG200')
pl.plot(-2,-2,'-',color='g', label='PEG200 Added')
#pl.plot(-2,-2,'-',color='g', label='PEG1000')
pl.plot(-2,-2,'-',color='r', label='PEG1000 Added')
pl.plot(-2,-2,'-',color='c', label='PEG3400')
pl.plot(-2,-2,'-',color='m', label='PEG10000')

pl.axis([-2.0,31.0,0.05,1.2])
pl.legend(loc='lower left')
pl.title('Conductance vs. % Weight PEG')
pl.xlabel('Percent PEG')
pl.ylabel('Conductance nS')

pl.savefig('1206121_conductance_vs_peg.eps', dpi=600)
pl.show()
