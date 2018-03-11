import numpy as np
from numpy import linspace, exp, sin, pi, sqrt, log
from mpl_toolkits.mplot3d import Axes3D
from readfiles import readlookup
import matplotlib.pyplot as plt
from matplotlib import pyplot as plt
import matplotlib.ticker as mtick
import pylab as p


filename='../Output/Clusters.dat'
colnum=8
results1 = readlookup(filename,colnum)

#filename='../Output/Clusters09em.dat'
#results2 = readlookup(filename,colnum)

#filename='../Output/Clusters067emnocorr.dat'
#results3 = readlookup(filename,colnum)

#filename='../Output/Clusters09emnocorr.dat'
#results4 = readlookup(filename,colnum)


Fedens=8.493592929137993e+22

experimentaldatay=np.array([1.8E+23,1.0E+24,1.9E+24,5E+24])
experimentaldatax=np.array([8.E-05,9.E-04,9.E-03,1.E-01])


upper_error = np.array([4E+23,2.E+24,4.E+24,8.E+24])
lower_error =  np.array([7.E+22,6.E+23,8.E+23,2.4E+24])
asymmetric_error = [lower_error, upper_error]
#plt.plot([0.1,1.2E+07/sim], [48.8,48.8], 'k-',lw=3)


xarrn=0
yarrn=6
#plt.plot(results1[xarrn], results1[yarrn])
#plt.plot(results2[xarrn], results2[yarrn])


#ax3.set_title('different sample sizes')

#labels=['Without Grouping Method'
#,'With Grouping Method'
#]
lowlim=min(results1[yarrn])
uplim=max(results1[yarrn])
#plt.xscale('log')
#plt.yscale('log')
#plt.ylim([9673000300000000,9.5784429999999999e+22])
#results1[yarrn]=[log(x) for x in results1[yarrn] ]
fig, ax = plt.subplots()
ax.plot(results1[xarrn], results1[yarrn],label= 'Low He implantation rate')
#ax.plot(results2[xarrn], results2[yarrn],label= 'EM=0.9')
#ax.plot(results1[xarrn], results3[yarrn],label= 'Em=0.67, no sinkstrengthcorrection')
#ax.plot(results2[xarrn], results4[yarrn],label= 'EM=0.9, no sinkstrengthcorrection')
#ax.plot(results3[xarrn], results3[yarrn],label= 'Max void size = 1E4')
#ax.plot(results4[xarrn], results4[yarrn],label= 'Max void size = 1E3')

#ax.errorbar(experimentaldatax, experimentaldatay, yerr=asymmetric_error, fmt='o',color='r',linewidth=2)
#ax.plot(experimentaldatax, experimentaldatay,linewidth=20, color='red',linestyle='o')

plt.rc('font', family='serif')
ax.set_xlabel(r'Dose (dpa)')
ax.set_ylabel(r'Total Cluster Density (m$^{-3}$)')

#plt.ylim(2.76E+21 , 1.E25)
#ax.set_ylim(round(lowlim),round(uplim))
#ax.set_ybound(lower=round(lowlim),upper=round(uplim))
ax.set_yscale('log')
ax.set_xscale('log')


legend = ax.legend(loc='lower right', shadow=True)
#plt.legend( labels,
#        loc = 'upper left',handlelength=1.5, handleheight=1,fontsize = 'small')

#plt.ylim(0, 0.05E-4)
plt.tight_layout()
plt.show()
