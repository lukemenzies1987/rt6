import numpy as np
from numpy import linspace, exp, sin, pi, sqrt, log
from mpl_toolkits.mplot3d import Axes3D
from readfiles import readlookup
import matplotlib.pyplot as plt
from matplotlib import pyplot as plt
import matplotlib.ticker as mtick
import pylab as p
from matplotlib import rc

#filename='./finalreportvoid1/Clusters.dat'
filename='../Output/Clusters.dat'
#linestart=3
colnum=8
results1 = readlookup(filename,colnum)


#filename='./finalreportvoid2/Clusters.dat'
#filename='./Clustersmid.dat'
#results2 = readlookup(filename,colnum)

#filename='./finalreportvoid3/Clusters.dat'
#filename='./Clustershigh.dat'
#results3 = readlookup(filename,colnum)


Fedens=8.493592929137993e+22


#plt.plot([0.1,1.2E+07/sim], [48.8,48.8], 'k-',lw=3)
rc('text', usetex=True)

xarrn=0
yarrn=6
#plt.plot(results1[xarrn], results1[yarrn])
#plt.plot(results2[xarrn], results2[yarrn])

#fig, ax1 = plt.subplots()
#t = np.arange(0.01, 10.0, 0.01)
#s1 = np.exp(t)
#ax1.plot(t, s1, 'b-')
#ax1.set_xlabel('time (s)')
# Make the y-axis label and tick labels match the line color.
#ax1.set_ylabel('exp', color='b')
#for tl in ax1.get_yticklabels():
#    tl.set_color('b')


#ax2 = ax1.twinx()
#s2 = np.sin(2*np.pi*t)
#ax2.plot(t, s2, 'r.')
#ax2.set_ylabel('sin', color='r')

experimentaldatay=np.array([3.0E+21,1.7E+22])
experimentaldatax=np.array([1.0E+01,1.0E+02])
#ax3.set_title('different sample sizes')

#labels=['Without Grouping Method'
#,'With Grouping Method'
#]
#lowlim=min(results1[yarrn])
#uplim=max(results1[yarrn])

# example error bar values that vary with x-position

# error bar values w/ different -/+ errors

upper_error = np.array([3.0E+21,1.5E+22])
lower_error =  np.array([1.5E+21,1.0E+22])
asymmetric_error = [lower_error, upper_error]
#plt.xscale('log')
#plt.yscale('log')
#plt.ylim([9673000300000000,9.5784429999999999e+22])
#results1[yarrn]=[log(x) for x in results1[yarrn] ]
fig, ax = plt.subplots()
ax.plot(results1[yarrn+1], results1[yarrn],linewidth=2, color='green',linestyle='--',label= 'Low')
#lines2=ax.plot(results2[yarrn+1], results2[yarrn],linewidth=2, color='blue',linestyle='--',label= 'Medium')
#lines3=ax.plot(results3[yarrn+1], results3[yarrn],linewidth=2, color='red',linestyle='--',label= 'High')

#ax.plot(experimentaldatax, experimentaldatay,linewidth=20, color='red',linestyle='o')

ax.text(3.E-3, 7.E21, r'He Implantation T=623K')
#ax.set_ylim(round(lowlim),round(uplim))
#ax.set_ybound(lower=round(lowlim),upper=round(uplim))
ax.set_yscale('log')
ax.set_xscale('log')



plt.rc('font', family='serif')
#ax.set_xlabel(r'Dose (dpa)')
ax.set_xlabel(r'He Concentration (appm)')
ax.set_ylabel(r'Total Cluster Density (m$^{-3}$)')
ax2 = ax.twiny()

ax2.errorbar(experimentaldatax, experimentaldatay, yerr=asymmetric_error, fmt='o',color='r',linewidth=2)
lines1=ax2.plot(results1[xarrn], results1[yarrn],linewidth=2, color='green',linestyle='--',label= 'Low')
#ax2.plot(results2[xarrn], results2[yarrn],linewidth=2, color='blue',linestyle='--',label= 'Medium')
#ax2.plot(results3[xarrn], results3[yarrn],linewidth=2, color='red',linestyle='--',label= 'High')


ax2.lines.remove(lines1[0])
#ax.lines.remove(lines2[0])
#ax.lines.remove(lines3[0])
ax2.set_xscale('log')
ax2.set_xlabel(r'Dose (dpa)')
legend = ax2.legend(loc='upper right', shadow=True)
#plt.legend( labels,
#        loc = 'upper left',handlelength=1.5, handleheight=1,fontsize = 'small')

#plt.ylim(2.76E+20 , 1.E22)
plt.tight_layout()
plt.show()
