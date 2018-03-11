"""
===========================
Plots with different scales
===========================

Demonstrate how to do two plots on the same axes with different left and
right scales.

The trick is to use *two different axes* that share the same *x* axis.
You can use separate `matplotlib.ticker` formatters and locators as
desired since the two axes are independent.

Such axes are generated by calling the `Axes.twinx` method.  Likewise,
`Axes.twiny` is available to generate axes that share a *y* axis but
have different top and bottom scales.

The twinx and twiny methods are also exposed as pyplot functions.

"""

import numpy as np
import matplotlib.pyplot as plt
from readfiles import readlookup

xarrn=0
yarrn=6
colnum=8
filename='../Output/Clusters.dat'
results1 = readlookup(filename,colnum)
Gdpa=1.84E-6
length=len(results1[0])
for i in range(length): 
  results1[0][i]=results1[0][i]*Gdpa

"""
filename='./Output1/ClustersMed.dat'
results2 = readlookup(filename,colnum)
Gdpa=5.52E-7
length=len(results2[0])
for i in range(length): 
  results2[0][i]=results2[0][i]*Gdpa


filename='./Output1/ClustersHigh.dat'
results3 = readlookup(filename,colnum)
Gdpa=1.84E-6
length=len(results3[0])
for i in range(length): 
  results3[0][i]=results3[0][i]*Gdpa

"""
fig, ax1 = plt.subplots()


ax1.plot(results1[yarrn+1], results1[xarrn],linewidth=2, color='green',linestyle='--',label= 'Low')
#ax1.plot(results2[yarrn+1], results2[xarrn],linewidth=2, color='blue',linestyle='--',label= 'Medium')
#ax1.plot(results3[yarrn+1], results3[xarrn],linewidth=2, color='red',linestyle='--',label= 'High')
ax1.cla()

# Make the y-axis label, ticks and tick labels match the line color.


ax2 = ax1.twiny()

ax2.plot(results1[xarrn], results1[yarrn],linewidth=2, color='green',linestyle='--',label= 'Low')
#ax2.plot(results2[xarrn], results2[yarrn],linewidth=2, color='blue',linestyle='--',label= 'Medium')
#ax2.plot(results3[xarrn], results3[yarrn],linewidth=2, color='red',linestyle='--',label= 'High')


ax1.text(2.E-3, 3.E22, r'He Implantation T=623K')
ax1.set_xlabel(r'He Concentration (appm)')
ax1.set_ylabel(r'Total Cluster Density (m$^{-3}$)')
ax2.set_xlabel(r'Dose (dpa)')


if yarrn==1:
  ax1.set_ylabel(r'Conc. of Vacancies (per atom)')
if yarrn==2:
  ax1.set_ylabel(r'Conc. of SIA (per atom)')
if yarrn==3:
  ax1.set_ylabel(r'Conc. of He (per atom)')
if yarrn==4:
  ax1.set_ylabel(r'Conc. of He at GB (appm)')
  plt.axhline(y=48.0, xmin=results1[xarrn][0], xmax=results1[xarrn][-1], linewidth=2, color = 'k')
  ax.annotate('Critical Concentration', xy=(0.0005,48.8), xytext=(0.0005, 2.0),
              arrowprops=dict(facecolor='black', shrink=0.0005),)
if yarrn==6:
  ax1.set_ylabel(r'Total Cluster Density (m$^{-3}$)')
  experimentaldatay=np.array([4.0E+21,1.7E+22])
  experimentaldatax=np.array([1.0E-03,1.0E-02])
  upper_error = np.array([3.0E+21,1.5E+22])
  lower_error =  np.array([2.5E+21,11.0E+21])
  asymmetric_error = [lower_error, upper_error]
  ax2.errorbar(experimentaldatax, experimentaldatay, yerr=asymmetric_error, fmt='o',color='r',linewidth=2)
  plt.xlim(8.E-7 , 3.E-2)
  plt.ylim(7.2E+19 , 4.E22)
legend = ax2.legend(loc='lower right', shadow=True)

ax1.set_yscale('log')
ax2.set_yscale('log')
ax2.set_xscale('log')
ax1.set_xscale('log')



fig.tight_layout()
plt.show()