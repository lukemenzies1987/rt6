import numpy as np
from numpy import linspace, exp, sin, pi, sqrt, log
from mpl_toolkits.mplot3d import Axes3D
from readfiles import readlookup
import matplotlib.pyplot as plt
from matplotlib import pyplot as plt
import matplotlib.ticker as mtick
import pylab as p



colnum=2
filename='./sdfforvoids.dat'
results1 = readlookup(filename,colnum)

filename='./sdfforclusters.dat'
#filename='./sdfforvoidsanogrouping.dat'
results2 = readlookup(filename,colnum)

Fedens=8.493592929137993e+22


#plt.plot([0.1,1.2E+07/sim], [48.8,48.8], 'k-',lw=3)


xarrn=0
yarrn=1
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
ax.plot(results1[xarrn], results1[yarrn],'.')
ax.plot(results2[xarrn], results2[yarrn])
#ax.set_ylim(1E16,round(uplim))
#ax.set_ybound(lower=round(lowlim),upper=round(uplim))
ax.set_yscale('log')
#ax.set_xscale('log')

#legend = ax.legend(loc='lower right', shadow=True)
#plt.legend( labels,
#        loc = 'upper left',handlelength=1.5, handleheight=1,fontsize = 'small')

#plt.ylim(0, 0.05E-4)
plt.tight_layout()
plt.show()
