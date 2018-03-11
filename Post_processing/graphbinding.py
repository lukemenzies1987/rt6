from readfiles import readlookup

import matplotlib.pyplot as plt
from matplotlib import pyplot as plt
from numpy import zeros
filename='../Output/BindingEnergyV-Cluster.dat'

#colnum=6

data=[]

f = open (filename) 
for line in f:
    char = line.split()
    for i in char:
       data.append(i)
f.close()

for i in range(0,len(data)):
   data[i]=float(data[i])

width=100
length=11
#data = readlookup(filename,colnum)

bdata=zeros((length,width))
#plt.plot([0.1,1.2E+07/sim], [48.8,48.8], 'k-',lw=3)
k=0
for i in range(length):
  for j in range(width):
    bdata[i][j]=data[k]
    k+=1
xarrn=0
yarrn=1

vacnum= [i+1 for i in xrange(width)]

#ax3.set_title('different sample sizes')

#labels=['Without Grouping Method'
#,'With Grouping Method'
#]
#lowlim=min(data[yarrn])
#uplim=max(data[yarrn])
#plt.xscale('log')
#plt.yscale('log')
#plt.xlim([0.1,10])
#results1[yarrn]=[log(x) for x in results1[yarrn] ]
fig, ax = plt.subplots()
for i in range(5):
  ax.plot(vacnum,bdata[i],'-',label= str(i)+' He')

#ax.plot(results4[xarrn], results4[yarrn],label= 'Max void size = 1E3')
#ax.set_xlim([0.1,10])
#ax.set_ylim([0.0,5.5])
#ax.set_ybound(lower=round(lowlim),upper=round(uplim))
#ax.set_yscale('log')
#ax.set_xscale('log',basex=10)
ax.text(25, 2.3, 'T=423K')
ax.set_ylabel("Binding Energy (eV)")
ax.set_xlabel("Number of vacancies in cluster")
legend = ax.legend(loc='lower right', shadow=True)
#plt.legend( labels,
#        loc = 'upper left',handlelength=1.5, handleheight=1,fontsize = 'small')

#plt.ylim(0, 0.05E-4)
plt.tight_layout()
plt.show()
