from readfiles import readlookup

import matplotlib.pyplot as plt
from matplotlib import pyplot as plt

filename='bindingengdata.dat'

colnum=6

#data = readlookup(filename,colnum)

data=[]

f = open (filename) 
for line in f:
    char = line.split()
    for i in char:
       data.append(i)
f.close()

for i in range(0,len(data)):
   data[i]=float(data[i])

He= range(len(data)/6)
V= range(len(data)/6)
HeVrat= range(len(data)/6)
EbSIA= range(len(data)/6)
EbHe= range(len(data)/6)
EbV= range(len(data)/6)

He=data[0::6]
V=data[1::6]
HeVrat=data[2::6]
EbSIA=data[3::6]
EbHe=data[4::6]
EbV=data[5::6]


#plt.plot([0.1,1.2E+07/sim], [48.8,48.8], 'k-',lw=3)


xarrn=0
yarrn=1
#plt.plot(results1[xarrn], results1[yarrn])
#plt.plot(results2[xarrn], results2[yarrn])


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
ax.plot(HeVrat,EbHe,'.')#,label= 'Max void size = 1E6')
#ax.plot(results2[xarrn], results2[yarrn],label= 'Max void size = 1E5')
#ax.plot(results3[xarrn], results3[yarrn],label= 'Max void size = 1E4')
#ax.plot(results4[xarrn], results4[yarrn],label= 'Max void size = 1E3')
#ax.set_xlim([0.1,10])
#ax.set_ylim([0.0,5.5])
#ax.set_ybound(lower=round(lowlim),upper=round(uplim))
#ax.set_yscale('log')
#ax.set_xscale('log',basex=10)

legend = ax.legend(loc='lower right', shadow=True)
#plt.legend( labels,
#        loc = 'upper left',handlelength=1.5, handleheight=1,fontsize = 'small')

#plt.ylim(0, 0.05E-4)
plt.tight_layout()
plt.show()
