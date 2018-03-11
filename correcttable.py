import numpy as np
from numpy import linspace, exp, sin, pi, sqrt, log
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import pyplot as plt
import matplotlib.ticker as mtick
import pylab as p

#This python script formats csv files to put into the correct format for tabled binding energy values. 

def readlookup(filename,colnum):
  from numpy import zeros

  def get_digits(str1):
      c = ""
      for i in str1:
          if i.isdigit() or i ==" " or i == "+" or i=="." or i=="E" or i=="-" or i==";" or i=="\t":
            if i==";":
              i=" "
            if i=="\t":
              i=" "
            c += i        
      return c
  #colnum=3
  lines =[]
  bufdata=[]
  data=[]
  #f = open ('W182 Total.csv') 
  f = open (filename) 
  for line in f:
      if line.find('inf')>0:
         continue
      if line.find('nan')>0:
         continue
      lines.append(line)
      line2= get_digits(line)
      char = line2.split()
      for i in char: 
         if len(i)==1:
            continue      
         bufdata.append(i)

  f.close()
  length=len(bufdata)/colnum
  data=zeros((colnum,length))

  for i in range(length):
    for j in range(colnum):
      data[j][i]= (float(bufdata[colnum*i+j]))
  return data


coln=6
filename='bindingengdata.dat'

table=readlookup(filename,coln)
col1 = ["%.4f" % member for member in table[0]]
col2 = ["%.4f" % member for member in table[1]]
col3 = ["%.4f" % member for member in table[2]]
col4 = ["%.4f" % member for member in table[3]]
col5 = ["%.4f" % member for member in table[4]]
col6 = ["%.4f" % member for member in table[5]]

noa=len(table[0])

f = open('bindingengdata2.dat','w')
f.write(str(noa))
f.write('\n')
f.write('\n')
for i in range(noa):
   f.write(str(col1[i]))
   f.write('  ')
   f.write(str(col2[i]))
   f.write('  ')
   f.write(str(col3[i]))
   f.write('  ')
   f.write(str(col4[i]))
   f.write('  ')
   f.write(str(col5[i]))
   f.write('  ')
   f.write(str(col6[i]))
   f.write('  ')
   f.write('\n')

f.close()
