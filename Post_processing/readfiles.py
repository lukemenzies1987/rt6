def readfile3(filename,ncolumns1,ncolumns2,ncolumns3):

  from numpy import zeros

  #filename='TungstenInput_det0.m'
  lines =[]
  data=[]
  f = open (filename) 
  for line in f:
    lines.append(line)
    char = line.split()
    for i in char:      
       data.append(i)
       #print i
  f.close()

  end=data.index('];')
  array1 =data[:end]
  data.remove('];')
  srt=end
  end=data.index('];')
  array2 =data[srt:end]
  data.remove('];')
  srt=end
  end=data.index('];')
  array3 =data[srt:end]

  #array1= array1[1:]
  #array2= array2[1:]
  #array3= array3[1:]

  array1 = [float(i) for i in array1[3:]]
  array2 = [float(i) for i in array2[3:]]
  array3 = [float(i) for i in array3[3:]]

  name1=zeros((ncolumns1,len(array1)/ncolumns1))
  for j in range(ncolumns1):
    for i in range(len(array1)/ncolumns1):
      name1[j][i]=array1[j+i*ncolumns1]  

  name2=zeros((ncolumns2,len(array2)/ncolumns2))
  for j in range(ncolumns2):
    for i in range(len(array2)/ncolumns2):
      name2[j][i]=array2[j+i*ncolumns2]

  name3=zeros((ncolumns3,len(array3)/ncolumns3))
  for j in range(ncolumns3):
    for i in range(len(array3)/ncolumns3):
      name3[j][i]=array3[j+i*ncolumns3]  

  return name1,name2,name3

def readfile2(filename,ncolumns1,ncolumns2):

  from numpy import zeros

  #filename='TungstenInput_det0.m'
  lines =[]
  data=[]
  f = open (filename) 
  for line in f:
    lines.append(line)
    char = line.split()
    for i in char:      
       data.append(i)
       #print i
  f.close()

  end=data.index('];')
  array1 =data[:end]
  data.remove('];')
  srt=end
  end=data.index('];')
  array2 =data[srt:end]
  data.remove('];')


  #array1= array1[1:]
  #array2= array2[1:]
  #array3= array3[1:]

  array1 = [float(i) for i in array1[3:]]
  array2 = [float(i) for i in array2[3:]]


  name1=zeros((ncolumns1,len(array1)/ncolumns1))
  for j in range(ncolumns1):
    for i in range(len(array1)/ncolumns1):
      name1[j][i]=array1[j+i*ncolumns1]  

  name2=zeros((ncolumns2,len(array2)/ncolumns2))
  for j in range(ncolumns2):
    for i in range(len(array2)/ncolumns2):
      name2[j][i]=array2[j+i*ncolumns2]


  return name1,name2
 
def readlookup(filename,colnum):
  from numpy import zeros

  def get_digits(str1):
      c = ""
      for i in str1:
          if i.isdigit() or i ==" " or i == "+" or i=="." or i=="E" or i=="-" or i==";":
            if i==";":
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
         #print i
  f.close()
  length=len(bufdata)/colnum
  data=zeros((colnum,length))

  for i in range(length):
    for j in range(colnum):
      data[j][i]= (float(bufdata[colnum*i+j]))
  return data
