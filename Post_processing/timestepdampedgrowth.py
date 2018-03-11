import matplotlib.pyplot as plt
from numpy import linspace, exp, sin, pi, sqrt


def gwithdamp(arg,dtmax,xmax,dtconvrt):
  from numpy import exp
  y=dtmax*exp(-xmax*dtconvrt/arg)
  return y

dtmax=10.
xmax=1000.
npts=1000
dtconvrt=0.01

x= linspace(0,xmax,npts)
plt.plot(x,gwithdamp(x,dtmax,xmax,dtconvrt))
plt.xlabel('time (s)')
plt.ylabel('dt (s)')
plt.show()
