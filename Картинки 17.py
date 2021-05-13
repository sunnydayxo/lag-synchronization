from numpy import *
import numpy as np
from matplotlib import *
from scipy import *
from pylab import figure, show, setp
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import scipy.io as sio
from scipy.signal import find_peaks

def step_rossler(x1,y1,z1,x2,y2,z2, dt, a, b, c, K,w1,w2):
    dx1 = dt*(-w1*y1-z1 + K*(x2-x1))
    dy1 = dt*(w1*x1 + a*y1)
    dz1 = dt*(b + z1*(x1 - c))

    dx2 = dt*(-w2*y2-z2 + K*(x1-x2))
    dy2 = dt*(w2*x2 + a*y2)
    dz2 = dt*(b + z2*(x2 - c))
    
    x1,y1,z1 = x1+dx1, y1+dy1, z1+dz1
    x2,y2,z2 = x2+dx2, y2+dy2, z2+dz2
    
    return (x1,y1,z1, x2,y2,z2)
def orbit_rossler(x1,y1,z1,x2,y2,z2, dt, a, b, c, K, w1,w2, numsteps):
    XYZ = np.zeros((numsteps+1)*6)
    XYZ = np.reshape(XYZ, (numsteps+1, 6))
    XYZ[0,:] = x1,y1,z1,x2,y2,z2
    for i in range(1, numsteps+1):
        x1,y1,z1,x2,y2,z2 = step_rossler(x1,y1,z1,x2,y2,z2, dt, a,b,c,K, w1, w2)

        XYZ[i,:] = x1,y1,z1,x2,y2,z2
    return XYZ
def simFn(x1,x2, skew):
    if skew == 0:
        diff_skew = x1 - x2
    else:
        diff_skew = x1[skew:] - x2[:-skew]
    
    diff_skew_avg = np.average(diff_skew*diff_skew)
    
    x1_sq_avg = np.average(x1*x1)
    x2_sq_avg = np.average(x2*x2)
    factor = np.sqrt(x1_sq_avg*x2_sq_avg)
    
    return diff_skew_avg/factor
x1_0, y1_0, z1_0 = [0.001]*3
x2_0, y2_0, z2_0 = [0.002]*3

a=0.165
b=0.2
c=10
w1=0.99
w2=0.95

dt = 0.01
t = np.arange(0,1000, dt)
numsteps = len(t) - 1
K = 0.2
XYZ_0_2 = orbit_rossler(x1_0, y1_0, z1_0, x2_0, y2_0, z2_0,dt,a,b,c, K, w1,w2, numsteps)
x1,y1,z1,x2,y2,z2 = [ XYZ_0_2[:,i] for i in range(6) ]
x1 = XYZ_0_2[:,0]
x2 = XYZ_0_2[:,3]

x1 = x1[10000:]
x2 = x2[10000:]


tau = np.arange(0,30,dt)
S = np.array([ simFn(x2,x1,int(_tau/dt)) for _tau in tau ])
minskew = np.argmin(S[:1000])
print(minskew)
#plt.plot(x1,x2)
#ax = plt.gca()
#ax.set_aspect(1.0)
#ax.set_xlabel('$x1(t)$')
#ax.set_ylabel('$x2(t)$')

#plt.plot(tau, S)
#ax = plt.gca()
#ax.set_xlabel('$\Delta t$')
#ax.set_ylabel('$S(\Delta t)$')

#plt.plot(x1,x2)
#ax = plt.gca()
#ax.set_aspect(1.0)
#ax.set_xlabel('$x1(t)$')
#ax.set_ylabel('$x2(t)$')

#plt.plot(tau, S)
#ax = plt.gca()
#ax.set_xlabel('$\Delta t$')
#ax.set_ylabel('$S(\Delta t)$')

plt.plot(x1[:-minskew], x2[minskew:])
ax = plt.gca()
ax.set_aspect(1.0)
ax.set_xlabel('$x1(t + \Delta t)$')
ax.set_ylabel('$x2(t)$')

show()

