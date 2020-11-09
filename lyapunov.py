import numpy as np
from numpy import genfromtxt
from matplotlib.pylab import *
import matplotlib.pyplot as plt
from fancyplot import *

def lyapunov(_t,_R1, _z1, _phi1, _vR1, _vz1, _vphi1, _R2, _z2, _phi2, _vR2, _vz2, _vphi2):
    _x1=_R1*np.cos(_phi1)
    _y1=_R1*np.sin(_phi1)
    _vx1=_vR1*np.cos(_phi1)-_R1*np.sin(_phi1)*_vphi1
    _vy1=_R1*np.cos(_phi1)*_vphi1+_vR1*np.sin(_phi1)

    _x2=_R2*np.cos(_phi2)
    _y2=_R2*np.sin(_phi2)
    _vx2=_vR2*np.cos(_phi2)-_R2*np.sin(_phi2)*_vphi2
    _vy2=_R2*np.cos(_phi2)*_vphi2+_vR2*np.sin(_phi2)

    dim=np.amin(np.array([_R1.size,_R2.size]))
    _t=_t[1:dim+1]

    _x1=_x1[:dim]
    _y1=_y1[:dim]
    _z1=_z1[:dim]
    _vx1=_vx1[:dim]
    _vy1=_vy1[:dim]
    _vz1=_vz1[:dim]

    _x2=_x2[:dim]
    _y2=_y2[:dim]
    _z2=_z2[:dim]
    _vx2=_vx2[:dim]
    _vy2=_vy2[:dim]
    _vz2=_vz2[:dim]

    _d=np.zeros(dim)

    _d=np.sqrt((_x1-_x2)**2+(_y1-_y2)**2+(_z1-_z2)**2+(_vx1-_vx2)**2+(_vy1-_vy2)**2+(_vz1-_vz2)**2)
    print("fist d ", _d[0], "last d ", _d[-1], "last t", _t[-1])

    _lambda= (np.log(_d/_d[0]))/(_t)
    return _lambda

    
#presa dati

data1 = genfromtxt('ComplexPlummer_reg_0.txt', delimiter='	', skip_header=1, usecols=(0,1,2,3,4,5,6))
t1, R1, z1, phi1, vR1, vz1, vphi1 = data1[:,0], data1[:,1], data1[:,2], data1[:,3], data1[:,4], data1[:,5], data1[:,6]

data2 = genfromtxt('ComplexPlummer_reg_1.txt', delimiter='	', skip_header=1, usecols=(0,1,2,3,4,5,6))
t2, R2, z2, phi2, vR2, vz2, vphi2 = data2[:,0], data2[:,1], data2[:,2], data2[:,3], data2[:,4], data2[:,5], data2[:,6]

lyap_exp1=lyapunov(t1, R1, z1, phi1, vR1, vz1, vphi1, R2, z2, phi2, vR2, vz2, vphi2)



data1 = genfromtxt('ComplexPlummer_ch_0.txt', delimiter='	', skip_header=1, usecols=(0,1,2,3,4,5,6))
t1, R1, z1, phi1, vR1, vz1, vphi1 = data1[:,0], data1[:,1], data1[:,2], data1[:,3], data1[:,4], data1[:,5], data1[:,6]

data2 = genfromtxt('ComplexPlummer_ch_1.txt', delimiter='	', skip_header=1, usecols=(0,1,2,3,4,5,6))
t2, R2, z2, phi2, vR2, vz2, vphi2 = data2[:,0], data2[:,1], data2[:,2], data2[:,3], data2[:,4], data2[:,5], data2[:,6]

lyap_exp2=lyapunov(t1, R1, z1, phi1, vR1, vz1, vphi1, R2, z2, phi2, vR2, vz2, vphi2)


c=np.amin(np.array([lyap_exp1.size,lyap_exp2.size]))

fancylayout()

plt.xscale('log')
plt.xlabel('log(t)')
plt.yscale('log')
plt.ylabel('log($\lambda$)')
plt.plot(t1[1:c+1], lyap_exp1[:c])
plt.plot(t1[1:c+1], lyap_exp2[:c])
plt.show()
