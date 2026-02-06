from numpy import fromfile, array, float64
from numpy.random import random
import sys,os
import numpy as np
import matplotlib.pyplot as plt

def read_bin(file='TIMESERIES.r16'):
    f = open(file, 'rb')
    conc = fromfile(f, dtype=float64)
    return conc

def write_bin(file='TIMESERIES.r16',c=None):
    f = open(file, 'wb')
    f.write(c)
    return

dia = np.logspace(np.log10(3e-9),np.log10(1e-5),100)

Nc=2209
component_index=np.arange(100)+Nc
Nc=2409
N=13
dataroot='cells'

def Rt(Ci,i,endfile='c_ts.r16'):
    C = np.zeros(Ci)
    N=C.shape[0]
    for iz in range(C.shape[1]):
        for iy in range(C.shape[2]):
            C[:,iz,iy] = read_bin(f'{dataroot}/{iz+1:02d}/{iy+1:02d}/{endfile}').reshape((N,Nc))[:,i]
    return C

if len(sys.argv)>1:
    Nc=int(sys.argv[2])
    N=int(sys.argv[3])
    # C=np.zeros((N,80,99))
    # component_index=np.array([int(i) for i in sys.argv[1].split(',')])
    # for ci in component_index:
    #     C = C + Rt(C,ci)
    # name='particles'


c0 = np.array([Rt((1,80,99),ci,endfile='CFINAL_0000.r16') for ci in component_index])
cf = np.array([Rt((1,80,99),ci,endfile=f'CFINAL_{N-1:04d}.r16') for ci in component_index])

plt.ion()


h = np.arange(1,79)*10
h = np.append([0,5],h)
y = np.arange(99)*10
t = np.arange(N)*30
ws = 5
x = t*ws/1000
