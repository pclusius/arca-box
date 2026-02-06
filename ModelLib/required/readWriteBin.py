from numpy import fromfile, array, float64
from numpy.random import random
import sys,os
import numpy as np
import matplotlib.pyplot as plt
# netCDF4
# import matplotlib.pyplot as plt
#
# spectra = False
# writeFiles = False
# readNames = False
# pathToBinaryFile = ''
# selectTime = nynynynyny
# chemistryPath = ''
def read_bin(file='TIMESERIES.r16'):
    f = open(file, 'rb')
    conc = fromfile(f, dtype=float64)
    return conc

def write_bin(file='TIMESERIES.r16',c=None):
    f = open(file, 'wb')
    f.write(c)
    return

# conc = conc*random(conc.shape[0])
# write_bin(sys.argv[1],conc)
# N = int(sys.argv[2])
# conc = read_bin(sys.argv[1])
# M = conc.reshape((len(conc)//N,N))
component_index=0
# Nc=2209
# N=101
# dataroot='cells'
def Rt(Ci,i,endfile='c_ts.r16'):
    C = np.zeros(Ci.shape)
    N=C.shape[0]
    for iz in range(C.shape[1]):
        for iy in range(C.shape[2]):
            C[:,iz,iy] = read_bin(f'{dataroot}/cells/{iz+1:02d}/{iy+1:02d}/{endfile}').reshape((N,Nc))[:,i]
    return C

dataroot = os.getenv("locationPlumeData")
arca = os.getenv("arca")

if dataroot is None:
    print("Do 'source wall2dSettings.sh' before applying this script!")
    exit()
pblh = float(os.getenv("pblh"))
boxRuntime = float(os.getenv("boxRuntime"))
ny = int(os.getenv("nHorisontalColumns"))
nz = int(os.getenv("nVerticalLayers"))
dx = float(os.getenv("dx"))
plumeZ = int(os.getenv("zIndexPlume"))-1
plumeY = int(os.getenv("yIndexPlume"))-1

Nc=read_bin(dataroot+'/initial/background/CFINAL_0000.r16').shape[0]

names=np.genfromtxt(f'{arca}/src/chemistry/Hyde_VAR_MCM/SPC_NAMES.txt',dtype=str)
N=read_bin(dataroot+'/cells/01/01/c_ts.r16').shape[0]//Nc
BLH = np.array([float(open(f'{dataroot}/k-values/K-values_{i:04d}.txt', 'r').readlines()[2].strip()) for i in range(N)])
C=np.zeros((N,nz,ny))
if len(sys.argv)>1:
    try:
        component_index=int(sys.argv[1])
        C = Rt(C,component_index)
        if component_index+1>len(names):
            name='particles'
        else:
            name = names[component_index]
    except:
        component_index=np.array([int(i) for i in sys.argv[1].split(',')])
        for ci in component_index:
            C = C + Rt(C,ci)
        if component_index[0]+1>len(names):
            name='particles'
        else:
            name = ','.join([names[ci] for ci in component_index])

h = np.arange(1,nz+1)*dx
#h=np.append([0,5],h)
y=np.arange(ny)*dx
t = np.arange(N)*boxRuntime
ws=5
x = t*ws/1000

f,ax=plt.subplots(4,figsize=(14,12))
p1 = ax[0].pcolormesh(x,y,C[:,:,:].sum(1).T)#,norm='log')#,vmin=C[0,:,:].sum(0).min())
p2 = ax[1].pcolormesh(x,h,C[:,:,:].sum(2).T)#,norm='log')#,vmin=C[0,:,:].sum(1).min())
ax[0].set_title(f'{name}, total time: {t[-1]} s')
ax[0].set_ylabel('Width (m)')
ax[1].set_ylabel('Altitude (m)')
ax[1].plot(x,BLH, lw=0.5, ls='--',c='w')
# ax[0].set_xlabel(f'distance in {ws} m/s (km)')
# ax[1].set_xlabel(f'distance in {ws} m/s (km)')
plt.colorbar(p1)
plt.colorbar(p2)


p1 = ax[2].pcolormesh(x,y,C[:,plumeZ,:].T)
p2 = ax[3].pcolormesh(x,h,C[:,:,plumeY].T)
ax[2].set_ylabel('Width (m)')
ax[3].set_ylabel('Altitude (m)')
# ax[2].set_xlabel(f'distance in {ws} m/s (km)')
ax[3].set_xlabel(f'distance in {ws} m/s (km)')
ax[3].plot(x,BLH, lw=0.5, ls='--',c='w')
plt.colorbar(p1)
plt.colorbar(p2)
ax[3].set_ylim(h[0],h[-1])
ax[1].set_ylim(h[0],h[-1])

plt.figure()
plt.plot(x,C[:,plumeZ,plumeY])

f3 = plt.figure(figsize=(6,6))
os.makedirs(f'{dataroot}/wallAni',exist_ok=True)
os.chdir(dataroot)
os.system(f'rm wallAni/*.png')
for i in range(N):
    mx,my = np.unravel_index(np.argmax(C[i,:,:]),(nz,ny))
    plt.title(f'{name}, total time: {t[-1]} s')
    plt.pcolormesh(y,h,C[i,:,:])
    plt.gca().annotate(f'{C[i,:,:].max():12.3e}', xy=(y[my],h[mx]), xytext=(100,100),color='w',
            arrowprops=dict(facecolor='white', shrink=0.001))
    plt.savefig(f'wallAni/{i:04d}.png',dpi=150)
    plt.gca().clear()
plt.close(f3)
os.system("ffmpeg -y -framerate 5 -pattern_type glob -i 'wallAni/*.png'  -vf \"framerate=fps=30:interp_start=64:interp_end=164:scene=100\" wallani.mp4")


# plt.figure()
# plt.plot(C[:,:70,:].sum((1,2)) -C[0,:70,:].sum(),  label='bottom')
# plt.plot(C[:,70:,:].sum((1,2)) -C[0,70:,:].sum(),  label='top')
# plt.legend()


plt.figure()
# plt.ion()
[plt.plot(C[i,:,:].sum(1),h) for i in np.arange(0,N,10)]

plt.figure()
plt.plot(np.array([np.trapz(C[i,:,:].sum(1),h) for i in range(N)])/np.trapz(C[0,:,:].sum(1),h))

plt.show()







#
