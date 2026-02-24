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

def Rt(Ci,i,endfile='c_ts.r16'):
    C = np.zeros(Ci.shape)
    N=C.shape[0]
    for iz in range(C.shape[1]):
        for iy in range(C.shape[2]):
            C[:,iz,iy] = read_bin(f'{dataroot}/{case}/{iz+1:03d}/{iy+1:03d}/{endfile}').reshape((N,Nc))[:,i].sum(1)
    return C

dataroot = os.getenv("locationPlumeData")
case = os.getenv("CaseCells")
arca = os.getenv("arca")

if dataroot is None:
    print("Do 'source wall2dSettings.sh' before applying this script!")
    exit()
pblh = float(os.getenv("pblh"))
boxRuntime = float(os.getenv("boxRuntime"))
ny = int(os.getenv("nHorisontalColumns"))
nz = int(os.getenv("nVerticalLayers"))
dx = float(os.getenv("dx"))
ddz = float(os.getenv("ddz"))
if os.getenv("zIndexPlumeStr") != '':
    plumeZ = np.array([int(i)-1 for i in os.getenv("zIndexPlumeStr").split(',')]) # int(os.getenv("zIndexPlume"))-1
else:
    plumeZ = range(nz)
if os.getenv("yIndexPlumeStr") !='':
    plumeY = np.array([int(i)-1 for i in os.getenv("yIndexPlumeStr").split(',')]) # int(os.getenv("zIndexPlume"))-1
else:
    plumeY = range(ny)

Nc=read_bin(dataroot+'/initial/background/CFINAL_0000.r16').shape[0]

names=np.genfromtxt(f'{arca}/src/chemistry/Hyde_VAR_MCM/SPC_NAMES.txt',dtype=str)
nd = {n:i for i,n in enumerate(names)}

N=read_bin(dataroot+f'/{case}/001/001/c_ts.r16').shape[0]//Nc
Kz = np.zeros((N,nz))
if ny>2:
    Ky = np.zeros((N,ny))
for i_t in range(N):
    Kz[i_t,:] = [float(ik) for ik in  open(f'{dataroot}/{case}/k-values/K-values_{i_t:04d}.txt', 'r').readlines()[0].strip().split()]
    if ny>2:
        Ky[i_t,:] = [float(ik) for ik in  open(f'{dataroot}/{case}/k-values/K-values_{i_t:04d}.txt', 'r').readlines()[1].strip().split()]


BLH = np.array([float(open(f'{dataroot}/{case}/k-values/K-values_{i:04d}.txt', 'r').readlines()[2].strip()) for i in range(N)])/1000
Ustar = np.array([float(open(f'{dataroot}/{case}/k-values/K-values_{i:04d}.txt', 'r').readlines()[3].strip()) for i in range(N)])
KyKz = np.array([float(open(f'{dataroot}/{case}/k-values/K-values_{i:04d}.txt', 'r').readlines()[4].strip()) for i in range(N)])
C=np.zeros((N,nz,ny))

if len(sys.argv)>1:
    tt = sys.argv[1].split(',')
    try:
        component_index = [int(t) for t in tt]
    except:
        component_index = [int(nd[t]) for t in tt]
    if len(component_index)==1:
        component_index = component_index[0]
    else:
        component_index = ','.join([f'{i}' for i in component_index])
    try:
        component_index=int(component_index)
        C = Rt(C,np.array([component_index]))
        if component_index+1>len(names):
            name='particles'
        else:
            name = names[component_index]
    except:
        component_index=np.array([int(i) for i in component_index.split(',')])
        # for ci in component_index:
        C = Rt(C,component_index)
        if component_index[0]+1>len(names):
            name='particles'
        else:
            name = ','.join([names[ci] for ci in component_index])

# h = np.arange(1,nz+1)*dx/1000
# h =

dz=np.linspace(0,(nz-2)*ddz,nz-1)
dz0 = dx
h = np.append(0,np.cumsum(dz+dz0))/1000

y = np.arange(ny)*dx
y = (y - y[plumeY].mean())/1000
t = np.arange(N)*boxRuntime
ws = 5
x = t*ws/1000
#fig 1
f,(ax)=plt.subplots(2,2,figsize=(14,12))
ax=ax.flatten()
if ny>2:
    p1 = ax[0].pcolormesh(x,y,C[:,:,:].sum(1).T,cmap='nipy_spectral',norm='log')#,vmin=C[0,:,:].sum(0).min())
    plt.colorbar(p1, orientation='horizontal')
p2 = ax[1].pcolormesh(x,h,C[:,:,:].sum(2).T,cmap='nipy_spectral',norm='log')#,vmin=C[0,:,:].sum(1).min())
ax[0].set_title(f'{case}: {name}, total time: {t[-1]} s')
ax[0].set_ylabel('Width (km)')
ax[1].set_ylabel('Altitude (km)')
ax[1].plot(x,BLH, lw=0.5, ls='--',c='w')
# ax[0].set_xlabel(f'distance in {ws} m/s (km)')
# ax[1].set_xlabel(f'distance in {ws} m/s (km)')
plt.colorbar(p2, orientation='horizontal')


if ny>2:
    p1 = ax[2].pcolormesh(x,y,C[:,plumeZ,:].mean(1).squeeze().T,cmap='nipy_spectral')
    plt.colorbar(p1, orientation='horizontal')
p2 = ax[3].pcolormesh(x,h,C[:,:,plumeY].mean(2).squeeze().T,cmap='nipy_spectral')
ax[2].set_ylabel('Width (km)')
ax[3].set_ylabel('Altitude (km)')
# ax[2].set_xlabel(f'distance in {ws} m/s (km)')
ax[3].set_xlabel(f'distance in {ws} m/s (km)')
ax[3].plot(x,BLH, lw=0.5, ls='--',c='w')
plt.colorbar(p2, orientation='horizontal')
ax[3].set_ylim(h[0],h[-1])
ax[1].set_ylim(h[0],h[-1])

# [x.set_aspect('equal',anchor='SW',share=True) for x in ax]

#fig 2
plt.figure()
plt.title(f'[c] on the "plume axis" {name}')
plt.xlabel('km')
plt.ylabel('cm^-3')
plt.plot(x,C[:,plumeZ,:].mean(1)[:,plumeY].mean(1))

f3 = plt.figure(figsize=(6,4))
os.makedirs(f'{dataroot}/wallAni',exist_ok=True)
os.chdir(dataroot)
os.system(f'rm -f wallAni/*.png')

for i in range(N):
    mx,my = np.unravel_index(np.argmax(C[i,:,:]),(nz,ny))
    plt.title(f'{os.path.split(dataroot)[1]}, {case}: {name}, time: {t[i]}s, kY/kZ:{KyKz[i]:0.2f}, u*:{Ustar[i]:0.2f}')
    mmx = C[i,:,:].max()
    plt.pcolormesh(y,h,C[i,:,:], cmap='nipy_spectral')
    plt.gca().annotate(f'{C[i,:,:].max():12.3e}', xy=(y[my],h[mx]), xytext=(y[int(ny/10)],h[int(nz/10)]),color='w',
            arrowprops=dict(facecolor='white', shrink=0.001))
    [plt.scatter(0,h[pZ], marker='+',s=20) for pZ in plumeZ]
    plt.gca().set_aspect('equal')
    plt.savefig(f'wallAni/{i:04d}.png',dpi=150)
    plt.gca().clear()
plt.close(f3)
os.system("ffmpeg -y -framerate 5 -pattern_type glob -i 'wallAni/*.png'  -vf \"framerate=fps=30:interp_start=64:interp_end=164:scene=100\" wallani.mp4")

# fig 3
plt.figure()
plt.title(f'Vertical [c] profiles for time instances {name}')
plt.xlabel('cm^-3')
plt.ylabel('m')
[plt.plot(C[i,:,:].sum(1),h,label=t[i]) for i in np.arange(0,N,10)]
plt.legend()
# fig 4
plt.figure()
plt.title(f'C(t)/C(t=0) {name}')
plt.xlabel('km')
plt.ylabel('[-]')
plt.plot(x,np.array([np.trapz(C[i,:,:].sum(1),h) for i in range(N)])/np.trapz(C[0,:,:].sum(1),h))

if ny>2:
    f,ax=plt.subplots(2,figsize=(14,12))
    p1 = ax[0].pcolormesh(x,y,Ky[:,:].T,cmap='nipy_spectral',vmin=0)
    p2 = ax[1].pcolormesh(x,h,Kz[:,:].T,cmap='nipy_spectral',vmin=0)
    plt.colorbar(p1)
else:
    f,ax=plt.subplots(1,figsize=(14,12))
    p2 = ax.pcolormesh(x,h,Kz[:,:].T,cmap='nipy_spectral',vmin=0)
plt.colorbar(p2)

plt.show()


if ny>2:
    import evtk as ev
    nx = N
    parts = np.ones((nz,nx,ny))
    # parts = np.zeros((nz,nlat,nlon))

    xo = np.ones(parts.shape)
    yo = np.ones(parts.shape)
    zo = np.ones(parts.shape)

    for k in range(nz):
        xo[k,:,:] = np.meshgrid(y,x)[0]
        yo[k,:,:] = np.meshgrid(y,x)[1]
        zo[k,:,:] = h[k]
        parts[k,:,:] = C[:,k,:]*1e9/2.5e19

    ev.hl.gridToVTK('/home/pecl/3dsolid',xo,yo,zo,pointData = {name : parts})

    #
