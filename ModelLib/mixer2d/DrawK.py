import matplotlib.pyplot as plt
import numpy as np
import sys,os

import time

fig = plt.figure()
ax = fig.add_subplot(111)

c = 'r'
size = 10
fileout = 'settings/Kdraw.txt'
locked = False
Kdraw = []
Zdraw = []
points = [ax.scatter([],[],marker='.')]
dataK = []
current = None

if os.getenv("nVerticalLayers") is not None:
    nz = int(os.getenv("nVerticalLayers"))
    dx = float(os.getenv("dx"))
    pblh = float(os.getenv("pblh"))
    root = os.getenv("locationPlumeData")
    fileout = root+'/'+fileout
else:
    print("Do 'source wall2dSettings.sh' before applying this script!")
    exit()

print(f'Modifying K-values in {root}, nz: {nz}, dz: {dx} m')
Zmodel = np.linspace(0,nz*dx,nz)
ax.set_ylim(0,nz*dx)
ax.set_xlim(0,50)
plt.axhline(pblh, ls=':')

try:
    Kprev = np.genfromtxt(fileout)
    plt.plot(Kprev,Zmodel)
except:
    pass

def save(filename=fileout, index=-1):
    try:
        with open (fileout, 'w+') as file: file.write(' '.join([f'{i:8.3f}' for i in dataK[index]])+'\n\n')
    except:
        print('*failed to write file')

def lock_on(event):
    global locked , current, points#,dx,dy,dx0,dy0
    locked=True
    current = ax.scatter([],[],marker='.', label=f'{len(points)-1:01d}')
    # points.append(ax.scatter([],[],marker='.', label=f'{len(points)-1:01d}'))
    print('Recording ...')

def lock_off(event):
    global locked, Kdraw, Zdraw , points, current #,dx,dy,dx0,dy0,lengths,points
    locked=False
    if len(Zdraw)<3:
        current.remove()
        current = None
        print("Draw slower, more points than 2 needed...")
        return
    else:
        points.append(current)
    sorter = np.argsort(Zdraw)
    kint = np.interp(Zmodel, np.array(Zdraw)[sorter],np.array(Kdraw)[sorter])
    dataK.append(kint)
    Kdraw, Zdraw = [],[]
    points[-1].set_offsets(np.array([kint,Zmodel]).T)
    ax.legend()


def record(event):
    global locked, Kdraw, Zdraw  , points#,dx,dy,dx0,dy0,lengths,points
    if event.inaxes != ax: return
    if locked:
        Kdraw.append(event.xdata)
        Zdraw.append(event.ydata)
        # print(event.xdata,event.ydata)
        current.set_offsets(np.array([Kdraw,Zdraw]).T)
        time.sleep(50/1000.0)

cid3 = fig.canvas.mpl_connect('button_press_event', lock_on)
cid4 = fig.canvas.mpl_connect('button_release_event', lock_off)
cid4 = fig.canvas.mpl_connect('motion_notify_event', record)

plt.ion()
plt.tight_layout()
plt.show()
