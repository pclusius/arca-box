import numpy as np
import matplotlib.pyplot as plt
import netCDF4
from matplotlib.widgets import RadioButtons, Button
# import sys; import os; sys.path.append(os.path.abspath("/home/pecl/p_modules"))
# import custom_lib as c

class compound:
    def __init__(self, var):
        self.name = var
    def pl(self):
        plot_var(self.name)
        return


def roundradios(radioax, bt, scale):
    rpos = radioax.get_position().get_points()
    fh = fig.get_figheight()
    fw = fig.get_figwidth()
    rscale = (rpos[:,1].ptp() / rpos[:,0].ptp()) * (fh / fw)
    for circ in bt.circles:
        circ.height /= rscale/scale
        circ.width *= scale

def uniq(list):
    chk = []
    for element in list:
        if element not in chk:
            chk.append(element)
    return chk
pth = '/media/pecl/EE7B-B587/PC180416/N1/'
pth = ''
pth = '../output/'
files = (['Test_01_general.nc'])

import sys
cmdline = False
if len(sys.argv) > 1:
    files = [str(sys.argv[1])]
    pth=''
    cmdline = True


N_f = len(files)
ncs = []
first_file = 0
fig, ax = plt.subplots(figsize=(11,7))

for file in files:
    ncs.append(netCDF4.Dataset(pth+file, 'r'))

hnames = [i for i in ncs[0].variables]
dims = [i for i in ncs[0].dimensions]

# remove string variables and constants
cache = []
for i,v in enumerate(hnames):
    try:
        int(ncs[0].variables[v][0])
        cache.append(v)
    except:
        pass
hnames = []
for i,v in enumerate(cache):
    if (len(np.shape(ncs[0].variables[v][:])) >1):
        hnames.append(v)
    # elif (len(uniq(ncs[0].variables[v][:]))>1):
    if not 'Shifter' in v and not 'Multipl' in v:
        hnames.append(v)

for var in hnames:
    strg = var+"=compound('%s')"%(var)
    exec(strg)

ax.set_title(file)

if len(hnames)<60:
    plt.subplots_adjust(left=0.35, bottom=0.1)
else:
    plt.subplots_adjust(left=0.1, bottom=0.1)

for jj in range(N_f):
    try:
        plt.plot(ncs[0].variables[ncs[0].variables[hnames[0]].dimensions[0]],
            ncs[jj].variables[hnames[0]],lw=3, label=hnames[0]+' in '+files[jj][:-3])
    except:
        plt.plot(ncs[jj].variables[hnames[0]],lw=3, label=hnames[0]+' in '+files[jj][:-3])

plt.legend(loc='upper right')

if len(hnames)<60:
    rax = plt.axes([0.36, 0.78, 0.1, 0.1], facecolor='w', alpha=0.1)
else:
    rax = plt.axes([0.1, 0.78, 0.1, 0.1], facecolor='w', alpha=0.1)

radio = RadioButtons(rax, ('linear','log'), active=0)

if len(hnames)<60:
    axrads = plt.axes([0.05, 0.1, 0.2, 0.8], facecolor='w')
    radplt = RadioButtons(axrads, (hnames), active=0)
    roundradios(axrads, radplt,1)

roundradios(rax, radio,1)

def scalefunc(scale):
    ax.set_yscale(scale)

def plot_var(plott):
    sc = ax.get_yscale()
    ax.clear()
    for jj in range(N_f):
        try:
            ax.plot( ncs[jj].variables[ncs[jj].variables[plott].dimensions[0]], ncs[jj].variables[plott],lw=3, label=plott+' in '+files[jj][:-3])
        except:
            ax.plot(ncs[jj].variables[plott],lw=3, label=plott+' in '+files[jj][:-3])
    ax.set_title(files)
    ax.set_yscale(sc)
    ax.grid(True)
    ax.legend(loc='upper right')

def update(val):
    fig.canvas.draw_idle()

radio.on_clicked(scalefunc)
radio.on_clicked(update)

if len(hnames)<60:
    radplt.on_clicked(plot_var)
    radplt.on_clicked(update)

def NAMES():
    for i in range(len(hnames)//6 + 1):
        print(hnames[i*5+0],'  ',hnames[i*5+1],'  ',hnames[i*5+2],'  ',hnames[i*5+3],'  ',hnames[i*5+4],'  ',hnames[i*5+5])
    return
print("use <variables>.pl() to plot any variable. E.g. O3.plot() for Ozone. NAMES() to show all variables")


ax.grid(True)
if (not cmdline):
    plt.ion()
plt.show()

if (cmdline):
    raise SystemExit
