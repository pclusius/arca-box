import matplotlib
#matplotlib.use('TKAgg')
import pdb; bp=pdb.set_trace

import matplotlib.pyplot as plt
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import zipfile as zip
import gzip
import sys; import os; sys.path.append(os.path.abspath("/home/pecl/p_modules"))
import custom_lib as c
# from netCDF4 import Dataset
import netCDF4

import scipy as sc
import scipy.signal
import scipy.interpolate
from scipy import optimize
from scipy.integrate import quad, simps, quadrature
from scipy.misc import derivative
# from matplotlib import style
# style.use('seaborn-dark')
# plt.rcParams['image.cmap'] = 'jet'

# import scipy as sc
from scipy import optimize
# from scipy.integrate import quad, simps, quadrature
# from scipy.misc import derivative

import matplotlib.gridspec as gridspec
from matplotlib.colorbar import Colorbar

def filter_nan(mod_data, obs_data):
    """
    this functions removed the data  from modeled and observed data
    whereever the observed data contains nan

    this is used by all other functions, otherwise they will produce nan as
    output
    """
    data = np.array([mod_data.flatten(), obs_data.flatten()])
    data = np.transpose(data)
    data = data[~np.isnan(data).any(1)]
    return (data[:, 0], data[:, 1])

def rmse(mod_data, obs_data):
    """Root Mean Squared Error"""
    mod_data, obs_data = filter_nan(mod_data, obs_data)
    return np.sqrt(np.mean((mod_data - obs_data) ** 2))


def mae(mod_data, obs_data):
    """Mean Absolute Error"""
    mod_data, obs_data = filter_nan(mod_data, obs_data)
    return np.mean(abs(mod_data - obs_data))


def bias(mod_data, obs_data):
    """Bias"""
    mod_data, obs_data = filter_nan(mod_data, obs_data)
    return np.mean(mod_data - obs_data)


def KK_eq(J,CS, GR):
    gamma = 0.23
    CSp = CS / 4 / np.pi / 1e-5
    J3 = J*np.exp( gamma * (1/3 - 1/1.17) * CSp / GR)
    return J3

def KK_J1(J,CS, GR):
    gamma = 0.23
    CSp = CS / 4 / np.pi / 1e-5
    J1 = J*np.exp( - gamma * (1/3 - 1/1.17) * CSp / GR)
    return J1

col = ("#E69F00", "#3771C8", "#009E73", "#§4247ef", "#0072B2", "#D55E00", "#CC79A7")

root ='/home/pecl/05-ARCA/supermodel-phase-1/INOUT/GRADU'

distro = 0

show = False
show = True

selected_days = ([
# '2016-06-12',
'2016-10-02',
# '2017-05-28',
# '2017-06-03',
# '2017-07-15',
# '2018-05-27',
# '2019-07-31',
])

selected_days = np.flip(selected_days)

casefolders = ([
'J_SECOND',
'J_THIRD',
 ])

dates = selected_days.copy()

if len(sys.argv) > 1:
    casefolders = [str(sys.argv[2])]
    selected_days = [str(sys.argv[1])]
    root ='/home/pecl/05-ARCA/supermodel-phase-1/INOUT/GRADU'

Gases = ([
# '/gas_SO2.dat',
])

col_meas = 'tomato'
col_mod  = 'deepskyblue'

for i,d in enumerate(selected_days):
    dates[i] = 'dm'+d[2:4]+d[5:7]+d[8:]

dmpath = '/home/pecl/01-TUTKIMUS/Gradu/Casestudy/filleddmps'

def idiam(diam, target):
    return np.argmin(abs(diam-target))

def log_to_lin(diam, sum):
    log = np.log10(diam)
    FF = (log[2:]-log[1:-1])/2 + (log[1:-1]-log[0:-2])/2
    dm_nconc_f = sum.copy()
    dm_nconc_f[:,1:-1] = dm_nconc_f[:,1:-1]*FF
    dm_nconc_f[:,0] = dm_nconc_f[:,0]*(log[1]-log[0])
    dm_nconc_f[:,-1] = dm_nconc_f[:,-1]*(log[-1]-log[-2])
    return dm_nconc_f

def lin_to_log(diam, sum):
    log = np.log10(diam)
    FF = (log[2:]-log[1:-1])/2 + (log[1:-1]-log[0:-2])/2
    dm_nconc_f = sum.copy()
    dm_nconc_f[:,1:-1] = dm_nconc_f[:,1:-1]/FF
    dm_nconc_f[:,0] = dm_nconc_f[:,0]/(log[1]-log[0])
    dm_nconc_f[:,-1] = dm_nconc_f[:,-1]/(log[-1]-log[-2])
    return dm_nconc_f

for id, date in enumerate(selected_days):
    for ic, casefolder in enumerate(casefolders):

        title = date
        polku = root+'/HYDE_'+date+'/'+casefolder
        JD = c.date_to_julian(date, show=1)


        """ from here on the MODEL output is processed"""
        try:
            nc = netCDF4.Dataset(polku+'/Particles.nc', 'r')
            num_conc = nc.variables['NUMBER_CONCENTRATION'][:,:]
            diam     = nc.variables['DIAMETER'][0,:]
            time     = nc.variables['time_in_hrs'][:]
            vol_mod_um  = diam**3 * np.pi / 6.0 * 1e18
        except:
            continue

        """ from here on the dmps measurements are processed"""

        file = dmpath +'/'+ dates[id]+'.sum'
        dmdata = np.genfromtxt(file)
        dmdiam = dmdata[0,2:]
        dmtime = dmdata[1:,0]
        dmtime = (dmtime-dmtime[0])*24
        dm_conc = dmdata[1:,2:]
        dm_sum =  dmdata[1:,1]

        vol_mes_um  = dmdiam**3 * np.pi / 6.0 * 1e18

        dm_linconc = log_to_lin(dmdiam, dm_conc)
        mo_logconc = lin_to_log(diam, num_conc)


        """define the range of particle diamters considewreds in volume plot"""
        imimod = idiam(dmdiam,3e-9)
        imames = idiam(dmdiam,35e-9)
        imamod = idiam(diam,35e-9)

        """ Figure """
        fig, axes = plt.subplots(4,1, figsize=(8,8), sharex=False)

        """volume plot"""
        vol = np.sum(dm_linconc[:,:imames+1]*vol_mes_um[:imames+1],1)
        VOL = np.sum(num_conc[:,imimod:imamod+1]*vol_mod_um[imimod:imamod+1],1)

        N = 9
        vol[N-1:] = np.convolve(vol, np.ones((N,))/N, mode='valid')

        axes[0].plot(dmtime, vol, label='measured')
        axes[0].plot(time, VOL, label='modelled')

        axes[0].set_title(title)
        axes[0].set_xlabel('Time (hrs)')
        axes[0].set_ylabel('Volume (µm³)')

        """number conc plot"""
        NUM = np.sum(dm_conc[:,:imames+1],1)
        NUM[N-1:] = np.convolve(NUM, np.ones((N,))/N, mode='valid')

        axes[1].plot(dmtime, NUM, label='measured')
        axes[1].plot(time, np.sum(mo_logconc[:,imimod:imamod+1],1), label='modelled')

        axes[1].set_xlabel('Time (hrs)')
        axes[1].set_ylabel('Total concentration (µm³)')

        X,Y,Z = c.prep_surf(dmdiam, dmtime, dm_conc)
        plt1 = axes[2].pcolormesh(X+JD, Y, Z, vmin=0, vmax=4)
        axes[2].set_yscale('log')
        axes[2].set_ylim(3,1000)

        X,Y,Z = c.prep_surf(diam, time, mo_logconc)
        plt2 = axes[3].pcolormesh(X+JD, Y, Z, vmin=0, vmax=4)
        axes[3].set_yscale('log')
        axes[3].set_ylim(3,1000)
        axes[3].set_xlim(time[0]+JD,time[0]+JD+1)

        # cb = Colorbar(ax = axes[4], mappable = plt1, orientation = 'vertical', ticklocation = 'left')
        # plt.tight_layout()
        # plt.savefig('/home/pecl/01-TUTKIMUS/Gradu/Casestudy/plots/bar.png', dpi=230)
        """Decorate plots"""
        # for a in axes:
        #     a.legend()
        #     a.grid()

plt.show()



















#
