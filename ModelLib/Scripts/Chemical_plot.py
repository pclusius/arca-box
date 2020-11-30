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

showParPlot = False
showParPlot = True

showChemPlot = False
showChemPlot = True

selected_days = ([
'2016-06-12',
'2016-10-02',
'2017-05-28',
'2017-06-03',
'2017-07-15',
'2018-05-27',
'2019-07-31',
])

ammOulanka = ([
1.663575e+08,
2.079469e+09,
1.393244e+09,
1.871522e+08,
2.079469e+08,
9.357608e+08,
5.094698e+08,
])
ammHyy = ([
1.663575e+08,
2.079469e+09,
2.162647e+09,
8.317874e+08,
2.474568e+09,
4.366884e+09,
4.782778e+08,
])
# 2019-07-31 212
# 2018-05-27 147
# 2017-07-15 196
# 2017-06-03 154
# 2017-05-28 148
# 2016-10-02 276
# 2016-06-12 164

sels = [0,1,1,2,1,2,2]
sels = [0,0,0,0,0,0,0]
sels = np.flip(sels)

selected_days = np.flip(selected_days)

casefolders = ([
'J_FIN',
# 'J_CHASE',
# 'J_SECOND',
# 'J_THIRD',
 ])

gencomp = [
'RESOLVED_BASE'
'RESOLVED_J_FACTR'
]

compounds = [
'H2SO4',
'SO2',
'OH'
]

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

fulltimeSA = []
fulltimeOH = []
fulltimeSO2 = []
fulltimeNH3 = []
fulltime = []


fiz, ax7 = plt.subplots(2,7, figsize=(12,3.5), sharey=False, sharex=True)
fiz.subplots_adjust(hspace=0.1, wspace=0.0, bottom = 0.2)
fiz.add_subplot(111, frameon=False)
plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
plt.xlabel('time of day (hrs)', fontsize=13)
plt.ylabel('Concentration (1/cm³)', fontsize=13)


fioh, axoh = plt.subplots(1,7, figsize=(12,3), sharey=False, sharex=True)
fioh.subplots_adjust(hspace=0.1, wspace=0.0, bottom = 0.25)
fioh.add_subplot(111, frameon=False)
plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
plt.xlabel('time of day (hrs)', fontsize=13)
plt.ylabel('Concentration (1/cm³)', fontsize=13)


fiamm, axamm = plt.subplots(1,7, figsize=(12,3), sharey=False, sharex=True)
fiamm.subplots_adjust(hspace=0.1, wspace=0.0, bottom = 0.25)
fiamm.add_subplot(111, frameon=False)
plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
plt.xlabel('time of day (hrs)', fontsize=13)
plt.ylabel('Concentration (1/cm³)', fontsize=13)

fifac, axfac = plt.subplots(1,7, figsize=(12,3), sharey=False, sharex=True)
fifac.subplots_adjust(hspace=0.1, wspace=0.0, bottom = 0.25)
fifac.add_subplot(111, frameon=False)
plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
plt.xlabel('time of day (hrs)', fontsize=13)
plt.ylabel('Jdma/Jtot', fontsize=13)


for id, date in enumerate(selected_days):
    for ic, casefolder in enumerate(casefolders):
        if not sels[id] == ic:
            continue

        title = date
        polku = root+'/HYDE_'+date+'/'+casefolder
        JD = c.date_to_julian(date, show=1)

        """ from here on the MODEL output is processed"""
        try:
            ncp = netCDF4.Dataset(polku+'/Particles.nc', 'r')
            num_conc = ncp.variables['NUMBER_CONCENTRATION'][:,:]
            diam     = ncp.variables['DIAMETER'][0,:]
            time     = ncp.variables['time_in_hrs'][:]
            vol_mod_um  = diam**3 * np.pi / 6.0 * 1e18
        except:
            print('Reading Particles failed')
            continue

        """ from here on the MODEL CHEMISTRY is processed"""
        try:
            ncc = netCDF4.Dataset(polku+'/Chemistry.nc', 'r')
            chtime = ncc.variables['time_in_hrs'][:]
        except:
            print('Reading chemistry failed')
            continue

        """ from here on the FILE GENERAL.NC is processed"""
        try:
            ncg = netCDF4.Dataset(polku+'/General.nc', 'r')
            gentime = ncg.variables['time_in_hrs'][:]
            nucl = ncg.variables['NUC_RATE_IN'][:]
            resolvedbase = True
        except:
            print('Reading General output failed')
            resolvedbase = False
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


        """ Plot Chemistry related """
        if showChemPlot:
            chemConc = np.empty((len(chtime), len(compounds)))
            # fig, ax = plt.subplots(len(compounds)+resolvedbase,1, figsize=(8,1.5*(resolvedbase+len(compounds))), sharex=True)
            for iii,ccc in enumerate(compounds):
                chemConc[:,iii] = ncc.variables[ccc][:]
                # ax[iii].plot(chtime, chemConc[:,iii], label=ccc)
                # ax[iii].legend()
                # ax[iii].grid()
            if resolvedbase:
                ammDMA = ncg.variables['RESOLVED_BASE'][:]
                FacDMA = ncg.variables['RESOLVED_J_FACTR'][:]
                filter = np.logical_and(nucl>0.001, np.logical_and(ammDMA<1e17, ammDMA>1))
                # ax[iii+1].semilogy(gentime[filter], ammDMA[filter], label='NH3',c='tomato')
                # ax[iii+1].semilogy(gentime[filter], ammHyy[id]*np.ones(len(gentime))[filter], lw=2,linestyle=':',c='tomato',label='EMEP weakly mean')

                # ax[iii+1].legend()
                # ax[iii+1].grid()
            # plt.savefig('/home/pecl/01-TUTKIMUS/Gradu/Casestudy/plots/'+date+'_chem.png', dpi=230)
            # fulltimeSA = np.append(fulltimeSA, chemConc[:,0])
            # fulltimeSO2 = np.append(fulltimeSO2, chemConc[:,1])
            # fulltimeOH = np.append(fulltimeOH, chemConc[:,2])
            # fulltimeNH3 = np.append(fulltimeNH3, chemConc[:,3])
            # fulltime = np.append(gentime, chtime)

            N = 7
            smoothSA = chemConc[:,0].copy()
            smoothSA[N-1:] = np.convolve(chemConc[:,0], np.ones((N,))/N, mode='valid')
            smoothSO2 = chemConc[:,1].copy()
            smoothSO2[N-1:] = np.convolve(chemConc[:,1], np.ones((N,))/N, mode='valid')

            ax7[1,6-id].semilogy(chtime, smoothSA, lw=2, c='tomato', label='H2SO4')
            ax7[0,6-id].set_ylim(1e4,5e7)
            ax7[0,6-id].grid(axis='y')
            ax7[0,6-id].set_title(date)
            ax7[0,6-id].semilogy(chtime, smoothSO2, lw=2, c='orange', label='SO2')
            ax7[0,6-id].set_ylim(5e7,5e10)
            ax7[1,6-id].set_ylim(5e4,5e7)
            axoh[6-id].set_ylim(2e5,2e7)
            axamm[6-id].set_ylim(5e7,1e12)
            axfac[6-id].set_ylim(0,100)
            ax7[1,6-id].grid(axis='y')
            ax7[1,6-id].set_xticks(np.arange(0,24,5))
            axamm[6-id].grid(axis='y')
            axamm[6-id].set_xticks(np.arange(0,24,5))
            axfac[6-id].grid(axis='y')
            axfac[6-id].set_xticks(np.arange(0,24,5))

            # if id==3: ax7[1,6-id].set_xlabel('time of day (hrs)', fontsize=13)
            # if id==6: ax7[1,6-id].set_ylabel('Concentration (1/cm³)', fontsize=13)
            for jee in range(7):
                if jee<6: ax7[0,6-jee].set_yticklabels([])
                if jee<6: ax7[1,6-jee].set_yticklabels([])
                if jee<6: axoh[6-jee].set_yticklabels([])
                if jee<6: axamm[6-jee].set_yticklabels([])
                if jee<6: axfac[6-jee].set_yticklabels([])

            axoh[6-id].set_title(date)
            axoh[6-id].semilogy(chtime, chemConc[:,2], lw=3, label='OH')
            axoh[6-id].grid(axis='y')
            axoh[6-id].set_xticks(np.arange(0,24,5))


            axamm[6-id].set_title(date)
            axamm[6-id].semilogy(gentime[filter], ammDMA[filter], label='Base',c='deepskyblue', lw=3)
            axamm[6-id].semilogy(gentime[filter], ammHyy[id]*np.ones(len(gentime))[filter], lw=2,linestyle=':',c='k',label='EMEP')

            axfac[6-id].set_title(date)
            axfac[6-id].plot(gentime[filter], FacDMA[filter]*100, label='dma factor',c='red', lw=3)








        """ Plot particle volume, surface plot, number conc """
        if showParPlot:
            fig, axes = plt.subplots(4,1, figsize=(3,8), sharex=True)

            """volume plot"""
            vol = np.sum(dm_linconc[:,:imames+1]*vol_mes_um[:imames+1],1)
            VOL = np.sum(num_conc[:,imimod:imamod+1]*vol_mod_um[imimod:imamod+1],1)

            N = 9
            vol[N-1:] = np.convolve(vol, np.ones((N,))/N, mode='valid')

            axes[0].plot(dmtime, vol, label='measured', markersize=3, marker='o', lw=0,markerfacecolor='None')
            axes[0].plot(time, VOL, label='modelled', lw=2, c='tomato')

            axes[0].set_title(title)
            axes[0].set_ylabel('Volume (µm³)')
            axes[0].set_ylim(0,0.1)

            """number conc plot"""
            NUM = np.sum(dm_conc[:,:imames+1],1)
            NUM[N-1:] = np.convolve(NUM, np.ones((N,))/N, mode='valid')

            axes[1].plot(dmtime, NUM, label='measured', marker='x', markersize=3, lw=0 )
            axes[1].plot(time, np.sum(mo_logconc[:,imimod:imamod+1],1), label='modelled',c='orange', lw=2)
            axes[1].set_ylabel('Total concentration (µm³)')
            axes[1].set_ylim(0,200000)


            X,Y,Z = c.prep_surf(dmdiam, dmtime, dm_conc)
            plt1 = axes[2].pcolormesh(X*24, Y, Z, vmin=0, vmax=4)
            axes[2].set_yscale('log')
            axes[2].set_ylim(3,1000)

            X,Y,Z = c.prep_surf(diam, time, mo_logconc)
            plt2 = axes[3].pcolormesh((23*X), Y, Z, vmin=0, vmax=4)
            secondax = axes[3].twinx()
            secondax.plot(gentime, nucl, lw=2, c='w')
            secondax.set_ylim(0,1.6)
            axes[3].set_yscale('log')
            axes[3].set_ylim(3,1000)
            axes[3].set_xlim(0,24)


            """Decorate plots"""
            # axes[0].set_xlim(0,24)
            # axes[1].set_xlim(0,24)
            # axes[0].grid()
            # axes[0].legend()
            # axes[1].grid()
            # axes[1].legend()
            # axes[3].set_xlabel('Time (hrs)')
            plt.tight_layout()
            # plt.savefig('/home/pecl/01-TUTKIMUS/Gradu/Casestudy/plots/'+date+casefolders[ic]+'_parts.png', dpi=230)
# fiz.tight_layout()
ax7[0,0].legend()
ax7[1,0].legend()
axoh[0].legend()
axamm[6].legend()
axfac[6].legend()
# fiz.savefig('/home/pecl/01-TUTKIMUS/Gradu/Casestudy/plots/SA_full.png', dpi=230)
# fioh.savefig('/home/pecl/01-TUTKIMUS/Gradu/Casestudy/plots/OH_full.png', dpi=230)
# fiamm.savefig('/home/pecl/01-TUTKIMUS/Gradu/Casestudy/plots/NH3_full_new.png', dpi=230)
# fifac.savefig('/home/pecl/01-TUTKIMUS/Gradu/Casestudy/plots/DMAFAC.png', dpi=230)
plt.show()



















#
