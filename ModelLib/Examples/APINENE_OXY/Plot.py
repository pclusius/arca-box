import matplotlib.pyplot as plt
import numpy as np
import sys
import os

kb  =  1.38064852e-23 #[m2 kg /s/K]
Na  =  6.02214086e+23 # [1/mol]
Rg  =  8.31446
C0_ap = 38.3
t_obs, m_obs = np.genfromtxt('MeasuredMass.txt', unpack=True)
t_cal, m_cal = np.genfromtxt('calculatedMass.txt', unpack=True)
raw_x, raw_y = np.genfromtxt('plot.txt', unpack=True)

yyy =  raw_y[:len(t_obs)]
d_raw = raw_y[len(t_obs)]-raw_y[len(t_obs)+10]
d_y   = m_cal[0]-m_cal[10]
yyy = yyy*(d_y/d_raw)
yyy = yyy - (yyy[0] - m_cal[0])
fontsize=14
sfx = ''
from netCDF4 import Dataset as ncd

files = ['../../../INOUT/PAPER_0009/EXAMPLE_APINENE_OXY/Particles.nc',]

shift = 1
save_fig=True
plot_losses=True

legends=[
'Simulated total particle mass',
]
ls=[
'-',
'--',
]

f,ax = plt.subplots(figsize=(8,5))
if plot_losses: ax2=ax.twinx()
for i,run in enumerate(files):
    print(run)
    nc = ncd(run)
    time = nc.variables['TIME_IN_HRS'][:]+0.317*shift
    number = nc.variables['NUMBER_CONCENTRATION'][:,:]*1e6
    names=[]
    for n in nc.variables.keys():
        try: t = int(nc.variables[n].type)
        except: t=0
        if t>0: names.append(n)


    diam = nc.variables['DIAMETER'][0,:]
    mass = np.sum(np.sum(nc.variables['PARTICLE_COMPOSITION'][:,:,:], axis=(2))*number, axis=1)*1e9
    depmass = np.sum(nc.variables['DEPOSITED_PAR_COMP'][:,:],1)*1e9
    meanloss = ((nc.variables['PARTICLE_LOSS_RATE'][:,:])*3600)
    meanflux = np.sum(nc.variables['MASS_FLUX_ON_PAR'][:,:], 1)*1e9 # kg -> µg
    nc2=ncd(run.replace('Particles', 'General'))
    diam_a = np.argmin(abs(diam-50e-9))
    diam_b = np.argmin(abs(diam-500e-9))
    meanloss = np.mean((nc.variables['PARTICLE_LOSS_RATE'][:,diam_a:diam_b+1])*3600, 1)

    evaporation = np.where(meanflux<0, -1*meanflux, 0.0 )
    condensation = np.where(meanflux>0, meanflux, 0.0 )
    cumulative_vaporation = np.zeros(evaporation.shape)
    cumulative_condn = np.zeros(condensation.shape)
    cumulative_vaporation[1:] = np.cumsum(evaporation[1:]*np.diff(time)*3600)
    cumulative_condn[1:] = np.cumsum(condensation[1:]*np.diff(time)*3600)
    try:
        meanChemflux = np.sum(nc.variables['MASS_FLUX_ON_WALLS'][:,:], 1)*1e9 # kg -> µg
        Ch_Vol = nc.Chamber_volume
    except:
        meanChemflux = depmass * 0.0e0
        Ch_Vol = 1
    evapFromWalls = np.where(meanChemflux<0, -1*meanChemflux, 0.0 )
    condToWalls = np.where(meanChemflux>0, meanChemflux, 0.0 )
    cumulative_ChemVaporation = np.zeros(evapFromWalls.shape)
    cumulative_ChemCondn = np.zeros(condToWalls.shape)
    cumulative_ChemVaporation[1:] = np.cumsum(evapFromWalls[1:]*np.diff(time)*3600)
    cumulative_ChemCondn[1:] = np.cumsum(condToWalls[1:]*np.diff(time)*3600)

    # ax.plot(time, minmass+cumulative_condn,label='Condensed mass', linestyle=':',color='red', lw=2)

    # the simulated particle mass measurement
    ax.plot(time, mass,label=legends[i], linestyle=ls[i], lw=2)

    # the evaporation from suspended particles
    ax.fill_between(time, mass,mass+cumulative_vaporation,label='Mass evaporated from suspended particles', color='blue', alpha=0.5)

    # Particle wall losses
    ax.fill_between(time, mass, mass+depmass,label='Particle wall losses', color='grey',alpha=0.3)

    # Particle wall loss rate
    if plot_losses: ax2.plot(time, meanloss,label='Mean loss rate', linestyle='--',c='grey')

    # Particle wall loss rate
    ax.fill_between(time, mass+depmass, mass+depmass+cumulative_ChemCondn/Ch_Vol, \
                        linewidth=0, facecolor='none',edgecolor='forestgreen',hatch='\\\\\\',\
                        label='Chemical wall loss of vapors')

    APMass = C0_ap*1e-9*(1e5/(273.15+20)/Rg)*136.238e6 #µg/m³
    print('AMF from condensation/evaporation:', (cumulative_condn[-1]-cumulative_vaporation[-1])/APMass)
    print('AMF from condensation only:', (cumulative_condn[-1])/APMass)

zero_time_obs = t_obs#-cit
zero_time_cal = t_cal#-cit
ax.scatter(zero_time_obs[zero_time_obs>0], yyy[zero_time_obs>0],s=9,edgecolor='k',marker='s',facecolor='none', label='Measured')
ax.scatter(zero_time_cal[zero_time_cal>0], m_cal[zero_time_cal>0], s=9,edgecolor='k',marker='^',facecolor='none',label='Meas, par. losses removed')
ax.set_xlim(0,6)
ax.set_ylim(0,100)
ax.set_xlabel('Time [h]',fontsize=fontsize)
ax.set_ylabel('Mass [µg/m³]',fontsize=fontsize)
if plot_losses: ax2.set_ylabel('Particle loss rate [/h]',fontsize=fontsize)

handles1, labels1 = ax.get_legend_handles_labels()
if plot_losses: handles2, labels2 = ax2.get_legend_handles_labels()
if plot_losses:
    ax.legend(handles1+handles2, labels1+labels2)
else:
    ax.legend(handles1, labels1)
ax.tick_params(axis='both', which='major', labelsize=fontsize*0.75)
if plot_losses: ax2.tick_params(axis='both', which='major', labelsize=fontsize*0.75)

plt.tight_layout()
if save_fig: plt.savefig(f'PathakChamber_{shift:d}_{sfx:s}.png', dpi=300)
plt.show()
