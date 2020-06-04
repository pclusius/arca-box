import netCDF4
import numpy as np
import matplotlib.pyplot as plt
import sys
from pathlib import Path
import scipy as sc
from scipy import stats

kb  =  1.38064852e-23

# if len(sys.argv) > 1:
#     file = sys.argv[1]
# else:
#     file = '/home/pecl/05-APCAD/supermodel-phase-1/INOUT/HYDE/PC_2018-04-04/V100/Particle.nc'
# try:
#     nc = netCDF4.Dataset(Path(file), 'r')
# except:
#     print( 'Could not open the file, is it accessible?' )
files = [
'/home/pecl/05-APCAD/supermodel-phase-1/INOUT/HYDE/PC_2018-04-04/V50/Particle.nc',
'/home/pecl/05-APCAD/supermodel-phase-1/INOUT/HYDE/PC_2018-04-04/V100/Particle.nc',
'/home/pecl/05-APCAD/supermodel-phase-1/INOUT/HYDE/PC_2018-04-04/V150/Particle.nc',
'/home/pecl/05-APCAD/supermodel-phase-1/INOUT/HYDE/PC_2018-04-04/V200/Particle.nc',
]
for file in files:
    nc = netCDF4.Dataset(Path(file), 'r')
    masses = np.sum(np.sum(nc.variables['PARTICLE_COMPOSITION'], axis=1)[:,:-2], axis=1)
    time = np.linspace(0,20,masses.shape[0])
    plt.plot(time, masses, label=file[-16:-12])
plt.ylim(0,)
plt.legend()
plt.show()

#
