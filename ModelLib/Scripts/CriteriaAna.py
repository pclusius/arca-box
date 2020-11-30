import netCDF4
import numpy as np
import matplotlib.pyplot as plt
import sys
from pathlib import Path
import scipy as sc
from scipy import stats

parfile = '/home/pecl/05-ARCA/supermodel-phase-1/INOUT/HYDE/PC_2018-04-04/OPTIFULL/Particles.nc'
original_list = '/home/pecl/05-ARCA/supermodel-phase-1/ModelLib/Vapour_names_massive.dat'

testCases = [
'/home/pecl/05-ARCA/supermodel-phase-1/ModelLib/Vapour_names_1860.dat',
'/home/pecl/05-ARCA/supermodel-phase-1/ModelLib/Vapour_names_C2.dat',
'/home/pecl/05-ARCA/supermodel-phase-1/ModelLib/Vapour_names_massive.dat',
'/home/pecl/05-ARCA/supermodel-phase-1/ModelLib/Vapour_names_MASSIVX_04.dat',
]
c = ['r.', 'b.', 'g.', 'y.']
labels = ['x=y=1','x=2,y=1','Aero-opt, other day', 'Aero-opt, same day']
nc = netCDF4.Dataset(parfile, 'r')
# org_compo = np.zeros((np.shape(nc.variables['PARTICLE_COMPOSITION'][:,:,:-2])[1:]))

n_org = np.shape(nc.variables['PARTICLE_COMPOSITION'][:,:,:])[2]
org_compo = np.sum(np.dstack([nc.variables['NUMBER_CONCENTRATION'][-1,:]]*n_org)[0,:,:-2]*nc.variables['PARTICLE_COMPOSITION'][-1,:,:-2], 0)

# for j in range(np.shape(nc.variables['PARTICLE_COMPOSITION'][:,:,:-2])[2]):
#     org_compo[:,j] = nc.variables['PARTICLE_COMPOSITION'][-1,:,j]*nc.variables['NUMBER_CONCENTRATION'][-1,:]
#
# org_compo = np.sum(org_compo, axis=0)
sorter = np.argsort(-org_compo)




vap_names = np.genfromtxt(original_list, dtype=str)


for i,file in enumerate(testCases):
    testvapours = np.genfromtxt(file, dtype=str)
    for n in range(1,800):
        if n==1:
            plt.plot(n,len(list(set(testvapours[:n]).intersection(vap_names[sorter][:n])))/n*100, c[i], label=labels[i])
        else:
            plt.plot(n,len(list(set(testvapours[:n]).intersection(vap_names[sorter][:n])))/n*100, c[i])
        # if n == 30:
        #     print('sorted by actual', vap_names[sorter][:n])
        #     print('criteria order', testvapours[:n])
plt.xlabel('Number of vapours considered')
plt.ylabel('% of selected vapours in actual contribution')
plt.title('Criteria: conc^x/Psat^y and Aerosol optimized')
plt.grid()
plt.legend()
plt.show()

    # compo = np.sum(np.sum(nc.variables['PARTICLE_COMPOSITION'][:,:,:-2], axis=2)*nc.variables['NUMBER_CONCENTRATION'][:,:], axis=1)



#
