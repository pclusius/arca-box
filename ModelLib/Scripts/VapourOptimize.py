import netCDF4
import numpy as np
import sys
from pathlib import Path


def whatever(aa, bb):
    return aa/bb

N_vapours=200

opt_names = open('../Vapour_names_opt.dat','w')
opt_props = open('../Vapour_properties_opt.dat','w')

kb  =  1.38064852e-23
if len(sys.argv) > 1:
    file = sys.argv[1]
else:
    file = '/home/pecl/05-APCAD/supermodel-phase-1/INOUT/HYDE/PC_2018-04-01/BASE/Chemistry.nc'

try:
    nc = netCDF4.Dataset(Path(file), 'r')
except:
    print( 'Could not open the file, is it accessible?' )

C_all = {}
for comp in nc.variables.keys():
    C_all[comp] = np.mean(nc.variables[comp][50:])

m,a,b = np.genfromtxt('/home/pecl/05-APCAD/ChemistryPackage/Vapour_properties_all.dat', unpack=True)
orgs = np.genfromtxt('/home/pecl/05-APCAD/ChemistryPackage/Vapour_names_all.dat', dtype=str)
conc_sat = 10**(a-b/298.)

HOAindex = list(orgs).index('HOA')
testvalues = np.zeros(len(orgs))

for i,comp in enumerate(orgs):
    testvalues[i] = whatever(C_all[comp], conc_sat[i])

order = np.argsort(-testvalues)


for i in range(N_vapours):
    print(orgs[order][i], testvalues[order][i])
    opt_names.write(orgs[order][i]+'\n')
    opt_props.write('%24.12f\t%24.12f\t%24.12f\n' %(m[order][i],a[order][i],b[order][i]))

opt_names.write('HOA\n')
opt_props.write('%24.12f\t%24.12f\t%24.12f\n' %(m[order][HOAindex],a[order][HOAindex],b[order][HOAindex]))

opt_names.close()
opt_props.close()




#
