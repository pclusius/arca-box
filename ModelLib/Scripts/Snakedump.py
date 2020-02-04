import numpy as np
import netCDF4
import sys; import os;

relpath = '../../INOUT/PC_0000/ODDFUNCTS/'
files = os.listdir(relpath)

ncfiles = []
for f in files:
    if f == 'General.nc':
        # print(f)
        ncfiles.append(f)

# use_first = input('Apply first file to all (y/n)?: ')
# if use_first == 'y': use_first = True
# else: use_first = False
# print(use_first)

first = netCDF4.Dataset(relpath+ncfiles[0], 'r')
vars = first.variables
chosenones = []
for v in vars:
    q = input('Include '+v+',with dimension: '+first[v].dimensions[0]+' (y/n): ')
    if q == 'y' or q == '' and q != 'n': chosenones.append(v)

if len(chosenones)>0:
    pass
else:
    xxx

for i,file in enumerate(ncfiles):
    nc = netCDF4.Dataset(relpath+file, 'r')
    matrix = np.zeros(  (first.dimensions['time'].size, len(chosenones)) )
    header=[]
    for j,variable in enumerate(chosenones):
        try:
            header.append(variable+'('+nc.variables[variable].units+')')
        except:
            header.append(variable+'([])')
        matrix[:,j] = nc.variables[variable][:]

    for z,u in enumerate(header):
        if z==0: headerstr = '%24s'%u
        else: headerstr = headerstr + '%30s'%u
    np.savetxt(relpath+file[:-3]+'_ncd.txt', matrix, delimiter='    ', fmt="%1.20e", header=headerstr)
