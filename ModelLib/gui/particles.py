from numpy import log10, genfromtxt
from pathlib import Path

def log_to_lin(diam, sum):
    log = log10(diam)
    FF = (log[2:]-log[1:-1])/2 + (log[1:-1]-log[0:-2])/2
    dm_nconc_f = sum.copy()
    dm_nconc_f[:,1:-1] = dm_nconc_f[:,1:-1]*FF
    dm_nconc_f[:,0] = dm_nconc_f[:,0]*(log[1]-log[0])
    dm_nconc_f[:,-1] = dm_nconc_f[:,-1]*(log[-1]-log[-2])
    return dm_nconc_f

def lin_to_log(diam, sum):
    log = log10(diam)
    FF = (log[2:]-log[1:-1])/2 + (log[1:-1]-log[0:-2])/2
    dm_nconc_f = sum.copy()
    dm_nconc_f[:,1:-1] = dm_nconc_f[:,1:-1]/FF
    dm_nconc_f[:,0] = dm_nconc_f[:,0]/(log[1]-log[0])
    dm_nconc_f[:,-1] = dm_nconc_f[:,-1]/(log[-1]-log[-2])
    return dm_nconc_f


def parseSum(file):
    try:
        data = genfromtxt(Path(file))
    except:
        return 'Could not open the file: '+data
    if data[0,1] == 0:
        n_conc = data[1:,2:]
        diam = data[0,2:]
    else:
        n_conc = data[1:,1:]
        diam = data[0,1:]
    time = data[1:,0]

    n_conc[n_conc<=0] = 1.01
    n_conc = log10(n_conc)
    n_conc[n_conc<=0] = 0.0

    return time, diam, n_conc

def loadNC(file):
    if '.sum' in file:
        return parseSum(file)
    try:
        import netCDF4
    except:
        return 'This feature need netCDF for Python'
    try:
        nc = netCDF4.Dataset(Path(file), 'r')
    except:
        return 'Could not open the file, is it accessible?'
    vnames = [i for i in nc.variables]
    for n in vnames:
        if 'number_concentration' in n.lower():
            n_conc_str = n
        if 'diameter' in n.lower():
            diam_str = n
        if 'radius' in n.lower():
            radius_str = n
        if 'time_in' in n.lower():
            time_str = n

    try:
        number_concentration = nc.variables[n_conc_str][:]
        time = nc.variables[time_str][:]
    except:
        return 'Could not find necessary variables from the file'
    try:
        diameter = nc.variables[diam_str][:]
    except:
        diameter = nc.variables[radius_str][:]*2

    try:
        number_concentration = lin_to_log(diameter[1,:],number_concentration*1e-6)
    except:
        return 'The file shape is not as it should be, is it complete run?'

    number_concentration[number_concentration<=0] = 1.01
    number_concentration = log10(number_concentration)
    number_concentration[number_concentration<=0] = 0.0


    return time, diameter[1,:], number_concentration
