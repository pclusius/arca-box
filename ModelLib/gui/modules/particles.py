"""
=============================================================================
Copyright (C) 2021  Multi-Scale Modelling group
Institute for Atmospheric and Earth System Research (INAR), University of Helsinki
Contact information arca@helsinki.fi

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
=============================================================================
"""


from numpy import log10, genfromtxt,nan_to_num
from pathlib import Path
from os import stat
from os.path import exists

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


def parseSum(file, assume_log):
    lin = False
    if not exists(file):
        return 'File not found'
    if stat(file).st_size == 0: return 'File '+file+' is empty.'
    try:
        data = genfromtxt(Path(file),invalid_raise=False)
        if len(data.shape)<2: return 'Data in file '+file+' has not enough dimensions'
    except IOError:
        return 'Could not open the file: '+data
    if data[0,1] == 0:
        n_conc = data[1:,2:]
        diam = data[0,2:]
    elif assume_log:
        n_conc = data[1:,1:]
        diam = data[0,1:]
        print('Assuming lognormalized input (menu: Plot->Assume lognormal input)')
    else:
        lin = True
        n_conc = data[1:,1:]
        diam = data[0,1:]
        print('Assuming raw input')
    time = data[1:,0]
    nan_to_num(n_conc,copy=False)
    if lin:
        n_conc = n_conc / log10(diam[3]/diam[2])
    n_conc[n_conc<=0] = 1.01
    n_conc = log10(n_conc)
    n_conc[n_conc<=0] = 0.0
    return time-time[0], diam, n_conc

def loadNC(file, passon=False):
    if ('.sum' in file) or ('.dat' in file):
        return parseSum(file, passon)
    try:
        import netCDF4
    except:
        return 'This feature needs netCDF4 for Python'
    try:
        nc = netCDF4.Dataset(Path(file), 'r')
    except:
        return 'Could not open the file, is it accessible?'
    vnames = [i for i in nc.variables]
    for n in vnames:
        if 'number_concentration' == n.lower():
            n_conc_str = n
        if 'diameter' == n.lower():
            diam_str = n
        if 'radius' == n.lower():
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
    nc.close()
    try:
        number_concentration = number_concentration/log10(diameter[1,3]/diameter[1,2])
    except:
        return 'The file shape is not as it should be, is it complete run?'

    number_concentration[number_concentration<=0] = 1.01
    number_concentration = log10(number_concentration)
    number_concentration[number_concentration<=0] = 0.0


    return time, diameter[1,:], number_concentration
