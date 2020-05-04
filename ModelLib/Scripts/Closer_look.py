import numpy as np
import matplotlib.pyplot as plt
import netCDF4
from pathlib import Path


def comp(file, time, bin):
    nc = netCDF4.Dataset(Path(file), 'r')
    composition = nc.variables['PARTICLE_COMPOSITION'][:]
    nc.close()
    focus = composition[time,bin,:]
    return focus

def part(file, time):
    nc = netCDF4.Dataset(Path(file), 'r')
    particles = nc.variables['NUMBER_CONCENTRATION'][:]
    nc.close()
    particles = particles[time,:]
    return particles


# file = '/home/pecl/05-APCAD/supermodel-phase-1/INOUT/HYDE/PC_2018-04-10/AREF/Particle.nc'
# print('ref vika', part(file, -1)[:12])
# # print('ref vika', comp(file, -1, 0)[-12:])
# file = '/home/pecl/05-APCAD/supermodel-phase-1/INOUT/HYDE/PC_2018-04-10/ACONT/Particle.nc'
# # print('cont eka', comp(file, 0, 0)[-12:])
# print('cont eka', part(file, 0)[:12])
