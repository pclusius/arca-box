# -*- coding: utf-8 -*-
from os.path import split as ossplit
from numpy import sum as npsum
from numpy import ravel, array
import netCDF4
import numpy.ma as ma
from modules.head import units

class NcPlot:
    """Class for plot file contents"""
    def __init__(self, file):
        self.path = file
        self.masterfile = ossplit(file)[1]
        self.legend = ossplit(ossplit(file)[0])[1]
        ncs = netCDF4.Dataset(file, 'r')
        self.getaircc(file, ncs)
        self.parvars = {}
        self.convars = {}
        self.invvars = {}
        self.par = False
        if self.masterfile == 'Particles.nc':
            self.par = True
            b = ravel(ncs.variables['CONDENSABLES'][:,:].astype(str),'C')
            names=b.reshape(ncs.variables['CONDENSABLES'].shape)
            m3 = ncs.variables['NUMBER_CONCENTRATION'][:]*1e18 # to convert to m3 and nanograms/m3
            self.composition_ng = npsum(ncs.variables['PARTICLE_COMPOSITION'][:,:,:]*m3[:,:,newaxis], 1)
            for i,word in enumerate(names):
                self.parvars[(''.join(list(word))).strip()] = i

        for timedim in ncs.dimensions:
            if ncs.dimensions[timedim].isunlimited():
                break

        checker = lambda v,n: v.lower() in timedim and 'Shifter' not in n and 'Multipl' not in n and 'TIME_IN' not in n.upper()
        cache = array([i.name for i in ncs.get_variables_by_attributes(ndim=1)])
        timevars = [checker(i.dimensions[0], i.name) for i in ncs.get_variables_by_attributes(ndim=1)]
        self.varnames = cache[timevars]

        # Unfortunately these early version files are still somewhere out there
        try:
            self.time = ncs.variables['TIME_IN_SEC'][:]/3600
        except:
            try:
                self.time = ncs.variables['time_in_sec'][:]/3600
            except:
                self.time = ncs.variables['Time_in_sec'][:]/3600

        if ma.is_masked(self.time):
            self.is_masked = True
            self.mask = ~self.time.mask
        else:
            self.is_masked = False
            self.mask = self.time == self.time

        self.time = self.time[self.mask]
        # self.conc_matrix = zeros((len(self.time),len(self.varnames[self.mask])))
        for i,n in enumerate(cache[timevars]):
            self.convars[n] = i
            self.invvars[i] = n
            # self.conc_matrix[:,i] = ncs.variables[n][self.mask]
        self.nc = ncs

    def close(self):
        self.nc.close()

    def getconc(self,n):
        if n in self.convars:
            return self.nc.variables[n][self.mask], units.get(n,units['REST'])[0]
        else: return

    def getloc(self,i):
        if i in self.invvars:
            return self.nc.variables[self.invvars[i]][self.mask], units.get(self.invvars[i],units['REST'])[0]
        else: return

    def getcom(self,n):
        if self.par and n in self.parvars:
            return self.composition_ng[:,self.parvars[n]][self.mask], units.get(n,units['REST'])[0]
        else: return

    def getcomsum(self,names):
        retarr = zeros(len(self.mask))
        for i,n in enumerate(names):
            if i==0: u = units.get(n,units['REST'])[0]
            if self.par and n in self.parvars:
                if u == units.get(n,units['REST'])[0]:
                    unit = True
                else:
                    unit = False
                retarr += self.composition_ng[:,self.parvars[n]][self.mask]
        if not unit: u='[-]'
        return retarr, u

    def getconcsum(self,names):
        retarr = zeros(len(self.mask))
        for i,n in enumerate(names):
            if i==0: u = units.get(n,units['REST'])[0]
            if n in self.convars:
                if u == units.get(n,units['REST'])[0]:
                    unit = True
                else:
                    unit = False
                retarr += self.getconc(n)[0]
        if not unit: u='[-]'
        return retarr, u

    def getaircc(self, file, ncs):
        try:
            if self.masterfile != 'General.nc':
                air_nc = netCDF4.Dataset(osjoin(ossplit(file)[0],'General.nc'), 'r')
                temp = air_nc.variables['TEMPK'][:]
                pres = air_nc.variables['PRESSURE'][:]
                air_nc.close()
            else:
                temp = ncs.variables['TEMPK'][:]
                pres = ncs.variables['PRESSURE'][:]
            self.aircc = 1e-6 * pres / temp /1.38064852e-23
            self.have_aircc = True
        except:
            # qt_box.popup('Missing air concentration', 'File "General.nc" was not found from this directory. Cannot calculate mixing ratios but will show values.')
            self.have_aircc = False
