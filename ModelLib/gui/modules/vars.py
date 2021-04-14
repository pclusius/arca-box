# -*- coding: utf-8 -*-
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


class INITFILE:
    def __init__(self, names=None):
        super(INITFILE,self).__init__()
        self.PATH = self._PATH()
        self.FLAG = self._FLAG()
        self.TIME = self._TIME()
        self.PARTICLE = self._PARTICLE()
        self.ENV = self._ENV()
        self.MCM = self._MCM()
        self.MODS = self._MODS(names)
        self.MISC = self._MISC()
        self.VAP = self._VAP()
        self.PRECISION = self._PRECISION()
        self.CUSTOM = self._CUSTOM()
        self.RAW = self._RAW()
        self.SETTINGS = self._SETTINGS()

    def printall(self, mods=None, target='p', f=None):
        if target == 'p':
            eol = ''
            cmd = 'print'
        if target == 'f':
            eol = '\\n'
            cmd = 'f.write'

        self.PATH.printall(cmd,f,eol)
        self.FLAG.printall(cmd,f,eol)
        self.TIME.printall(cmd,f,eol)
        self.PARTICLE.printall(cmd,f,eol)
        self.ENV.printall(cmd,f,eol)
        self.MCM.printall(cmd,f,eol)
        self.MODS.printall(cmd,f,eol,mods)
        self.MISC.printall(cmd,f,eol)
        self.VAP.printall(cmd,f,eol)
        self.PRECISION.printall(cmd,f,eol)
        self.CUSTOM.printall(cmd,f,eol)
        self.RAW.printall(cmd,f,eol)
        self.SETTINGS.printall(cmd,f,eol)

    class _PATH:
        def __init__(self):
            self.INOUT_DIR=0
            self.CASE_NAME=0
            self.RUN_NAME=0

        def printall(self,cmd,f,eol):
            exec("%s('&NML_PATH%s')"%(cmd, eol))
            exec("%s(' INOUT_DIR = \\'%s\\'%s')"%(cmd,self.INOUT_DIR,eol))
            exec("%s(' CASE_NAME = \\'%s\\'%s')"%(cmd,self.CASE_NAME,eol))
            exec("%s(' RUN_NAME = \\'%s\\'%s')"%(cmd,self.RUN_NAME,eol))
            exec("%s('/ \\n%s')"%(cmd, eol))

    class _FLAG:
        def __init__(self):
            self.CHEMISTRY_FLAG=0
            self.AEROSOL_FLAG=0
            self.ACDC_SOLVE_SS=0
            self.ACDC=0
            self.CONDENSATION=0
            self.COAGULATION=0
            self.DEPOSITION=0
            self.CHEM_DEPOSITION=0
            self.MODEL_H2SO4=0
            self.RESOLVE_BASE=0
            self.ORG_NUCL=0
            self.PRINT_ACDC=0
            self.USE_SPEED=0
            self.AFTER_CHEM_ON=0
            self.AFTER_NUCL_ON=0

        def printall(self,cmd,f,eol):
            exec("%s('&NML_FLAG%s')"%(cmd, eol))
            exec("%s(' CHEMISTRY_FLAG = %s%s')"%(cmd,self.CHEMISTRY_FLAG,eol))
            exec("%s(' AEROSOL_FLAG = %s%s')"%(cmd,self.AEROSOL_FLAG,eol))
            exec("%s(' ACDC_SOLVE_SS = %s%s')"%(cmd,self.ACDC_SOLVE_SS,eol))
            exec("%s(' ACDC = %s%s')"%(cmd,self.ACDC,eol))
            exec("%s(' CONDENSATION = %s%s')"%(cmd,self.CONDENSATION,eol))
            exec("%s(' COAGULATION = %s%s')"%(cmd,self.COAGULATION,eol))
            exec("%s(' DEPOSITION = %s%s')"%(cmd,self.DEPOSITION,eol))
            exec("%s(' CHEM_DEPOSITION = %s%s')"%(cmd,self.CHEM_DEPOSITION,eol))
            exec("%s(' MODEL_H2SO4 = %s%s')"%(cmd,self.MODEL_H2SO4,eol))
            exec("%s(' RESOLVE_BASE = %s%s')"%(cmd,self.RESOLVE_BASE,eol))
            exec("%s(' ORG_NUCL = %s%s')"%(cmd,self.ORG_NUCL,eol))
            exec("%s(' PRINT_ACDC = %s%s')"%(cmd,self.PRINT_ACDC,eol))
            exec("%s(' USE_SPEED = %s%s')"%(cmd,self.USE_SPEED,eol))
            exec("%s(' AFTER_CHEM_ON = %s%s')"%(cmd,self.AFTER_CHEM_ON,eol))
            exec("%s(' AFTER_NUCL_ON = %s%s')"%(cmd,self.AFTER_NUCL_ON,eol))
            exec("%s('/ \\n%s')"%(cmd, eol))

    class _TIME:
        def __init__(self):
            self.RUNTIME=0
            self.DT=0
            self.FSAVE_INTERVAL=0
            self.PRINT_INTERVAL=0
            self.FSAVE_DIVISION=0
            self.DATE='2000-01-01'
            self.INDEX=''

        def printall(self,cmd,f,eol):
            exec("%s('&NML_TIME%s')"%(cmd, eol))
            exec("%s(' RUNTIME = %s%s')"%(cmd,self.RUNTIME,eol))
            exec("%s(' DT = %s%s')"%(cmd,self.DT,eol))
            exec("%s(' FSAVE_INTERVAL = %s%s')"%(cmd,self.FSAVE_INTERVAL,eol))
            exec("%s(' PRINT_INTERVAL = %s%s')"%(cmd,self.PRINT_INTERVAL,eol))
            exec("%s(' FSAVE_DIVISION = %s%s')"%(cmd,self.FSAVE_DIVISION,eol))
            exec("%s(' DATE = \\'%s\\'%s')"%(cmd,self.DATE,eol))
            exec("%s(' INDEX = \\'%s\\'%s')"%(cmd,self.INDEX,eol))
            exec("%s('/ \\n%s')"%(cmd, eol))

    class _PARTICLE:
        def __init__(self):
            self.PSD_MODE=0
            self.N_BINS_PAR=0
            self.MIN_PARTICLE_DIAM=0
            self.MAX_PARTICLE_DIAM=0
            self.N_MODAL=0
            self.MMODAL_INPUT_INUSE=-1
            self.DMPS_FILE=0
            self.EXTRA_PARTICLES=0
            self.MMODAL_INPUT=0
            self.DMPS_READ_IN_TIME=0
            self.DMPS_HIGHBAND_LOWER_LIMIT=0
            self.DMPS_LOWBAND_UPPER_LIMIT=0
            self.USE_DMPS=0
            self.USE_DMPS_PARTIAL=0

        def printall(self,cmd,f,eol):
            exec("%s('&NML_PARTICLE%s')"%(cmd, eol))
            exec("%s(' PSD_MODE = %s%s')"%(cmd,self.PSD_MODE,eol))
            exec("%s(' N_BINS_PAR = %s%s')"%(cmd,self.N_BINS_PAR,eol))
            exec("%s(' MIN_PARTICLE_DIAM = %s%s')"%(cmd,self.MIN_PARTICLE_DIAM,eol))
            exec("%s(' MAX_PARTICLE_DIAM = %s%s')"%(cmd,self.MAX_PARTICLE_DIAM,eol))
            exec("%s(' N_MODAL = %s%s')"%(cmd,self.N_MODAL,eol))
            exec("%s(' MMODAL_INPUT_INUSE = %s%s')"%(cmd,self.MMODAL_INPUT_INUSE,eol))
            exec("%s(' DMPS_FILE = \\'%s\\'%s')"%(cmd,self.DMPS_FILE,eol))
            exec("%s(' EXTRA_PARTICLES = \\'%s\\'%s')"%(cmd,self.EXTRA_PARTICLES,eol))
            exec("%s(' MMODAL_INPUT = \\'%s\\'%s')"%(cmd,self.MMODAL_INPUT,eol))
            exec("%s(' DMPS_READ_IN_TIME = %s%s')"%(cmd,self.DMPS_READ_IN_TIME,eol))
            exec("%s(' DMPS_HIGHBAND_LOWER_LIMIT = %s%s')"%(cmd,self.DMPS_HIGHBAND_LOWER_LIMIT,eol))
            exec("%s(' DMPS_LOWBAND_UPPER_LIMIT = %s%s')"%(cmd,self.DMPS_LOWBAND_UPPER_LIMIT,eol))
            exec("%s(' USE_DMPS = %s%s')"%(cmd,self.USE_DMPS,eol))
            exec("%s(' USE_DMPS_PARTIAL = %s%s')"%(cmd,self.USE_DMPS_PARTIAL,eol))
            exec("%s('/ \\n%s')"%(cmd, eol))

    class _ENV:
        def __init__(self):
            self.ENV_FILE=0
            self.LOSSES_FILE=0
            self.CHAMBER_FLOOR_AREA=0
            self.CHAMBER_CHAMBER_HEIGHTAREA=0
            self.CHAMBER_HEIGHT=0
            self.EDDYK=0
            self.USTAR=0
            self.ALPHAWALL=0

        def printall(self,cmd,f,eol):
            exec("%s('&NML_ENV%s')"%(cmd, eol))
            exec("%s(' ENV_FILE = \\'%s\\'%s')"%(cmd,self.ENV_FILE,eol))
            exec("%s(' LOSSES_FILE = \\'%s\\'%s')"%(cmd,self.LOSSES_FILE,eol))
            exec("%s(' CHAMBER_FLOOR_AREA = %s%s')"%(cmd,self.CHAMBER_FLOOR_AREA,eol))
            exec("%s(' CHAMBER_CIRCUMFENCE = %s%s')"%(cmd,self.CHAMBER_CIRCUMFENCE,eol))
            exec("%s(' CHAMBER_HEIGHT = %s%s')"%(cmd,self.CHAMBER_HEIGHT,eol))
            exec("%s(' EDDYK = %s%s')"%(cmd,self.EDDYK,eol))
            exec("%s(' USTAR = %s%s')"%(cmd,self.USTAR,eol))
            exec("%s(' ALPHAWALL = %s%s')"%(cmd,self.ALPHAWALL,eol))
            exec("%s('/ \\n%s')"%(cmd, eol))

    class _MCM:
        def __init__(self):
            self.MCM_FILE=0

        def printall(self,cmd,f,eol):
            exec("%s('&NML_MCM%s')"%(cmd, eol))
            exec("%s(' MCM_FILE = \\'%s\\'%s')"%(cmd,self.MCM_FILE,eol))
            exec("%s('/ \\n%s')"%(cmd, eol))

    class _MODS:
        def __init__(self,names):
            self.names = names
            pass
        def printall(self,cmd,f,eol, mods=None):
            units =['C','K','Pa','hPa','bar','kPa','mbar','#','ppm','ppb','ppt','ppq']

            exec("%s('&NML_MODS%s')"%(cmd, eol))
            if mods != None:
                for v in self.names:
                    if v in mods:
                        m = mods[v]
                        unit = m.unit
                        if unit in units:
                            pass
                        else:
                            unit='#'
                        multistr = '%12.5e'%(m.multi)
                        multistr = multistr.replace('e', 'd', 1)
                        shiftstr = '%12.5e'%(m.shift)
                        shiftstr = shiftstr.replace('e', 'd', 1)
                        minstr = '%12.5e'%(m.min)
                        minstr = minstr.replace('e', 'd', 1)
                        maxstr = '%12.5e'%(m.max)
                        maxstr = maxstr.replace('e', 'd', 1)
                        if m.pmInUse == 'Yes' or m.pmInUse == 'yes':
                            mode = max(1,abs(m.mode))
                        else: mode = abs(m.mode) * -1
                        if 'str' in str(type(m.col)) :
                            m.col = -1
                        strr = "MODS(%d)%s= %2d %3d %s %s %s %s %fd0 %0fd0 %fd0 %fd0 %fd0 %s%s%s %s%s%s ! %s"%(
                        m.Find,' '*(4-len(str(m.Find))),mode,m.col, multistr,shiftstr,minstr, maxstr, m.sig,m.mju, m.fv,m.ph,m.am, "\\'", unit,"\\'", "\\'", m.tied,"\\'", m.name)
                        exec("%s(' %s%s')"%(cmd,strr,eol))
            exec("%s('/ \\n%s')"%(cmd, eol))

    class _MISC:
        def __init__(self):
            self.LAT=0
            self.LON=0
            self.WAIT_FOR=0
            self.DESCRIPTION=0
            self.CH_ALBEDO=0
            self.DMA_F=0
            self.RESOLVE_BASE_PRECISION=0
            self.FILL_FORMATION_WITH=0
            self.SKIP_ACDC=0
            self.GR_SIZES=0

        def printall(self,cmd,f,eol):
            exec("%s('&NML_MISC%s')"%(cmd, eol))
            exec("%s(' LAT = %s%s')"%(cmd,self.LAT,eol))
            exec("%s(' LON = %s%s')"%(cmd,self.LON,eol))
            exec("%s(' WAIT_FOR = %s%s')"%(cmd,self.WAIT_FOR,eol))
            exec("%s(' DESCRIPTION = \\'%s\\'%s')"%(cmd,self.DESCRIPTION,eol))
            exec("%s(' CH_ALBEDO = %s%s')"%(cmd,self.CH_ALBEDO,eol))
            exec("%s(' DMA_F = %s%s')"%(cmd,self.DMA_F,eol))
            exec("%s(' RESOLVE_BASE_PRECISION = %s%s')"%(cmd,self.RESOLVE_BASE_PRECISION,eol))
            exec("%s(' FILL_FORMATION_WITH = \\'%s\\'%s')"%(cmd,self.FILL_FORMATION_WITH,eol))
            exec("%s(' SKIP_ACDC = %s%s')"%(cmd,self.SKIP_ACDC,eol))
            exec("%s(' GR_SIZES = \\'%s\\'%s')"%(cmd,self.GR_SIZES,eol))
            exec("%s('/ \\n%s')"%(cmd, eol))

    class _VAP:
        def __init__(self):
            self.USE_ATOMS=0
            self.VAP_NAMES=0
            self.VAP_ATOMS=0

        def printall(self,cmd,f,eol):
            exec("%s('&NML_VAP%s')"%(cmd, eol))
            exec("%s(' USE_ATOMS = %s%s')"%(cmd,self.USE_ATOMS,eol))
            exec("%s(' VAP_NAMES = \\'%s\\'%s')"%(cmd,self.VAP_NAMES,eol))
            exec("%s(' VAP_ATOMS = \\'%s\\'%s')"%(cmd,self.VAP_ATOMS,eol))
            exec("%s('/ \\n%s')"%(cmd, eol))

    class _PRECISION:
        def __init__(self):
            self.DDIAM_RANGE=''
            self.DPNUM_RANGE=''
            self.DVAPO_RANGE=''

        def printall(self,cmd,f,eol):
            exec("%s('&NML_PRECISION%s')"%(cmd, eol))
            exec("%s(' DDIAM_RANGE = %s%s')"%(cmd,self.DDIAM_RANGE,eol))
            exec("%s(' DPNUM_RANGE = %s%s')"%(cmd,self.DPNUM_RANGE,eol))
            exec("%s(' DVAPO_RANGE = %s%s')"%(cmd,self.DVAPO_RANGE,eol))
            exec("%s('/ \\n%s')"%(cmd, eol))

    class _CUSTOM:
        def __init__(self):
            self.CUSTOMS=[]

        def printall(self,cmd,f,eol):
            exec("%s('&NML_CUSTOM%s')"%(cmd, eol))
            for key, value in self.CUSTOMS:
                if key.strip() != '':
                    try:
                        float(value)
                        exec("%s(' %s = %s%s')"%(cmd,key,value,eol))
                    except:
                        if value.upper() == '.TRUE.' or value.upper() == '.FALSE.':
                            exec("%s(' %s = %s%s')"%(cmd,key,value,eol))
                        else:
                            try:
                                _ = [float(iii) for iii in value.split(',')]
                                exec("%s(' %s = %s%s')"%(cmd,key,value,eol))
                            except:
                                exec("%s(' %s = \\'%s\\'%s')"%(cmd,key,value.replace("'",'').replace('"',''),eol))
            exec("%s('/ \\n%s')"%(cmd, eol))

    class _RAW:
        def __init__(self):
            self.RAW=''

        def printall(self,cmd,f,eol):
            exec("%s('%s%s')"%(cmd,self.RAW.replace('\n','\\n').replace("'", "\\'").replace('"', '\\"'),eol))
            exec("%s('\\n%s')"%(cmd, eol))
            exec("%s('# Following settings are for the GUI and not directly used by the model ----- %s')"%(cmd, eol))
            exec("%s('# RAW_INPUT = %s%s')"%(cmd,self.RAW.replace('\n','<br>').replace("'", "\\'").replace('"', '\\"'),eol))


    class _SETTINGS:
        def __init__(self):
            self.BATCH=0
            self.INPUT=0

        def printall(self,cmd,f,eol):
            exec("%s('# INPUT_SETTINGS = \\'%s\\'%s')"%(cmd,self.INPUT,eol))
            exec("%s('# BATCH_SETTINGS = \\'%s\\'%s')"%(cmd,self.BATCH,eol))
            exec("%s('#---------------------------------------------------------------------------- \\n%s')"%(cmd, eol))


mods = {}
