#

# class Comp:
#     def __init__(self):
#         self.index  = 0
#         self.mode  = 0
#         self.col  = -1
#         self.multi = 1e0    # Multiplication factor in MODE0
#         self.shift = 0e0    # Constant to be addded in MODE0
#         self.min = 1e1      # Minimum value for the parametrized concentration OR constant value if max <= min
#         self.max = 1e5      # Peak value
#         self.sig = 1e0      # Standard deviation for the Gaussian=sig of the bell curve
#         self.mju = 12e0     # Time of peak value
#         self.fv  = 0e0      # Angular frequency [hours] of modifying sine function
#         self.ph  = 0e0      # Angular frequency [hours] of modifying sine function
#         self.am  = 1e0      # Amplitude of modificaion
#         self.name = 'NONAME'# Human readable name for modified variable
#         self.unit = '[-]'     # unit name
#         self.Find = 1
#         self.pmInUse = 'no'


class INITFILE:
    def __init__(self, names=None):
        super(INITFILE,self).__init__()
        # self.names = names
        self.PATH = self._PATH()
        self.FLAG = self._FLAG()
        self.TIME = self._TIME()
        self.PARTICLE = self._PARTICLE()
        self.ENV = self._ENV()
        self.MCM = self._MCM()
        self.MODS = self._MODS(names)
        self.MISC = self._MISC()
        self.VAP = self._VAP()

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

    class _PATH:
        def __init__(self):
            # self.WORK_DIR=0
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
            self.NUCLEATION=0
            self.ACDC=0
            self.EXTRA_DATA=0
            self.CURRENT_CASE=0
            self.CONDENSATION=0
            self.COAGULATION=0
            self.MODEL_H2SO4=0
            self.RESOLVE_BASE=0
            self.PRINT_ACDC=0

        def printall(self,cmd,f,eol):
            exec("%s('&NML_FLAG%s')"%(cmd, eol))
            exec("%s(' CHEMISTRY_FLAG = %s%s')"%(cmd,self.CHEMISTRY_FLAG,eol))
            exec("%s(' AEROSOL_FLAG = %s%s')"%(cmd,self.AEROSOL_FLAG,eol))
            exec("%s(' ACDC_SOLVE_SS = %s%s')"%(cmd,self.ACDC_SOLVE_SS,eol))
            exec("%s(' NUCLEATION = %s%s')"%(cmd,self.NUCLEATION,eol))
            exec("%s(' ACDC = %s%s')"%(cmd,self.ACDC,eol))
            exec("%s(' EXTRA_DATA = %s%s')"%(cmd,self.EXTRA_DATA,eol))
            exec("%s(' CURRENT_CASE = %s%s')"%(cmd,self.CURRENT_CASE,eol))
            exec("%s(' CONDENSATION = %s%s')"%(cmd,self.CONDENSATION,eol))
            exec("%s(' COAGULATION = %s%s')"%(cmd,self.COAGULATION,eol))
            exec("%s(' MODEL_H2SO4 = %s%s')"%(cmd,self.MODEL_H2SO4,eol))
            exec("%s(' RESOLVE_BASE = %s%s')"%(cmd,self.RESOLVE_BASE,eol))
            exec("%s(' PRINT_ACDC = %s%s')"%(cmd,self.PRINT_ACDC,eol))
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
            self.N_BINS_PARTICLE=0
            self.MIN_PARTICLE_DIAM=0
            self.MAX_PARTICLE_DIAM=0
            # self.DMPS_DIR=0
            # self.EXTRA_P_DIR=0
            self.DMPS_FILE=0
            self.EXTRA_PARTICLES=0
            self.DMPS_READ_IN_TIME=0
            self.DMPS_HIGHBAND_LOWER_LIMIT=0
            self.DMPS_LOWBAND_UPPER_LIMIT=0
            self.USE_DMPS=0
            self.USE_DMPS_SPECIAL=0

        def printall(self,cmd,f,eol):
            exec("%s('&NML_PARTICLE%s')"%(cmd, eol))
            exec("%s(' PSD_MODE = %s%s')"%(cmd,self.PSD_MODE,eol))
            exec("%s(' N_BINS_PARTICLE = %s%s')"%(cmd,self.N_BINS_PARTICLE,eol))
            exec("%s(' MIN_PARTICLE_DIAM = %s%s')"%(cmd,self.MIN_PARTICLE_DIAM,eol))
            exec("%s(' MAX_PARTICLE_DIAM = %s%s')"%(cmd,self.MAX_PARTICLE_DIAM,eol))
            # exec("%s(' DMPS_DIR = \\'%s\\'%s')"%(cmd,self.DMPS_DIR,eol))
            # exec("%s(' EXTRA_P_DIR = \\'%s\\'%s')"%(cmd,self.EXTRA_P_DIR,eol))
            exec("%s(' DMPS_FILE = \\'%s\\'%s')"%(cmd,self.DMPS_FILE,eol))
            exec("%s(' EXTRA_PARTICLES = \\'%s\\'%s')"%(cmd,self.EXTRA_PARTICLES,eol))
            exec("%s(' DMPS_READ_IN_TIME = %s%s')"%(cmd,self.DMPS_READ_IN_TIME,eol))
            exec("%s(' DMPS_HIGHBAND_LOWER_LIMIT = %s%s')"%(cmd,self.DMPS_HIGHBAND_LOWER_LIMIT,eol))
            exec("%s(' DMPS_LOWBAND_UPPER_LIMIT = %s%s')"%(cmd,self.DMPS_LOWBAND_UPPER_LIMIT,eol))
            exec("%s(' USE_DMPS = %s%s')"%(cmd,self.USE_DMPS,eol))
            exec("%s(' USE_DMPS_SPECIAL = %s%s')"%(cmd,self.USE_DMPS_SPECIAL,eol))
            exec("%s('/ \\n%s')"%(cmd, eol))

    class _ENV:
        def __init__(self):
            # self.ENV_PATH=0
            self.ENV_FILE=0
            self.TEMPUNIT='C'

        def printall(self,cmd,f,eol):
            exec("%s('&NML_ENV%s')"%(cmd, eol))
            # exec("%s(' ENV_PATH = \\'%s\\'%s')"%(cmd,self.ENV_PATH,eol))
            exec("%s(' ENV_FILE = \\'%s\\'%s')"%(cmd,self.ENV_FILE,eol))
            exec("%s(' TEMPUNIT = \\'%s\\'%s')"%(cmd,self.TEMPUNIT,eol))
            exec("%s('/ \\n%s')"%(cmd, eol))

    class _MCM:
        def __init__(self):
            # self.MCM_PATH=0
            self.MCM_FILE=0

        def printall(self,cmd,f,eol):
            exec("%s('&NML_MCM%s')"%(cmd, eol))
            # exec("%s(' MCM_PATH = \\'%s\\'%s')"%(cmd,self.MCM_PATH,eol))
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
                            mode = m.mode
                        else: mode = 0
                        strr = "MODS(%d)%s= %d %3d %s %s %s %s %fd0 %0fd0 %fd0 %fd0 %fd0 %s%s%s ! %s"%(
                        m.Find,' '*(4-len(str(m.Find))),mode,m.col, multistr,shiftstr,minstr, maxstr, m.sig,m.mju, m.fv,m.ph,m.am, "\\'", unit,"\\'", m.name)
                        exec("%s(' %s%s')"%(cmd,strr,eol))
            exec("%s('/ \\n%s')"%(cmd, eol))

    class _MISC:
        def __init__(self):
            # self.JD=0
            self.LAT=0
            self.LON=0
            self.WAIT_FOR=0
            self.PYTHON=0
            self.DESCRIPTION=0
            # self.SOLVER=0
            self.CH_ALBEDO=0
            self.DMA_F=0
            self.RESOLVE_BASE_PRECISION=0
            self.FILL_FORMATION_WITH=0

        def printall(self,cmd,f,eol):
            exec("%s('&NML_MISC%s')"%(cmd, eol))
            # exec("%s(' JD = %s%s')"%(cmd,self.JD,eol))
            exec("%s(' LAT = %s%s')"%(cmd,self.LAT,eol))
            exec("%s(' LON = %s%s')"%(cmd,self.LON,eol))
            exec("%s(' WAIT_FOR = %s%s')"%(cmd,self.WAIT_FOR,eol))
            exec("%s(' PYTHON = %s%s')"%(cmd,self.PYTHON,eol))
            exec("%s(' DESCRIPTION = \\'%s\\'%s')"%(cmd,self.DESCRIPTION,eol))
            # exec("%s(' SOLVER = \\'%s\\'%s')"%(cmd,self.SOLVER,eol))
            exec("%s(' CH_ALBEDO = %s%s')"%(cmd,self.CH_ALBEDO,eol))
            exec("%s(' DMA_F = %s%s')"%(cmd,self.DMA_F,eol))
            exec("%s(' RESOLVE_BASE_PRECISION = %s%s')"%(cmd,self.RESOLVE_BASE_PRECISION,eol))
            exec("%s(' FILL_FORMATION_WITH = \\'%s\\'%s')"%(cmd,self.FILL_FORMATION_WITH,eol))
            exec("%s('/ \\n%s')"%(cmd, eol))

    class _VAP:
        def __init__(self):
            self.VAP_LOGICAL=0
            self.VAP_NAMES=0
            self.VAP_PROPS=0

        def printall(self,cmd,f,eol):
            exec("%s('&NML_VAP%s')"%(cmd, eol))
            exec("%s(' VAP_LOGICAL = %s%s')"%(cmd,self.VAP_LOGICAL,eol))
            exec("%s(' VAP_NAMES = \\'%s\\'%s')"%(cmd,self.VAP_NAMES,eol))
            exec("%s(' VAP_PROPS = \\'%s\\'%s')"%(cmd,self.VAP_PROPS,eol))
            exec("%s('/ \\n%s')"%(cmd, eol))


mods = {}
