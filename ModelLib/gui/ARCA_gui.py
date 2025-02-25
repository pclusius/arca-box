#!/usr/bin/env python3
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

from PyQt5 import QtCore, QtWidgets, QtGui, uic
import pyqtgraph as pg
from layouts import varWin,gui10,batchDialog1,batchDialog2,batchDialog3,vdialog,cc,about,input,t_editor
from modules import variations,vars,batch,GetVapourPressures as gvp
from modules.grepunit import grepunit
from subprocess import Popen, PIPE, STDOUT
from numpy import linspace,log10,sqrt,log,exp,pi,sin,shape,unique,array,ndarray,where,newaxis,flip,zeros, sum as npsum, ravel, mean, round as npround
from numpy import argsort
import numpy.ma as ma
from os import walk, mkdir, getcwd, chdir, chmod, environ, system, name as osname, remove as osremove, rename as osrename
from os.path import exists, dirname, getmtime, abspath, split as ossplit, join as osjoin, relpath as osrelpath
from shutil import copyfile as cpf
from shutil import move as mvf
from re import sub,IGNORECASE, findall, match, finditer
import time
import re
import sys,os
import pickle
from urllib.parse import urljoin
import csv
from modules.config import get_config, set_config, remove_config

try:
    from scipy.ndimage import gaussian_filter
    from scipy.signal import savgol_filter
    scipyIs = True
except:
    print('Consider adding SciPy to your Python')
    scipyIs = False
try:
    from modules import particles as par
    import netCDF4
    netcdf = True
except:
    print('NetCDF4 for Python is essential for full functionality.')
    netcdf = False

try:
    import platform
    operatingsystem = platform.system()
    # "Windows"/"Linux"/"Darwin"
except:
    operatingsystem = 'Linux'

# -----------------------------------------------------------------------------
# Generally these default settings should not changed, do so with your own risk
# -----------------------------------------------------------------------------

environ["QT_AUTO_SCREEN_SCALE_FACTOR"] = "1"
NAMES_override = ''
QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling, True) # enable highdpi scaling
QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_UseHighDpiPixmaps, True)    # use highdpi icons
# See if scaling is necessary, currently only on Windows
args = []
startupInitfile = None
# Get the dark/light mode settings
_dark_mode = get_config("style", "darkmode", fallback='False')
if _dark_mode.lower() == 'false' or _dark_mode.lower() == 'true':
    dark_mode = eval('%s'%(_dark_mode.capitalize()))
    set_config("style", "darkmode", str(dark_mode))
if dark_mode:
    icondir = 'dark'
    darken_by = float(get_config("style", "darkness", fallback=0.35))
    set_config("style", "darkness", str(darken_by))

else:
    icondir = 'light'

if len(sys.argv)>1:
    for ia,a in enumerate(sys.argv):
        if '--scaling_' in a:
            environ["QT_SCALE_FACTOR"] = "%3.2f"%float(a.replace('--scaling_',''))
            args.append('-NS')
        elif a == '--names':
            NAMES_override = sys.argv[ia+1]
            print(f'WARNING: using list of available input from {NAMES_override}')
        elif a == '--initfile':
            startupInitfile = sys.argv[ia+1]
            print(f'Loading settings from: {startupInitfile}')
        elif '--scaleall_' in a:
            sfs = a.replace('--scaleall_','').split('_')
            environ["QT_SCREEN_SCALE_FACTORS"] = "%s"%(';'.join(sfs))
            args.append('-NS')

if osname.upper() == 'NT' or operatingsystem == 'Windows' and not '-NS' in args: # this or is only till I sort out if platform usually works for people
    try:
        import ctypes
        sf = (ctypes.windll.shcore.GetScaleFactorForDevice(0) / 100)
        if sf<1.3:
            sf = 1
        sf = 1/sf
        u32 = ctypes.windll.user32
        scrhgt = u32.GetSystemMetrics(1)
        if scrhgt < 850:
            sf = sf * scrhgt / 850.
    except:
        sf=1
        print("Could not get the scaling factor of the screen, using %3.2f. Let's hope the GUI looks ok"%sf)
    environ["QT_SCALE_FACTOR"] = "%3.2f"%sf

if '-NS' in args: print('Overriding scaling from cmdline')

## Some constants --------------------------------------------
# widths of the columns in "Input variables" tab
column_widths = [150,80,70,70,55,80,50,3]

# available units for variables, used to fill the tables and graphs with appropriate units
units = {
'TEMPK': ['K','C'],
'TIME_IN_HRS': ['hrs'],
'TIME_IN_SEC': ['s'],
'PRESSURE': ['Pa','hPa','mbar','kPa','bar','atm'],
'REL_HUMIDITY': ['%'],
'CONDENS_SINK':['s⁻¹'],
'CS_CALC':['s⁻¹'],
'CON_SIN_NITR':['s⁻¹'],
'SW_RADIATION':['W m⁻²'],
'ION_PROD_RATE':['ip cm⁻³s⁻¹'],
'NUC_RATE_IN':['cm⁻³s⁻¹'],
'J_ACDC_1_CM3':['cm⁻³s⁻¹'],
'J_ACDC_2_CM3':['cm⁻³s⁻¹'],
'J_ACDC_3_CM3':['cm⁻³s⁻¹'],
'J_ACDC_4_CM3':['cm⁻³s⁻¹'],
'J_ACDC_5_CM3':['cm⁻³s⁻¹'],
'J_ACDC_NH3_CM3':['cm⁻³s⁻¹'],
'J_ACDC_DMA_CM3':['cm⁻³s⁻¹'],
'J_ACDC_SUM_CM3':['cm⁻³s⁻¹'],
'J_TOTAL_CM3':['cm⁻³s⁻¹'],
'emi':['cm⁻³s⁻¹'],
'GMD':['m'],
'GSD':['[-]'],
'REST':['cm⁻³','ppm','ppb','ppt','ppq']
}
# Name of the executable -------------------------------------------
exe_name = 'arcabox.exe'
# Path to variable names -------------------------------------------
if NAMES_override != '':
    path_to_names = NAMES_override
else:
    path_to_names = 'ModelLib/required/NAMES.dat'
# Path to extra variables ------------------------------------------
path_to_xtras = 'ModelLib/required/AEMS.dat'
# GUI root
gui_path = 'ModelLib/gui/'
# GUI defaults are saved into this file. If it doesn't exist, it gets created in first start
defaults_file_path = gui_path+'conf/defaults.init'
minimal_settings_path = gui_path+'conf/minimal.init'
monitorfont_pickle = 'conf/monitorfont.pickle'
globalfont_pickle = 'conf/globalfont.pickle'

# Default global shortwave radiation spectrum file
# (located in ModelLib/Photolyse/Spectra) -------------------------
defaultSpectrum = 'glob_swr_distr.txt'
# This path will be added to Common out if no other option is given
default_inout = 'INOUT'
# This case name will be used as case name if no other name is given
default_case = 'DEFAULTCASE'.upper()
# This run name will be used as run name if no other name is given
default_run  = 'DEFAULTRUN'.upper()
# name and location of the temporary settings file used for test runs
tempfile = gui_path+'tmp/GUI_INIT.tmp'
# initial maximum for function creator tab sliders
slMxs = (200,190,220,100,200)
# 10 colors for plots
colors = [(120,0,0),(180,0,0),(220,0,0),(255,10,0),(255,85,0),
(255,85,140),(255,85,255),(180,0,255),(110,0,255),(0,0,255)]
# linetypes
linetypes = [QtCore.Qt.SolidLine, QtCore.Qt.DashLine]
# BG colour for ENV vars
env_no  = array((215,238,244))
env_yes = array((128,179,255))
# BG colour for ORG vars
org_no  = array((243,240,239))
org_yes = array((172,147,147))

# BG colour for XTR vars
xtr_no  = array((255, 230, 128))
xtr_yes = array((212, 170, 0))

if dark_mode:
    org_yes = array(darken_by*org_yes,dtype=int)
    org_no = array(darken_by*org_no,dtype=int)
    env_yes = array(darken_by*env_yes,dtype=int)
    env_no = array(darken_by*env_no,dtype=int)
    xtr_yes = array(darken_by*xtr_yes,dtype=int)
    xtr_no = array(darken_by*xtr_no,dtype=int)

# icon
modellogo = gui_path+"/icons/%s/ArcaLogo.png"%icondir
boxicon = gui_path+"/icons/%s/thebox_ico.png"%icondir
CurrentVersion = "ARCA Box Model "
# Some messages
netcdfMissinnMes = ('Please note:',
'To open NetCDF-files you need netCDF4 for Python.\nYou can istall it with pip, package manager (or perhaps: python3 -m pip install --user netCDF4.')

# get current directory (to render relative paths) ----------
currentdir   = getcwd()
currentdir   = currentdir.replace('/ModelLib/gui', '')
currentdir_l = len(currentdir)
chdir(currentdir)

with open("ModelLib/required/version.txt", encoding="utf-8") as f:
    for line in f: break
CurrentVersion += line.strip('\n\r').strip('\n')

helpd = {}
with open(osjoin(gui_path,'conf','helplinks.txt'), 'r', encoding="utf-8") as b:
    for l in b:
        ll = l.replace('\n','').replace('\r\n','')
        if ll != '':
            k,v = ll.split(',')
            if not 'http://' in v and not 'https://' in v:
                v = urljoin('file:///', osjoin(getcwd(),v))
            helpd[k] = v.strip('\n')


# files that can be modified with the Editor
nucl_homs = "ModelLib/required/nucl_homs.txt"
custom_functions = "src/custom_functions.f90"
SCREENPRINT_NML = "ModelLib/required/NML_SCREENPRINTS.def"
# Create chemistry script location
ccloc = 'ModelLib/gui/chemistry_package_PZ'

# Guess initially the current python version and check the calling script for precise version
currentPythonVer = 'python'
try:
    with open('run_arca.sh') as f:
        for line in f:
            if gui_path in line:
                currentPythonVer = line.split()[0]
except:
    pass

NAMES = []
namesPyInds = {}
namesFoInds = {}
bold = QtGui.QFont()
bold.setBold(True)
bold.setWeight(75)
roman = QtGui.QFont()
roman.setBold(False)
roman.setWeight(50)


## Create lists and dictionaries related to NAMES -------------
i = 0
j = 0
with open(path_to_names) as f:
    for line in f:
        j += 1
        name = line[:-1]
        if name == '':
            continue
        if i+1!=j:
            print('\n\nEmpty line in the middle of NAMES.DAT, it will not work!\n\n')
            exit('ARCA will not start, check out '+path_to_names+' and remove empty lines.')
        if '#' in line:
            divider_i=i
            name = 'MCM compounds start here'
        NAMES.append(name)
        namesPyInds[name] = i
        namesFoInds[name] = i+1
        i += 1

divider_xtr_i = i
with open(path_to_xtras) as f:
    for line in f:
        j += 1
        name = line.strip()
        if name == '':
            print('\n\nEmpty line in the middle of AEMS.dat, it will not work!\n\n')
            exit('ARCA will not start, check out '+path_to_xtras+' and remove empty lines.')
        NAMES.append(name)
        namesPyInds[name] = i
        namesFoInds[name] = i+1
        i += 1

## -----------------------------------------------------------
nml = vars.INITFILE(NAMES)
mmc = {'True': 'k', 'False': 0.95}
# The popup window for batch file preview
class batchW(QtWidgets.QDialog):
    def __init__(self, parent = None, n=0):
        super(batchW, self).__init__(parent)
        if n==3:
            self.ui = batchDialog3.Ui_batchDialog()
        elif n==2:
            self.ui = batchDialog2.Ui_batchDialog()
        else:
            self.ui = batchDialog1.Ui_batchDialog()
        self.ui.setupUi(self)
        self.ui.bDialogbuttonBox.accepted.connect(self.accept)
        self.ui.bDialogbuttonBox.rejected.connect(self.reject)
        self.ui.bDialogbuttonBox.helpRequested.connect(lambda: qt_box.helplink('batch'))

    def settext(self,a):
        """Setter for window text"""
        c = 1
        for i in range(3):
            if a[0][i]>0:
                exec('self.ui.label_%d.setText(a[1][%d])'%(c,i))
                exec('self.ui.tb_%d.appendPlainText(\'\'.join(a[2][%d]))'%(c,i))
                c +=1

# The popup window for About ARCA
class About(QtWidgets.QDialog):
    def __init__(self, parent = None):
        super(About, self).__init__(parent)
        self.ab = about.Ui_Dialog()
        self.ab.setupUi(self)
        self.ab.okgreat.clicked.connect(self.reject)
        self.ab.logo.setPixmap(QtGui.QPixmap(modellogo.replace('.png', 'HR.png')))

# The popup window for Create KPP files
class CCWin(QtWidgets.QDialog):
    def __init__(self, parent = None):
        super(CCWin, self).__init__(parent)
        self.ccw = cc.Ui_Dialog()
        self.ccw.setupUi(self)
        self.ccw.ccClose.clicked.connect(self.reject)
        self.ccw.createKPPsettings.clicked.connect(self.kpp)
        self.ccw.openOutput.clicked.connect(lambda: qt_box.openOutputDir(None, self.ccw.outDir.text()))
        self.ccw.browseOut.clicked.connect(lambda: qt_box.browse_path(self.ccw.outDir, 'dir'))
        self.ccw.browseSourceFile.clicked.connect(lambda: qt_box.browse_path(self.ccw.sourceFile, 'file'))
        self.ccw.browseReactFile.clicked.connect(lambda: qt_box.browse_path(self.ccw.react_file, 'file'))
        self.ccw.browseIncludes.clicked.connect(lambda: qt_box.browse_path(self.ccw.includedFiles, 'append'))
        self.ccw.mainFrame.setFont(qt_box.font)
        self.ccw.manualCC.clicked.connect(lambda: qt_box.helplink('ccmanual_CC'))
        self.ccw.pickDeffix.clicked.connect(self.ListFixed)
        self.ccw.justReact.toggled.connect(self.ch_txt)
        self.ccw.deffix.setPlainText("SO2\nNO\nNO2\nCO\nH2\nO3")

    def ch_txt(self):
        if self.ccw.justReact.isChecked():
            self.ccw.createKPPsettings.setText("Write reactivity file")
        else:
            self.ccw.createKPPsettings.setText("Create KPP settings")

    def ListFixed(self):
        file = qt_box.browse_path(None, 'fixed_cc', ftype="MCM mass subset file (*.txt)")
        list = qt_box.loadFixedFile(file, cc=True)
        self.ccw.deffix.appendPlainText('\n'.join(list))

    def kpp(self):
        cmds = self.ccw.sourceFile.text()
        if cmds == '':
            qt_box.popup('No input', 'Please give the main source file (e.g. mcm subset).',3)
            return
        if self.ccw.outDir.text() == '' or self.ccw.outDir.text() == './':
            qt_box.popup('No output directory', 'Please provide the output directory.',3)
            return
        includes = self.ccw.includedFiles.toPlainText().split()
        fixed = self.ccw.deffix.toPlainText().split()
        if self.ccw.inclPram.isChecked():
            # includes.append(osjoin(ccloc,'PRAM_v21.txt'))
            includes.append(osjoin('ModelLib','PRAM','PRAM_v21.txt'))
        out = osjoin(self.ccw.outDir.text(),'second.def')
        log = osjoin(self.ccw.outDir.text(),'second.log')
        if cmds[-4:] == '.kpp':
            newold_mcm = 'old'
        elif cmds[-4:] == '.eqn':
            newold_mcm = 'new'
        else:
            qt_box.popup('Cannot determine input syntax version', 'Assuming new type. Files with .kpp are old style mcm output, files with .eqn are new type',1)
            newold_mcm = 'new'
        commandstring = [currentPythonVer,ccloc+'/create_chemistry.py',cmds,'-o',out,'-l',log, '-v',newold_mcm]
        if len(fixed)>0:
            commandstring.append('-d')
            commandstring += fixed
        if len(includes)>0:
            commandstring.append('-f')
            commandstring += includes

        error = False
        warnings = False
        lines = True
        if self.ccw.justReact.isChecked():
            out = self.ccw.sourceFile.text()
        else:
            print( 'Calling chemistry script with:\n'+' '.join(commandstring))
            self.kppProcess = Popen([*commandstring], stdout=PIPE,stderr=STDOUT,stdin=None)
            output = ['Chemistry definitions were created.\n']
            boilerplate = '\n1) Run KPP in the output directory: "kpp second.kpp (or kpp second.kpp3 if kpp v.3 is used)"\n2) Recompile ARCA in tab "Chemistry".'
            while lines:
                self.ccout = self.kppProcess.stdout.readline().decode("utf-8")
                if '[WARNING' in self.ccout.upper():
                    warnings = True
                    output.append('Duplicate equations were found.')
                    output.append(self.ccout)
                if '[CRITICAL' in self.ccout.upper():
                    warnings = True
                    output.append('Included file was not found:')
                    output.append(self.ccout)
                if '[ERROR' in self.ccout.upper():
                    qt_box.popup('Script returned error', 'The script was unable to create chemistry definition, please see the log below.',3)
                    error = True
                self.ccw.ccMonitor.insertPlainText(self.ccout)
                if self.kppProcess.poll() != None and self.ccout == '':
                    lines= False
                    self.kppProcess.kill()
            #
        if error: return
        xml = osjoin(ccloc,'reactivity','reactivity_empty.xml')
        if self.ccw.react_file.text() != '':
            if not exists(self.ccw.react_file.text()):
                qt_box.popup('Reactivity file not found', 'Omitting reactivities.',2)
            else:
                xml = self.ccw.react_file.text()

        commandstring = [currentPythonVer,osjoin(ccloc,'reactivity','add_reactivity.py'),'-k',out,xml,
                        '-o',osjoin(self.ccw.outDir.text(),'second_reactivity.f90'),'--log', osjoin(self.ccw.outDir.text(),'reactivity.log')]
        self.kppProcess = Popen([*commandstring], stdout=PIPE,stderr=STDOUT,stdin=None)
        lines = True
        error = False
        while lines:
            self.ccout = self.kppProcess.stdout.readline().decode("utf-8")
            self.ccw.ccMonitor.insertPlainText(self.ccout)
            # if '[WARNING' in self.ccout.upper():
            #     warnings = True
            #     output.append('Duplicate equations were found.')
            #     output.append(self.ccout)
            # if '[CRITICAL' in self.ccout.upper():
            #     warnings = True
            #     output.append('Included file was not found:')
            #     output.append(self.ccout)
            if 'ERROR' in self.ccout.upper():
                qt_box.popup('Script returned error', 'The script was unable to create reactivity file, please see the log.',3)
                error = True
            if self.kppProcess.poll() != None and self.ccout == '':
                lines= False
                self.kppProcess.kill()

        if error: return
        if self.ccw.justReact.isChecked():
            qt_box.popup('Reactivity file created', 'Success, see the log for details.',0)
        else:
            cpf(osjoin(ccloc,'mcm_module.f90'),osjoin(self.ccw.outDir.text(),'mcm_module.f90'))
            cpf(osjoin(ccloc,'second_Constants.f90'),osjoin(self.ccw.outDir.text(),'second_Constants.f90'))
            cpf(osjoin(ccloc,'second.kpp'),osjoin(self.ccw.outDir.text(),'second.kpp'))
            cpf(osjoin(ccloc,'mcm_module_kpp3.f90'),osjoin(self.ccw.outDir.text(),'mcm_module_kpp3.f90'))
            cpf(osjoin(ccloc,'second.kpp3'),osjoin(self.ccw.outDir.text(),'second.kpp3'))
            if warnings: output.append('Read the log above, and resolve the problems in .def file, then:')
            qt_box.popup('Chemistry created', '\n'.join(output)+boilerplate,0)


# The popup window for variations
class Variation(QtWidgets.QDialog):
    def __init__(self, parent = None):
        super(Variation, self).__init__(parent)
        self.vary = varWin.Ui_Dialog()
        self.vary.setupUi(self)
        self.vary.addLine.clicked.connect(lambda: self.addL())
        self.vary.removeLine.clicked.connect(self.remL)
        self.vary.Close.clicked.connect(self.reject)
        self.vary.runVariations.clicked.connect(self.vars)
        self.vary.Browse.clicked.connect(self.br)
        self.vary.opsFromFile.clicked.connect(self.loadOps)
        self.vary.table.setColumnWidth(0, 90)
        for i in range(1,6):
            self.vary.table.setColumnWidth(i, 70)
        # self.vary.table.horizontalHeader().setStretchLastSection(True)
        self.vary.table.horizontalHeader().setSectionResizeMode(0, QtWidgets.QHeaderView.Stretch)
        self.vary.manualVar.clicked.connect(lambda: qt_box.helplink('variations'))
        self.vary.mainFrame.setFont(qt_box.font)



    def vars(self):
        p=self.vary.lineEdit.text()
        if self.vary.table.rowCount() == 0: return
        ops = zeros((self.vary.table.rowCount(), 7),dtype=object)
        for i in range(self.vary.table.rowCount()):
            for j in range(7):
                if j==0:
                    try:
                        ops[i,j] = float(self.vary.table.item(i,j).text())
                    except:
                        braprt = self.vary.table.item(i,j).text().split(',')
                        if len(braprt)>1:
                            grouped = []
                            for zz in braprt:
                                try:
                                    grouped.append(int(zz))
                                except:
                                    if zz.upper() in self.indices:
                                        grouped.append(self.indices[zz.upper()])
                                    else:
                                        print('Check the input, compound '+zz.upper()+' was not found.')
                                        return
                            print(grouped)
                            ops[i,j] = grouped
                        elif self.vary.table.item(i,j).text().upper() in self.indices:
                            ops[i,j] = self.indices[self.vary.table.item(i,j).text().upper()]
                        else:
                            print('Check the input, compound '+self.vary.table.item(i,j).text().upper()+' was not found.')
                            return
                else:
                    ops[i,j] = float(self.vary.table.item(i,j).text())

        print(variations.zzzz(p, ossplit(p)[0], ops, dryrun=False, nopause=True, replace_current=self.vary.replace_current.isChecked()))

    def loadOps(self):
        path = self.pickF(None)
        if exists(path):
            with open(path) as f:
                for line in f:
                    if not line[0] == '#':
                        if len(line.split()) <= 7:
                            self.addL(line.split())
                        else:
                            print('File has additional junk')

    def pickF(self, ftype):
        dialog = QtWidgets.QFileDialog()
        dialog.setNameFilter(ftype)
        options = dialog.Options()
        options |= dialog.DontUseNativeDialog
        dialog.setOption(QtWidgets.QFileDialog.DontUseNativeDialog)
        dialog.setWindowTitle('Choose File')
        if dialog.exec() == 1:
            path = dialog.selectedFiles()[0]
        else: path=''
        return path

    def br(self):
        path = self.pickF("Arca batch file (*.bash)")
        target = self.vary.lineEdit
        if path != '':
            target.clear()
            target.insert(path)
            self.indices = variations.zzzz(path, ossplit(path)[0],variations.ops, dryrun=True)
            jjj = 0
            print('\nFollowing variables and their indices are picked from chosen bash file\'s')
            print('first run. Use the indices to define the variables that are to be varied:')
            for k in self.indices.keys():
                if jjj%5==0: print('-'*35)
                print('%-16s: %-3d' %(k, self.indices[k]))
                jjj += 1
            print('-'*35)

    def addL(self, cols = ['1','1','1','1','0','0','0']):
        self.vary.table.insertRow(self.vary.table.rowCount())
        for i in range(len(cols)):
            tag = QtWidgets.QTableWidgetItem(cols[i])
            tag.setFont(bold)
            tag.setTextAlignment(QtCore.Qt.AlignCenter)
            self.vary.table.setItem(self.vary.table.rowCount()-1, i, tag)

    def remL(self):
        ii = self.vary.table.selectionModel().selectedIndexes()
        for i in ii:
            self.vary.table.removeRow(i.row())

# The popup window for Vapour pressure file
class VpressWin(QtWidgets.QDialog):
    def __init__(self, parent = None):
        super(VpressWin, self).__init__(parent)
        self.vp = vdialog.Ui_Dialog()
        self.vp.setupUi(self)
        self.vp.VapourClose.clicked.connect(self.reject)
        self.vp.UmanFrame.setEnabled(False)
        self.vp.useUMan.toggled.connect(lambda: qt_box.grayIfNotChecked(self.vp.useUMan,self.vp.UmanFrame))
        self.vp.massSmilesButton.clicked.connect(lambda: qt_box.browse_path(self.vp.lineEdit, 'file',
            ftype="*mcm_subset_mass.txt *mcm_export_species.tsv *.csv *.ssv *.tsv"))
        self.vp.PramButton.clicked.connect(lambda: qt_box.browse_path(self.vp.pramFile, 'file'))
        self.vp.browseVapourPath.clicked.connect(self.filename)
        self.vp.createVapourFileButton.clicked.connect(self.saveVapours)
        self.vp.UmanWWW.clicked.connect(lambda: qt_box.helplink('umanweb'))
        self.vp.manualVap.clicked.connect(lambda: qt_box.helplink('CreateVapourFile'))
        self.vp.mainFrame.setFont(qt_box.font)
        uMan_loc = get_config("paths", "uMan_loc", fallback='The path to local UManSYsProp needs to be specified here')
        self.vp.UManSys_location.setText(uMan_loc)
        self.vp.vpCombo.currentIndexChanged.connect(self.grayEvap)

    def grayEvap(self):
        if self.vp.vpCombo.currentIndex() == 2:
            self.vp.bpCombo.setEnabled(False)
        else:
            self.vp.bpCombo.setEnabled(True)
        return

    def filename(self):
        dialog = QtWidgets.QFileDialog()
        options = dialog.Options()
        options |= dialog.DontUseNativeDialog
        file = dialog.getSaveFileName(self, 'Save Vapours', options=options)[0]
        if file != '': self.vp.VapourPath.setText(file)

    def saveVapours(self):
        """This tool creates the Vapour and Elements files from user input or from the AMG server."""

        set_config("paths", "uMan_loc", self.vp.UManSys_location.text()  )

        filterlist = []
        if self.vp.filterWChem.isChecked():
            chemistry = qt_box.chemistryModules.currentText()
            try:
                count = 0
                with open('src/chemistry/'+chemistry+'/second_Parameters.f90','r') as f:
                    for line in f:
                        x = findall(r'INDF?_\w+',line.upper())
                        if x != []:
                            filterlist.append(x[0].replace('INDF_','').replace('IND_',''))
                            count += 1
            except:
                print('Could not parse current compounds from chemistry... saving all available vapours.')

        if self.vp.VapourPath.text() == '':
            qt_box.popup('Oops...', 'Output filename must be defined.',1)
            return
        if self.vp.useUMan.isChecked(): source = 'UMan'
        else: source = 'AMG'
        if self.vp.limPsat.text() == '' : plim = 1e-6
        else : plim = self.vp.limPsat.text()
        vp_method = {0:'nannoolal',1:'myrdal_and_yalkowsky',2:'evaporation'}
        bp_method = {0:'nannoolal',1:'joback_and_reid',2:'stein_and_brown'}
        if 'mcm_subset_mass.txt' in self.vp.lineEdit.text():
            newold_mcm = 'old'
        elif 'mcm_export_species.tsv' in self.vp.lineEdit.text():
            newold_mcm = 'new'
        elif 'mcm_export_species.tsv' in self.vp.lineEdit.text():
            newold_mcm = 'new'
        elif '.csv' in self.vp.lineEdit.text():
            newold_mcm = 'csv'

        message = gvp.getVaps(args={
        'vp_method':vp_method[self.vp.vpCombo.currentIndex()],
        'bp_method':bp_method[self.vp.bpCombo.currentIndex()],
        'server':source,
        'uMan_loc':self.vp.UManSys_location.text(),
        'smilesfile':self.vp.lineEdit.text(),
        'pram':self.vp.usePRAM.isChecked(),
        'pramfile':self.vp.pramFile.text(),
        'psat_lim':plim,
        'mcm_type':newold_mcm,
        'saveto':self.vp.VapourPath.text(),
        'filter':filterlist,
        'save_atoms':self.vp.saveElements.isChecked()
        })
        if len(message)==1: qt_box.popup('Oops...', 'No dice: '+message[0],3)
        if len(message)==2: qt_box.popup(message[0],message[1],0)

# The popup window for simple input
class Input(QtWidgets.QDialog):
    def __init__(self, parent = None, default = 0):
        super(Input, self).__init__(parent)
        self.inp = input.Ui_Dialog()
        self.inp.setupUi(self)
        self.inp.input.setValue(default)

# The popup window for editing simple text files
class Editor(QtWidgets.QDialog):
    def __init__(self, parent = None, file = None, cursorToStart=False):
        super(Editor, self).__init__(parent)
        self.editor = t_editor.Ui_Dialog()
        self.editor.setupUi(self)
        self.setWindowTitle("Editing "+file)
        f = open(file, 'r')
        t = f.read()
        self.editor.editedText.appendPlainText(t)
        self.editor.editedText.setFont(qt_box.mfont)
        if cursorToStart:
            crs=self.editor.editedText.textCursor()
            crs.movePosition(QtGui.QTextCursor.Start)
            self.editor.editedText.setTextCursor(crs);


# The Class for storing compound input
class Comp:
    """Class for input compounds/variables. Default values are used in Function creator"""
    def __init__(self):
        self.index  = 0
        self.mode   = 0
        self.col    = -1
        self.multi  = 1e0      # Multiplication factor in MODE0
        self.shift  = 0e0      # Constant to be addded in MODE0
        self.min    = 1e1      # Minimum value for the parametrized concentration OR constant value if max <= min
        self.max    = 1e5      # Peak value
        self.sig    = 2.34e0   # Standard deviation for the Gaussian=sig of the bell curve
        self.mju    = 12e0     # Time of peak value
        self.fv     = 0e0      # Angular frequency [hours] of modifying sine function
        self.ph     = 0e0      # Angular frequency [hours] of modifying sine function
        self.am     = 1e0      # Amplitude of modification
        self.tied   = ''       # Variable is tied to this compound
        self.name   = 'NONAME' # Human readable name for modified variable
        self.unit   = '#/cm3'  # unit name
        self.Find   = 1
        self.pmInUse = 'No'
        self.sliderVls = [39,84,0,0,20]
        self.sl_x = [1,1,1,1,1]

class NcPlot:
    """Class for plot file contents"""
    def __init__(self, file, mdim=False):
        self.path = file
        self.masterfile = ossplit(file)[1]
        ncs = netCDF4.Dataset(file, 'r')
        self.legend = getattr(ncs, 'experiment')+': '+self.masterfile
        self.getaircc(file, ncs)
        self.parvars = {}
        self.csat = {}
        self.convars = {}
        self.invvars = {}
        self.par = False
        self.measdmps = False
        if self.masterfile == 'Particles.nc':
            self.par = True
            try:
                # Better way to read netCDF4 strings:
                # netCDF4.chartostring(nc.variables['VAPOURS'][:])
                # Simplify this when time
                b = ravel(ncs.variables['VAPOURS'][:,:].astype(str),'C')
                names=b.reshape(ncs.variables['VAPOURS'].shape)
            except:
                b = ravel(ncs.variables['CONDENSABLES'][:,:].astype(str),'C') # the variable was renamed in v1.2
                names=b.reshape(ncs.variables['CONDENSABLES'].shape)
            self.nc_cm3 = ncs.variables['NUMBER_CONCENTRATION'][:] #
            self.composition_ng = npsum(ncs.variables['PARTICLE_COMPOSITION'][:,:,:]*self.nc_cm3[:,:,newaxis], 1)*1e18 # kg->grams,g->ng cm3->m3=1e18
            self.totalmass_ng_m = npsum(ncs.variables['PARTICLE_COMPOSITION'][:,:,:]*self.nc_cm3[:,:,newaxis], 2)*1e18
            self.diameter       = ncs.variables['DIAMETER'][0,:]
            self.lognorm_nc_cm3 = self.nc_cm3/log10(self.diameter[1]/self.diameter[0]) #
            for i,word in enumerate(names):
                comp = (''.join(list(word))).strip()
                self.parvars[comp] = i
                try:
                        self.csat[comp] = ncs.variables[comp].Psat_A - ncs.variables[comp].Psat_B/300
                        self.csat[comp] = 10**self.csat[comp] * (1e-6 * 101325 / 300 /1.38064852e-23)
                        self.csat_unit = 'cm⁻³'
                except:
                    try:
                        self.csat[comp] = 10** ncs.variables[comp].Csat_300
                        self.csat_unit = 'µg m⁻³'
                    except:
                        self.csat[comp] = - 9999
                        self.csat_unit = '?'
            self.sorter = argsort(list(self.csat.values()))
            try:
                DMPS_CONCENTRATION = ncs.variables['INPUT_CONCENTRATION'][:]
                self.massdmps = ncs.variables['MASS'][:]*DMPS_CONCENTRATION * 1e3 * 1e6 * 1e9 # kg->grams>ng, cm3->m3
                self.lognormdmps = DMPS_CONCENTRATION/log10(self.diameter[1]/self.diameter[0])
                self.measdmps = True
            except:
                print('File does not contain measured PSD')

        for timedim in ncs.dimensions:
            if ncs.dimensions[timedim].isunlimited():
                break

        checker = lambda v,n: v.lower() in timedim and 'Shifter' not in n and 'Multipl' not in n and 'TIME_IN' not in n.upper()
        if mdim:
            cache = array([i.name for i in ncs.get_variables_by_attributes(ndim=lambda d:d in range(1,3))])
            timevars = [checker(i.dimensions[0], i.name) for i in ncs.get_variables_by_attributes(ndim=lambda d:d in range(1,3))]
        else:
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
        try:
            self.timeint = float(getattr(ncs, 'Nominal_save_interval_s'))
        except:
            self.timeint = int(3600*(self.time[-2]-self.time[0])/len(self.time[:-2]))
        # self.conc_matrix = zeros((len(self.time),len(self.varnames[self.mask])))
        for i,n in enumerate(cache[timevars]):
            self.convars[n] = i
            self.invvars[i] = n
            # self.conc_matrix[:,i] = ncs.variables[n][self.mask]

        self.nc = ncs
        self.names = list(self.convars.keys())

    def closenc(self):
        try:
            self.nc.close()
        except:
            print('... Netcdf file was already closed')

    def getinfo(self,n):
        v = self.nc[n]
        sss = 'NetCDF attributes for %s\n\n'%n
        for a in v.ncattrs():
            sss += a+': '+getattr(v,a)+'\n'
        qt_box.popup('Variable info',sss,0)


    def getconc(self,n, return_unit=False):
        if n in self.convars:
            y = self.nc.variables[n][self.mask]
            if len(shape(y))>1:
                y = y[:,0]
            if return_unit:
                return y, '['+units.get(grepunit(n),units['REST'])[0]+']'
            else:
                return y
        else: return

    def getloc(self,i, return_unit=False):
        if i in self.invvars:
            if return_unit:
                return self.nc.variables[self.invvars[i]][self.mask], '['+units.get(grepunit(self.invvars[i]),units['REST'])[0]+']'
            else:
                return self.nc.variables[self.invvars[i]][self.mask]
        else: return

    def getcom(self,n, return_unit=False):
        if self.par and n in self.parvars:
            if return_unit:
                return self.composition_ng[:,self.parvars[n]][self.mask], '['+units.get(grepunit(n),units['REST'])[0]+']'
            else:
                return self.composition_ng[:,self.parvars[n]][self.mask]
        else: return

    def getcomsum(self,names, return_unit=False):
        retarr = zeros(len(self.mask))
        for i,n in enumerate(names):
            if i==0: u = units.get(grepunit(n),units['REST'])[0]
            if self.par and n in self.parvars:
                if u == units.get(grepunit(n),units['REST'])[0]:
                    unit = True
                else:
                    unit = False
                retarr += self.composition_ng[:,self.parvars[n]][self.mask]
        if not unit: u='[-]'
        if return_unit:
            return retarr, u
        else:
            return retarr

    def getconcsum(self,names, return_unit=False):
        retarr = zeros(len(self.mask))
        for i,n in enumerate(names):
            if i==0: u = units.get(grepunit(n),units['REST'])[0]
            if n in self.convars:
                if u == units.get(grepunit(n),units['REST'])[0]:
                    unit = True
                else:
                    unit = False
                retarr += self.getconc(n)
        if not unit: u='[-]'
        if return_unit:
            return retarr, u
        else:
            return retarr

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
            qt_box.popup('Missing air concentration', 'File "General.nc" was not found from this directory. Cannot calculate mixing ratios but will show values.')
            self.have_aircc = False


class QtBoxGui(gui10.Ui_MainWindow,QtWidgets.QMainWindow):
    """Main program window."""
    def __init__(self):
        super(QtBoxGui,self).__init__()
        self.setupUi(self)

    # -----------------------
    # Common stuff
    # -----------------------
        if dark_mode:
            self.saveButton.setStyleSheet("background-image: url(\'ModelLib/gui/icons/dark/saveas.png\'); background-repeat: no-repeat;")
            self.saveCurrentButton.setStyleSheet("background-image: url(\'ModelLib/gui/icons/dark/saveia.png\'); background-repeat: no-repeat;")
            self.loadButton.setStyleSheet("background-image: url(\'ModelLib/gui/icons/dark/load.png\'); background-repeat: no-repeat;")
            self.saveDefaults.setStyleSheet("background-image: url(\'ModelLib/gui/icons/dark/defaults.png\'); background-repeat: no-repeat;")
            self.recompile.setStyleSheet("background-image: url(\'ModelLib/gui/icons/dark/recompile.png\'); background-repeat: no-repeat;")

        self.setAcceptDrops(True)
        self.setWindowTitle(CurrentVersion)
        self.setWindowIcon(QtGui.QIcon(boxicon))
        self.inout_dir.setPlaceholderText("\""+default_inout+"\" if left empty")
        self.case_name.setPlaceholderText("\""+default_case+"\" if left empty")
        self.run_name.setPlaceholderText("\""+default_run+"\" if left empty")
        self.currentInitFileToSave = ''
        self.indir = ''
        self.fileLoadOngoing = False
        self.prints = 1
        self.plots  = 0
        self.wait_for = 0
        self.show_extra_plots = ''
        self.saveButton.clicked.connect(lambda: self.save_file())
        self.saveCurrentButton.clicked.connect(lambda: self.save_file(file=self.currentInitFileToSave, mode='noteOnTer'))
        self.actionSave_to_current.triggered.connect(lambda: self.save_file(file=self.currentInitFileToSave, mode='noteOnTer'))
        self.actionCreate_output_directories.triggered.connect(self.createCaseFolders)
        self.actionLoad_minimal_settings.triggered.connect(lambda: self.load_initfile(minimal_settings_path))
        self.loadButton.clicked.connect(lambda: self.browse_path(None, 'load'))
        self.actionSave_2.triggered.connect(lambda: self.save_file())
        self.actionPrint.triggered.connect(lambda: self.print_values())
        self.actionExport_current_case.triggered.connect(lambda: self.browse_path(None, 'export'))
        self.actionOpen.triggered.connect(lambda: self.browse_path(None, 'load'))
        self.actionQuit_Ctrl_Q.triggered.connect(self.close)
        self.actionSet_monitor_font_2.triggered.connect(lambda: self.guiSetFont(self.MonitorWindow,'monitor'))
        self.actionSet_Global_font.triggered.connect(lambda: self.guiSetFont(self.centralwidget,'global'))
        self.actionReset_fonts.triggered.connect(self.resetFont)
        self.actionCreate_Vapour_file.triggered.connect(self.vapours)
        self.actionCreateNewChemistry.triggered.connect(self.createCC)
        self.actionVariations.triggered.connect(self.createVAR)
        self.actionRecompile_model.triggered.connect(self.remake)
        self.actionSetDelay.triggered.connect(lambda: self.inputPopup("self.wait_for"))
        self.actionAbout_ARCA.triggered.connect(self.createAb)
        self.saveDefaults.clicked.connect(lambda: self.save_file(file=defaults_file_path))
        self.actionSave_as_defaults.triggered.connect(lambda: self.save_file(file=defaults_file_path))
        self.label_10.setPixmap(QtGui.QPixmap(modellogo))
        self.actionPrint_input_headers.triggered.connect(self.printHeaders)
        self.actionOpen_output_directory.triggered.connect(lambda: self.openOutputDir(None, self.currentAddressTb.text()))
        self.actionPrint_Custom_commands_cheat_sheet.triggered.connect(lambda: CustomCommandsCheatSheet())
        self.actionPlt_changes_from_current_dir.triggered.connect(self.plotChanges)
        self.actionShow_variable_attributes.triggered.connect(lambda: self.showOutputUpdate(info=True))
        if dark_mode: self.actionAcommodateForDarkMode.setChecked(True)
        self.actionAcommodateForDarkMode.triggered.connect(self.darkMode)
        self.actionRunARCA.triggered.connect(self.StartboxShortcut)
        self.actionStopCurrentRunAndIgnoreOutput.triggered.connect(self.QuitShortcut)
        self.actionStopCurrentRunClean.triggered.connect(self.QuitGracefullyShortcut)
    # -----------------------
    # tab General options
    # -----------------------
        self.namesdat.clear()
        self.namesdat.addItems(NAMES)
        for i in range(len(NAMES)):
            item = self.namesdat.item(i)
            if i<divider_i:
                # item.setForeground(QtGui.QColor(0, 70, 0))
                # item.setBackground(QtGui.QColor(200, 230, 200))
                item.setBackground(QtGui.QColor(*env_no))
            elif i>divider_xtr_i:
                # item.setForeground(QtGui.QColor(0, 70, 0))
                # item.setBackground(QtGui.QColor(200, 230, 200))
                item.setBackground(QtGui.QColor(*xtr_no))
            elif i==divider_i or i==divider_xtr_i:
                item.setFlags(item.flags() & ~QtCore.Qt.ItemIsEnabled & ~QtCore.Qt.ItemIsSelectable)
                item.setBackground(QtGui.QColor(*org_yes))
                item.setForeground(QtGui.QColor(250, 250, 250))
            else:
                # item.setBackground(QtGui.QColor(200, 200, 230))
                item.setBackground(QtGui.QColor(*org_no))
                # item.setForeground(QtGui.QColor(0, 0, 70))

        self.runtime.valueChanged.connect(lambda: self.updteGraph())
        self.runtime.valueChanged.connect(lambda: self.runtime_s.setValue(int(self.runtime.value()*3600)))
        self.runtime_s.editingFinished.connect(lambda: self.runtime.setValue(self.runtime_s.value()/3600))

        # Prepare the variable table
        for i in range(len(column_widths)):
            self.selected_vars.setColumnWidth(i, column_widths[i])

        # add minimum requirements
        self.add_new_line('TEMPK', 0)
        self.add_new_line('PRESSURE', 1)
        self.dateEdit.dateChanged.connect(lambda: self.indexRadioDate.setChecked(True))

        self.indexRadioDate.toggled.connect(lambda: self.lat.setEnabled(True))
        self.indexRadioIndex.toggled.connect(lambda: self.lat.setEnabled(False))
        self.indexRadioDate.toggled.connect(lambda: self.lon.setEnabled(True))
        self.indexRadioIndex.toggled.connect(lambda: self.lon.setEnabled(False))
        self.indexRadioDate.toggled.connect(lambda: self.label_9.setEnabled(True))
        self.indexRadioIndex.toggled.connect(lambda: self.label_9.setEnabled(False))
        self.indexRadioDate.toggled.connect(lambda: self.label_25.setEnabled(True))
        self.indexRadioIndex.toggled.connect(lambda: self.label_25.setEnabled(False))

        self.indexEdit.valueChanged.connect(lambda: self.indexRadioIndex.setChecked(True))
        self.browseCommonIn.clicked.connect(lambda: self.browse_path(self.inout_dir, 'dir'))
        self.browseEnv.clicked.connect(lambda: self.browse_path(self.env_file, 'file'))
        self.browseMcm.clicked.connect(lambda: self.browse_path(self.mcm_file, 'file'))
        self.browsePar.clicked.connect(lambda: self.browse_path(self.dmps_file, 'file'))
        self.browseXtr.clicked.connect(lambda: self.browse_path(self.extra_particles, 'file'))
        self.checkBox_aer.toggled.connect(lambda: self.grayIfNotChecked(self.checkBox_aer,None, 4))
        self.checkBox_che.stateChanged.connect(lambda: self.grayIfNotChecked(self.checkBox_che,None, 2))
        self.fsave_division.valueChanged.connect(self.toggle_printtime)
        self.checkBox_acd.stateChanged.connect(lambda: self.grayIfNotChecked(self.checkBox_acd,None, hide=3))
        self.dateEdit.dateChanged.connect(self.updatePath)
        self.indexEdit.valueChanged.connect(self.updatePath)
        self.case_name.textChanged.connect(self.updatePath)
        self.case_name.editingFinished.connect(lambda: self.allcaps(self.case_name))
        self.run_name.textChanged.connect(self.updatePath)
        self.run_name.editingFinished.connect(lambda: self.allcaps(self.run_name))
        self.inout_dir.textChanged.connect(self.updatePath)
        self.indexRadioDate.toggled.connect(self.updatePath)
        self.useSpeed.toggled.connect(self.adjust_dt)
        self.dateEdit.dateChanged.connect(self.updateEnvPath)
        self.dateEdit.dateChanged.connect(lambda: self.curDate.setText(self.dateEdit.text()))
        self.indexEdit.valueChanged.connect(self.updateEnvPath)
        self.indexRadioDate.toggled.connect(self.updateEnvPath)
        self.env_file.textChanged.connect(self.updateEnvPath)
        self.mcm_file.textChanged.connect(self.updateEnvPath)
        self.dmps_file.textChanged.connect(self.updateEnvPath)
        self.extra_particles.textChanged.connect(self.updateEnvPath)
        self.stripRoot_env.toggled.connect(self.updateEnvPath)
        self.stripRoot_mcm.toggled.connect(self.updateEnvPath)
        self.stripRoot_par.toggled.connect(self.updateEnvPath)
        self.stripRoot_xtr.toggled.connect(self.updateEnvPath)
        self.saveCurrentButton.setEnabled(False)
        self.actionSave_to_current.setEnabled(False)
        self.currentInitFile.setText('None loaded/saved')

        self.min_particle_diam.textChanged.connect(self.seeInAction)
        self.max_particle_diam.textChanged.connect(self.seeInAction)
        self.n_bins_particle.valueChanged.connect(self.seeInAction)
        self.fileTimeUnit_a.currentIndexChanged.connect(lambda: self.fileTimeUnit_b.setCurrentIndex(self.fileTimeUnit_a.currentIndex()))
        self.fileTimeUnit_b.currentIndexChanged.connect(lambda: self.fileTimeUnit_a.setCurrentIndex(self.fileTimeUnit_b.currentIndex()))
    # -----------------------
    # tab Input variables
    # -----------------------
        self.butMoveToSelVars.clicked.connect(self.select_compounds)
        self.markAll.clicked.connect(lambda: self.markReverseSelection('all'))
        self.invertMarks.clicked.connect(lambda: self.markReverseSelection('inv'))
        self.butRemoveSelVars.clicked.connect(self.remv_item)
        self.selected_vars.setColumnHidden(7, True)
        self.selected_vars.verticalHeader().setVisible(False);
        self.loadFixed.clicked.connect(lambda: self.browse_path(None, 'fixed', ftype="MCM mass subset file (*.txt)"))
        self.loadFixedChemistry.clicked.connect(self.loadFixedFromChemistry)
        self.findInput.textChanged.connect(self.filterListOfInput)

    # -----------------------
    # tab Function creator
    # -----------------------

        self.PLOT.setMenuEnabled(False)
        self.PLOT.setBackground('w')
        self.PLOT.setLabel('bottom','time',units='h')
        self.PLOT.getAxis('left').setStyle(autoExpandTextSpace=False)
        self.PLOT.getAxis('left').setWidth(w=60)
        self.confirm.hide()
        self.plotTo.clicked.connect(self.select_compounds_for_plot)
        self.saveParams.clicked.connect(self.saveParamValues)
        self.loadParams.clicked.connect(self.loadParamValues)
        self.fMin.editingFinished.connect(lambda: self.updteGraph())
        self.fMax.editingFinished.connect(lambda: self.updteGraph())
        self.fLog.clicked.connect(lambda: self.updteGraph())
        self.fLin.clicked.connect(lambda: self.updteGraph())
        self.fWidth.valueChanged.connect(lambda: self.updteGraph())
        self.scaleFs = [self.wScalei,self.peScale,self.anScale,self.phScale,self.amScale]
        self.sliders = [self.fWidth,self.fPeak,self.fFreq,self.fPhase,self.fAmp]
        for i in range(5):
            self.resetSlider(self.sliders[i], defCompound.sliderVls[i])

        self.wScalei.valueChanged.connect(lambda: self.fWidth.setMaximum(int(max(self.wScalei.value(), 0.5)*slMxs[0])))
        self.peScale.valueChanged.connect(lambda: self.fPeak.setMaximum(int(max(self.peScale.value(), 0.5)*slMxs[1])))
        self.anScale.valueChanged.connect(lambda: self.fFreq.setMaximum(int(max(self.anScale.value(), 0.5)*slMxs[2])))
        self.phScale.valueChanged.connect(lambda: self.fPhase.setMaximum(int(max(self.phScale.value(), 0.5)*slMxs[3])))
        self.amScale.valueChanged.connect(lambda: self.fAmp.setMaximum(int(max(self.amScale.value(), 0.5)*slMxs[4])))
        self.resW.clicked.connect(lambda: self.resetSlider(self.fWidth, defCompound.sliderVls[0]))
        self.resP.clicked.connect(lambda: self.resetSlider(self.fPeak,  defCompound.sliderVls[1]))
        self.resA.clicked.connect(lambda: self.resetSlider(self.fFreq,  defCompound.sliderVls[2]))
        self.resPh.clicked.connect(lambda: self.resetSlider(self.fPhase,defCompound.sliderVls[3]))
        self.resAm.clicked.connect(lambda: self.resetSlider(self.fAmp,  defCompound.sliderVls[4]))
        self.fPeak.valueChanged.connect(lambda: self.updteGraph())
        self.fFreq.valueChanged.connect(lambda: self.updteGraph())
        self.fPhase.valueChanged.connect(lambda: self.updteGraph())
        self.fAmp.valueChanged.connect(lambda: self.updteGraph())
        self.PLOT.showGrid(x=True,y=True)
        self.PLOT.showButtons()
        self.legend = self.PLOT.addLegend()
        self.skene = self.legend.scene()
        self.second = False
        self.updteGraph(first=True)

    # -----------------------
    # tab Chemistry
    # -----------------------
        self.ReplChem.setChecked(False)
        self.butSpectralData.clicked.connect(lambda: self.browse_path(self.spectralFunctions, 'file'))
        self.kppTool.clicked.connect(self.createCC)
        self.vapTool.clicked.connect(self.vapours)
        self.recompile.clicked.connect(self.remake)
        self.TimerCompile = QtCore.QTimer(self);
        self.TimerCompile.timeout.connect(self.progress)
        self.compileProgressBar.hide()
        self.running = 0
        self.get_available_chemistry(checkonly=False)
        self.chemLabel.setText('Current chemistry scheme in makefile: '+self.get_available_chemistry(checkonly=True))
        self.spectralFunctions.textChanged.connect(lambda: \
                                        self.frame_27.setEnabled(False) if (self.spectralFunctions.text()!='' \
                                        and ossplit(self.spectralFunctions.text())[1]!=defaultSpectrum) \
                                        else self.frame_27.setEnabled(True))

        self.SW_is_AF.toggled.connect(lambda: self.popup('Warning','If irradiation is interpreted as Actinic flux, default spectrum is invalid. \
Please provide valid spectral function.') \
                    if ((self.spectralFunctions.text()=='' or ossplit(self.spectralFunctions.text())[1]==defaultSpectrum) \
                    and self.SW_is_AF.isChecked()) else False)

        self.SW_is_AF.toggled.connect(lambda: self.grayIfChecked(self.SW_is_AF,self.groupBox_23))

    # -----------------------
    # tab Aerosols
    # -----------------------
        # self.frameBase.setEnabled(False)
        self.butVapourNames.clicked.connect(lambda: self.browse_path(self.vap_names, 'file'))
        self.butVapourAtoms.clicked.connect(lambda: self.browse_path(self.vap_atoms, 'file'))
        self.use_atoms.stateChanged.connect(lambda: self.grayIfNotChecked(self.use_atoms,self.vap_atoms))
        self.use_atoms.stateChanged.connect(lambda: self.grayIfNotChecked(self.use_atoms,self.butVapourAtoms))
        # self.use_dmps_partial.toggled.connect(lambda: self.toggle_gray(self.use_dmps_partial,self.gridLayout_11))
        self.saveBatch.clicked.connect(lambda: self.batchCaller())
        self.batchFrFile.clicked.connect(lambda: self.browse_path(self.ListbatchCaller, 'batchList'))
        self.batchRangeDayBegin.dateChanged.connect(lambda: self.batchRangeDay.setChecked(True))
        self.batchRangeDayEnd.dateChanged.connect(lambda: self.batchRangeDay.setChecked(True))
        self.batchRangeIndBegin.valueChanged.connect(lambda: self.batchRangeInd.setChecked(True))
        self.batchRangeIndEnd.valueChanged.connect(lambda: self.batchRangeInd.setChecked(True))
        # self.chemistryModules.setEnabled(False)
        self.ReplChem.clicked.connect(lambda: self.get_available_chemistry(checkonly=False))
        self.mmodal_input.textChanged.connect(self.seeInAction)
        self.n_modal.textChanged.connect(self.seeInAction)
        self.multiModalBox.toggled.connect(lambda: self.HPLotter.setBackground(mmc[str(self.multiModalBox.isChecked())]))
        self.mmodal_input.textChanged.connect(lambda: self.HPLotter.setBackground(mmc[str(self.multiModalBox.isChecked())]))
        self.OrgNuclVapFile.clicked.connect(lambda: self.editTxtFile(nucl_homs))
        self.viewCustomNucl.clicked.connect(lambda: self.editTxtFile(custom_functions))
        self.viewCustomChem.clicked.connect(lambda: self.editTxtFile(custom_functions))
        self.use_dmps.toggled.connect(lambda: self.grayIfChecked(self.use_dmps,self.multiModalBox))
        self.multiModalBox.toggled.connect(lambda: self.grayIfChecked(self.multiModalBox,self.use_dmps))
        self.EditAcdcGen.clicked.connect(lambda: self.editTxtFile('src/ACDC/'+self.CurrentAcdcSystem+"/settings.conf"))
        self.openPerlDir.clicked.connect(lambda: self.openOutputDir(None, 'src/ACDC/'+self.CurrentAcdcSystem+'/Perl_input'))
        self.runPerlScript.clicked.connect(lambda: self.rerunACDC(make=False))
        self.runPerlScript_and_make.clicked.connect(lambda: self.rerunACDC(make=True))
        self.ACDC_available_compounds = []
        self.acdc_systems_flags = [1,1,0,0,0]
        self.ACDC_current_links = [ [0, {'A': 'H2SO4', 'N': 'NH3'}  ],
                                    [0, {'A': 'H2SO4', 'D': 'DMA'}  ],
                                    [0, {'A': '-', 'N': '-'}        ],
                                    [0, {'A': '-', 'N': '-'}        ],
                                    [0, {'A': '-', 'N': '-'}        ] ]
        self.ACDC_n_systems = 5
        self.nickname = []
        self.parse_ACDC_systems()
        self.editACDC.setEnabled(False)
        self.CurrentAcdcSystem = ''
        self.ACDC_linker.setCurrentItem(self.ACDC_linker.topLevelItem(0))
        self.ACDC_linker.itemClicked.connect(self.setCurrentAcdcSystem)
        for i in range(2):
            self.ACDC_linker.setColumnWidth(i,[150,50][i]);

    # -----------------------
    # tab Losses
    # -----------------------
        self.floorArea.valueChanged.connect(self.update_title)
        self.chamberHeight.valueChanged.connect(self.update_title)
        self.browseLosses.clicked.connect(lambda: self.browse_path(self.losses_file, 'file'))
        self.losses_file.textChanged.connect(lambda: \
                                        self.grayIfChecked(True,self.deposition) if self.losses_file.text()!='' \
                                        else self.grayIfChecked(False,self.deposition))

    # -----------------------
    # tab Process Monitor
    # -----------------------
        self.currentEndLine = 0
        self.fulltext = ''
        fixedFont = QtGui.QFontDatabase.systemFont(QtGui.QFontDatabase.FixedFont)
        fixedFont.setPointSize(10)
        self.MonitorWindow.setFont(fixedFont)
        self.frameStop.setEnabled(False)
        self.startButton.clicked.connect(self.startBox)
        self.stopButton.clicked.connect(self.stopBox)
        self.softStopButton.clicked.connect(self.softStop)
        self.boxProcess = 0 # arcabox run handle
        self.monStatus  = 0
        self.Timer = QtCore.QTimer(self);
        self.pollTimer = QtCore.QTimer(self);
        self.Timer.timeout.connect(self.updateOutput)
        self.pollTimer.timeout.connect(self.pollMonitor)
        # self.pollTimer.timeout.connect(self.updateOutput)
        self.pauseScroll.clicked.connect(lambda: self.MonitorWindow.verticalScrollBar().setSliderPosition(self.MonitorWindow.verticalScrollBar().maximum()))
        self.mfont = self.MonitorWindow.font()
        self.mfont.setFamily(    get_config("fonts", "monitorfont",  fallback=str(self.mfont.family())   ) )
        self.mfont.setPointSize( int (get_config("fonts", "monitorsize",  fallback=self.mfont.pointSize()) ) )
        self.mfont.setBold(   eval("%s"%(get_config("fonts", "monitorbold",  fallback=str(self.mfont.bold()) ).capitalize()  ) ))
        self.mfont.setItalic( eval("%s"%(get_config("fonts", "monitoritalic",fallback=str(self.mfont.italic()) ).capitalize()  ) ))
        self.MonitorWindow.setFont(self.mfont)
        self.savefonts(self.mfont, 'monitor')
        self.font = self.centralwidget.font()
        self.font.setFamily(    get_config("fonts", "globalfont",  fallback=str(self.font.family())   ) )
        self.font.setPointSize( int (get_config("fonts", "globalsize",  fallback=self.font.pointSize()) ) )
        self.font.setBold(   eval("%s"%(get_config("fonts", "globalbold",  fallback=str(self.font.bold()) ).capitalize()  ) ))
        self.font.setItalic( eval("%s"%(get_config("fonts", "globalitalic",fallback=str(self.font.italic()) ).capitalize()  ) ))
        self.centralwidget.setFont(self.font)
        self.savefonts(self.font, 'global')
        self.viewPrintNML.clicked.connect(lambda: self.editTxtFile(SCREENPRINT_NML))
    # -----------------------
    # tab Output Graph
    # -----------------------
        if netcdf:
            self.show_netcdf.hide()
            self.show_netcdf_2.hide()
            self.fLog_2.clicked.connect(lambda: self.showOutputUpdate(info=False))
            self.fLin_2.clicked.connect(lambda: self.showOutputUpdate(info=False))
            self.findComp.textChanged.connect(self.filterListOfComp)
            self.loadNetcdf.clicked.connect(lambda: self.browse_path(None, 'addplot', ftype="NetCDF (*.nc)"))
            self.addSimilar.clicked.connect(self.addAnotherNC)
            self.loadNetcdf_mass.clicked.connect(lambda: self.browse_path(None, 'plot_mass', ftype="ARCA particle file (Particles.nc)"))
            self.loadNetcdf_massAdd.clicked.connect(lambda: self.browse_path(None, 'plot_mass_add', ftype="ARCA particle file (Particles.nc)"))
            self.loadNetcdf_massAdd.setEnabled(False)
            self.loadNetcdfPar.clicked.connect(lambda: self.browse_path(None, 'plotPar', ftype="NetCDF, sum (*.nc *.sum *.dat)",plWind=0))
            self.loadSumPar.clicked.connect(lambda: self.browse_path(None, 'plotPar', ftype="NetCDF, sum (*.nc *.sum *.dat)",plWind=1))
            # --------------------------------------------------------
            self.CloseLinePlotsButton.clicked.connect(self.closenetcdf)
            # self.CloseLinePlotsButton.clicked.connect(lambda: self.browse_path(None, 'addplot', ftype="NetCDF (*.nc)"))
            self.listOfplottedFiles = []
            self.LPD = []
            self.NC_lines = []
            self.ncleg = self.plotResultWindow.addLegend()
            self.ncleg_skene = self.ncleg.scene()
            self.massleg = self.plotResultWindow_2.addLegend()
            self.massleg_skene = self.massleg.scene()

        else:
            self.sumSelection.setEnabled(False)
            self.show_netcdf.show()
            self.fLin_2.setEnabled(False)
            self.fLog_2.setEnabled(False)
            self.loadNetcdf.clicked.connect(lambda: self.popup(*netcdfMissinnMes))
            self.loadNetcdf_mass.clicked.connect(lambda: self.popup(*netcdfMissinnMes))
            self.loadNetcdfPar.clicked.connect(lambda: self.browse_path(None, 'plotPar', ftype="sum (*.sum *.dat)",plWind=0))
            self.loadSumPar.clicked.connect(lambda: self.browse_path(None, 'plotPar', ftype="sum (*.sum *.dat)",plWind=1))

        pen = pg.mkPen(color=(0,0,0), width=1)
        for wnd in [self.plotResultWindow, self.plotResultWindow_2, self.plotResultWindow_3]:
            wnd.showGrid(x=True,y=True)
            wnd.setBackground('w')
            wnd.getAxis('left').setPen(pen)
            wnd.getAxis('bottom').setPen(pen)

        self.addSimilar.setEnabled(False)
        self.actionShow_variable_attributes.setEnabled(False)
        self.CloseLinePlotsButton.setEnabled(False)
        self.findComp.setEnabled(False)
        self.ppm.toggled.connect(lambda: self.showOutputUpdate(info=False))
        self.ppb.toggled.connect(lambda: self.showOutputUpdate(info=False))
        self.ppt.toggled.connect(lambda: self.showOutputUpdate(info=False))
        self.toggleppm('off')
        self.ShowPPC.setEnabled(False)
        self.sumSelection.stateChanged.connect(self.selectionMode)
        self.loadCurrentBg.clicked.connect(lambda: self.showParOutput('load current',1))
        self.oneDayFwd.clicked.connect(lambda: self.moveOneDay(1))
        self.oneDayBack.clicked.connect(lambda: self.moveOneDay(-1))
        self.ncs_mass = 0
        self.showAlsoMeasInMassConc.stateChanged.connect(self.updateMass)
        self.showAlsoMeasInMassConc.stateChanged.connect(self.updateNumbers)
        self.pollTimer.timeout.connect(self.livePlot)
        self.liveUpdate.setEnabled(False)
        self.lastModTime = 0
        self.firstParPlot = [0,0]
        self.resize(980, 840)

    # -----------------------
    # Help links
    # -----------------------
        self.helpGroupName.clicked.connect(lambda: self.helplink('groupName'))
        self.helpRunName.clicked.connect(lambda: self.helplink(  'runName'))
        self.helpBatch.clicked.connect(lambda: self.helplink(    'batch'))
        self.helpPrec.clicked.connect(lambda: self.helplink(     'precision'))
        self.helpInput.clicked.connect(lambda: self.helplink(    'input'))
        self.helpVariables.clicked.connect(lambda: self.helplink(     'variables'))
        self.helpParametric.clicked.connect(lambda: self.helplink(    'parametric'))
        self.helpSwitchChem.clicked.connect(lambda: self.helplink(    'switchChem'))
        self.helpCc.clicked.connect(lambda: self.helplink(            'ccmanual'))
        self.helpCustomf.clicked.connect(lambda: self.helplink(       'customFunc'))
        # self.helpACDC.clicked.connect(lambda: self.helplink(          'acdc1'))
        self.helpACDCAdv.clicked.connect(lambda: self.helplink(       'acdc2'))
        self.helpOrgNucl.clicked.connect(lambda: self.helplink(       'orgnucl'))
        self.helpCustomf2.clicked.connect(lambda: self.helplink(      'customFunc'))
        self.helpVapour.clicked.connect(lambda: self.helplink(        'CreateVapourFile'))
        self.helpPsdInit.clicked.connect(lambda: self.helplink(       'psdinit'))
        self.helpMultiModal.clicked.connect(lambda: self.helplink(    'mmodal'))
        self.helpDmpsSpecial.clicked.connect(lambda: self.helplink(   'dmpsSpec'))
        self.helpPsd.clicked.connect(lambda: self.helplink(           'psdScheme'))
        self.helpLosses.clicked.connect(lambda: self.helplink(        'losses'))
        self.helpCustom.clicked.connect(lambda: self.helplink(        'custominput'))
        self.actionARCA_webpage.triggered.connect(lambda: self.helplink('arcaweb'))
        self.actionOnline_manual.triggered.connect(lambda: self.helplink('manual'))
        self.actionFileHelp.triggered.connect(lambda: self.helplink('filehelp'))
        self.helpRadiation.clicked.connect(lambda: self.helplink('radiation'))
    # -----------------------
    # Load preferences, or create preferences if not found
    # -----------------------
        # try:
        #     # if defaults exist, use them
        if startupInitfile is not None and exists(startupInitfile):
            out = self.load_initfile(startupInitfile)
            self.show_currentInit(startupInitfile)
        else:
        if exists(defaults_file_path):
            out = self.load_initfile(defaults_file_path)
        else:
            try:
                # If defaults did not exist, load minimal working settings
                out = self.load_initfile(minimal_settings_path)
            except:
                # If they also were missing, use "factory settings"
                pass
            # Save the obtained settings as default
            self.save_file(file=defaults_file_path, mode='silent')
        if out==None: self.updateEnvPath()

    # -----------------------
    # Class methods
    # -----------------------

    def dragEnterEvent(self, e):
        if e.mimeData().hasText(): e.accept()
        else: e.ignore()
    def dragMoveEvent(self, e):
        if e.mimeData().hasText(): e.accept()
        else: e.ignore()
    def dropEvent(self, e):
        if operatingsystem == 'Windows': preF = 'file:///'
        else: preF = 'file://'
        file = e.mimeData().text().replace(preF,'').strip()
        out = self.load_initfile(file)
        if out==None: self.show_currentInit(file)


    def rerunACDC(self, make=False):
        bashscript = "src/ACDC/run_perl.sh"
        cdto = "src/ACDC/"+self.CurrentAcdcSystem
        self.rra = Popen(["bash", bashscript, "settings.conf",cdto])
        while True:
            xx = self.rra.poll()
            if xx != None: break
        num = int(self.CurrentAcdcSystem[-1])
        self.parse_ACDC_systems(num)
        if make:
            self.remake(syst=self.CurrentAcdcSystem[-1])


    def setCurrentAcdcSystem(self):
        title = 'Actions for the highlighted system'
        b1 = 'View/edit settings'
        b2 = 'Open system input directory'
        b3 = 'Run ACDC Perl script'
        b4 = 'Run ACDC Perl and recompile ARCA'
        if self.ACDC_linker.currentItem().parent() == None and self.ACDC_linker.currentIndex().data() != None:
            self.CurrentAcdcSystem = self.ACDC_linker.currentIndex().siblingAtColumn(0).data()
            self.editACDC.setEnabled(True)
            self.editACDC.setTitle(title+' ('+self.CurrentAcdcSystem+')')
            self.EditAcdcGen.setText(b1 + ' for ' + self.CurrentAcdcSystem)
            self.openPerlDir.setText(b2 + ' of ' + self.CurrentAcdcSystem)
            self.runPerlScript.setText('Update  *f.90 for ' + ' ' + self.CurrentAcdcSystem)
            self.runPerlScript_and_make.setText('Update  *f.90 for ' + ' ' + self.CurrentAcdcSystem + ' and recompile ARCA')
        else:
            self.CurrentAcdcSystem = ''
            self.editACDC.setEnabled(False)
            self.editACDC.setTitle(title)
            self.EditAcdcGen.setText(b1)
            self.openPerlDir.setText(b2)
            self.runPerlScript.setText(b3)
            self.runPerlScript_and_make.setText(b4)

    def openOutputDir(self, a, dir):
        if dir == '': dir = './'
        if not exists(dir):
            self.popup('Sorry','Directory does not exist.',2)
            return
        if operatingsystem == 'Windows':
            os.startfile(dir)
        if operatingsystem == 'Linux':
            os.system('xdg-open "%s"' % dir)
        if operatingsystem == 'Darwin':
            os.system('open "%s"' % dir)
        else:
            return


    def helplink(self, linkStr):
        if linkStr in helpd:
            QtGui.QDesktopServices.openUrl(QtCore.QUrl(helpd[linkStr]))
        else:
            self.popup('This is embarassing', 'Help link is missing. Try Toolbar->Help->Online manual.')


    def adjust_dt(self):
        if self.useSpeed.isChecked():
            self.dt.setValue(0.001)
        else:
            self.dt.setValue(10.0)

    def savefonts(self, font, name):
        set_config("fonts", "%sfont"%name,   str(font.family())   )
        set_config("fonts", "%ssize"%name,   str(font.pointSize()))
        set_config("fonts", "%sbold"%name,   str(font.bold())     )
        set_config("fonts", "%sitalic"%name, str(font.italic())   )

    def guiSetFont(self, wdgt, name, reset=False):
        if not reset:
            dialog = QtWidgets.QFontDialog()
            font, ok = dialog.getFont(wdgt.font(), parent=self)
            if ok:
                wdgt.setFont(font)
                self.savefonts(font, name)
        else:
            font = QtGui.QFont()
            font.setBold(False)
            font.setItalic(False)
            font.setWeight(50)
            if name == 'monitor':
                if osname.upper() == 'NT':
                    font.setFamily('Consolas')
                else:
                    font.setFamily('monospace')
                font.setStyleStrategy(QtGui.QFont.PreferDefault)
                self.mfont = font
            else:
                self.font = font
            wdgt.setFont(font)
            self.savefonts(font, name)

    def splot(self,vals,N, nb, x0,x1):
        """Harry Plotter"""
        def gaussian(x, mu, sig, A=1):
            return A*exp(-(x-mu)**2/(2*sig**2))/sqrt(2*pi*sig**2)
        try:
            N = float(N)
            luvut = array(vals.split()).astype(float)
            nb=int(nb)
            x0=float(x0)
            x1=float(x1)
            if x0 == 0 or x1 == 0:
                return
            if len(luvut) < 3:
                return
        except:
            return
        # X = diameters of PSD
        x = 10**linspace(log10(x0),log10(x1),nb)
        acl = zeros(len(x))
        k = log10(x[1]/x[0])
        for i in range(len(luvut)//3):
            if abs(luvut[3*i+1])>0:
                acl = acl + luvut[3*i+2]*(gaussian(log10(x), log10(luvut[3*i+0]),luvut[3*i+1]))
                # will change to this form asap
                # acl = acl + luvut[3*i+2]*(gaussian(log(x), log(luvut[3*i+0]),luvut[3*i+1]))
        if sum( acl )>0:
            Z = N * acl / (sum( acl ))
            ndel = -min(nb-1, 8)
            Z[ndel:] = where(Z[ndel:]>1e-12, 0,Z[ndel:])
            self.massLabel.setText('Total mass: %9.2e µg/m3 (Density: 1.4g/cm³)'%(npsum(Z*(pi*x**3)/6)*1.4e3*1e15))
            self.areaLabel.setText('Total Area: %9.2e µm2/cm3 '%(npsum(Z*(pi*x**2))*1e12))
            self.HPLotter.plot(x,Z/k,pen=pg.mkPen('w', width=4), clear=True, name='PSD')
            self.HPLotter.setLogMode(x=True)
            self.HPLotter.showGrid(x=True,y=True)


    def darkMode(self):
        msg = "The GUI is using system theme, but some colors are \
adjusted for MODE mode. The settings are applied after the next application launch."
        if self.actionAcommodateForDarkMode.isChecked():
            self.popup('Dark mode', msg.replace('MODE', 'dark'), 0)
            set_config("style", "darkmode",  'True')
        else:
            self.popup('Light mode', msg.replace('MODE', 'light'), 0)
            set_config("style", "darkmode",  'False')


    def resetFont(self):
        self.guiSetFont(self.MonitorWindow, 'monitor', reset=True)
        print(self.MonitorWindow.font().family())
        self.guiSetFont(self.centralwidget, 'global', reset=True)

    def allcaps(self, tbox):
        tbox.setText(tbox.text().upper())

    def exportCurrentCase(self, InitFileFull):
        path, Initfile = ossplit(InitFileFull)
        if path == '':
            self.popup('Sorry but no', 'Don\'t export to program root', 3)
            return
        if not exists(path+'/exportedData'): mkdir(path+'/exportedData')
        self.updatePath()
        if self.checkBox_aer.isChecked() and not exists(self.vap_names.text()):
            self.popup('Error', 'Vapour file is not found but aerosol module is ON.', icon=3)
            return
        if (self.use_atoms.isChecked() and self.checkBox_aer.isChecked()) and not exists(self.vap_atoms.text()):
            self.popup('Error', 'Vapour elemental composition file is not found but elemental composition is used.', icon=3)
            return
        curInit = self.currentInitFile.text()
        inout = self.inout_dir.text()
        indir = self.indir
        justInOut = ossplit(path)[1]
        self.inout_dir.setText('INOUT/'+justInOut)
        paths = [
            0,                              0,self.vap_names,
            0,                              0,self.vap_atoms,
            self.stripRoot_env.isChecked(), 1,self.env_file,
            self.stripRoot_mcm.isChecked(), 1,self.mcm_file,
            self.stripRoot_par.isChecked(), 1,self.dmps_file,
            False,                          1,self.losses_file
        ]
        original_names = []
        for i,f in enumerate(paths):
            if i%3==2:
                original_names.append(f.text())
                if paths[i-1]==0:
                    real = f.text()
                    if real != '':
                        if exists(real):
                            cpf(real,path+'/exportedData/'+ossplit(real)[1])
                        j = f.text().rstrip('/')
                        f.setText('INOUT/'+justInOut+'/exportedData/'+ossplit(real)[1])
                else:
                    real = self.pars(f.text(), file=indir, stripRoot=paths[i-2])
                    if real != '':
                        if exists(real):
                            cpf(real,path+'/exportedData/'+ossplit(real)[1])
                        t = self.pars(f.text(), file=indir, stripRoot=True)
                        f.setText('INOUT/'+justInOut+'/exportedData/'+ossplit(t)[1])

        self.updatePath()
        self.save_file(file=InitFileFull, mode='silent', changeTitle=False)
        for o,t in zip(original_names, paths[2::3]):
            t.setText(o)
        self.inout_dir.setText(inout)
        self.currentInitFile.setText(curInit)
        self.updatePath()
        self.popup('Success!', '''The current case was succesfully exported. When the exported settings
are loaded, make sure to update the "Common out" directory to correspond the new
location. In theory no other modifications should be necessary. Also note that if
the numerical model or chemistry scheme differs from the current, results may vary''', icon=0)

        return


    def livePlot(self):
        if self.liveUpdate.isChecked():
            if getmtime(osjoin(self.saveCurrentOutputDir,'particle_conc.sum')) != self.lastModTime:
                self.showParOutput(osjoin(self.saveCurrentOutputDir,'particle_conc.sum'),0)
                self.lastModTime = getmtime(osjoin(self.saveCurrentOutputDir,'particle_conc.sum'))


    def seeInAction(self):
        nb = self.n_bins_particle.value()
        try:
            x0 = float(self.min_particle_diam.text().replace('d','e'))
            x1 = float(self.max_particle_diam.text().replace('d','e'))
            self.splot(self.mmodal_input.text(),self.n_modal.text(),nb,x0,x1)
        except:
            pass


    def moveOneDay(self, days):
        day = self.dateEdit.date()
        day=day.addDays(days)
        self.dateEdit.setDate(day)
        self.showParOutput('load current',1)


    def show_currentInit(self,file):
        self.saveCurrentButton.setEnabled(True)
        self.saveCurrentButton.setStyleSheet("background-image: url('%s'); background-repeat: no-repeat" %(gui_path+"icons/%s/save.png"%icondir))
        self.actionSave_to_current.setEnabled(True)
        self.currentInitFile.setText(file)
        self.setWindowTitle(CurrentVersion+': '+file)
        self.saveCurrentButton.setToolTip('Save to '+ossplit(file)[1])
        self.currentInitFileToSave = file


    def inputPopup(self, v):
        """Envoke Question for main loop delay. (or any other question)"""
        exec('iii = '+v, locals(), globals())
        self.inwin = Input(default=iii)
        response = self.inwin.exec()
        if response != 0:
            exec("%s = %d" %(v,self.inwin.inp.input.value()))


    def updateEnvPath(self):
        if self.fileLoadOngoing:
            return
        status = self.update_nml()
        self.env_file.setToolTip('Location: "'+nml.ENV.ENV_FILE+'"')
        self.mcm_file.setToolTip('Location: "'+nml.MCM.MCM_FILE+'"')
        self.dmps_file.setToolTip('Location: "'+nml.PARTICLE.DMPS_FILE+'"')
        self.extra_particles.setToolTip('Location: "'+nml.PARTICLE.EXTRA_PARTICLES+'"')


    def get_case_kwargs(self,r):
        return {'begin':r[0],'end':r[1],'case':nml.PATH.CASE_NAME,'run':nml.PATH.RUN_NAME, 'common_root':nml.PATH.INOUT_DIR}


    def updatePath(self):
        if self.fileLoadOngoing:
            return
        status = self.update_nml()

        if self.indexRadioDate.isChecked():
            r = [self.dateEdit.text()]*2
        else:
            r = [self.indexEdit.value()]*2
        kwargs = self.get_case_kwargs(r)
        casedir = batch.batch(**kwargs)
        if len(casedir) == 8:
            # print(osjoin(casedir[-2]), osjoin(casedir[-1])) # xxx
            self.currentAddressTb.setText(casedir[-2]+'/')
            self.indir = casedir[-1]+'/'
        else:
            self.indir = '<Common root does not exist>/'
            self.currentAddressTb.setText(casedir)
        self.updateEnvPath()


    def ListbatchCaller(self,path):
        out = []
        with open(path, 'r') as f:
            for line in f:
                if line.strip() != '':
                    out.append(line.strip())
        self.batchCaller(listofinds=out)
        return


    def batchCaller(self, listofinds=[]):
        status = self.update_nml()
        if status != 0: return
        if self.batchRangeDay.isChecked():
            r = self.batchRangeDayBegin.text(),self.batchRangeDayEnd.text()
        else:
            r = self.batchRangeIndBegin.value(),self.batchRangeIndEnd.value()

        kwargs = self.get_case_kwargs(r)
        if len(listofinds)>0:
            kwargs['caselist'] = listofinds
        ret = batch.batch(**kwargs)
        if len(ret) == 8:
            dirs_to_create, conflicting_names, files_to_create, files_to_overwrite, existing_runs, dates,_,_ = ret
        else:
            self.popup('Error', ret, icon=2)
            return

        if conflicting_names != []:
            text = ''.join(['- %s\n' %i for i in conflicting_names])
            self.popup('File name conflicts',
            'The following (file)names are already in use:\n\n'+text
            +'\nPlease resolve the naming conflicts and try again.\n')
            return

        text_a = [[0,0,1],['','',''],['','','']]
        text_a[1][2] = 'The following files will be created:'

        if dirs_to_create != []:
            text_a[0][0] = 1
            text_a[1][0] = 'The following directories will be created:'
            text_a[2][0] = ''.join(['- %s\n' %i for i in dirs_to_create])

        if existing_runs != []:
            text_a[0][1] = 1
            text_a[1][1] = 'NOTE: There are some non-empty run directories and the output files would be overwritten by the Model:'
            text_a[2][1] = '\n'.join(['- %s' %i for i in existing_runs])
        # self.w.settext(0, text)

        if files_to_overwrite != []:
            for i,f in enumerate(files_to_create):
                if f in files_to_overwrite: files_to_create[i] = files_to_create[i]+' (current file will be overwritten)'

        text_a[2][2] = ''.join(['- %s\n' %i for i in files_to_create])
        if self.createBashFile.isChecked():
            bashfile = '%s/%s_%s_%s-%s.bash'%(kwargs['common_root'],kwargs['case'],kwargs['run'],dates[0],dates[-1])
            # If there are unnecessary trailing slashes, here will be double slashes, this should remove them
            bashfile = bashfile.replace('//','/')
            text_a[2][2] = text_a[2][2] + '\n- ' + bashfile
            if batch.paths(bashfile) == 1:
                text_a[2][2] = text_a[2][2] + ' (current file will be overwritten)'
        self.w = batchW(n=sum(text_a[0]))
        self.w.settext(text_a)
        response = self.w.exec()
        if response == 0:
            return

        dirs_to_create, _, files_to_create, _, _, dates,_,_= batch.batch(**kwargs)
        for path in dirs_to_create:
            mkdir(path)

        if self.createBashFile.isChecked():
            bf = open(bashfile, 'w')
            bf.write('#!/bin/bash\n')
            bf.write('# This file should be in %s\n'%kwargs['common_root'])
            bf.write('cd '+osrelpath(getcwd(), kwargs['common_root'] ) )
            bf.write('\n')

        for date,file in zip(dates,files_to_create):
            self.index_for_parser = date
            if self.createBashFile.isChecked():
                bf.write('./'+exe_name+' '+file+' |tee '+dirname(dirname(file)[:-1])+'/'+nml.PATH.RUN_NAME+'/runReport.txt'+'\n' )
            if self.batchRangeDay.isChecked():
                nml.TIME.DATE='%s'%(date)
                nml.TIME.NUMBER=''
            else:
                nml.TIME.DATE=''
                nml.TIME.NUMBER='%s'%(date)

            nml.ENV.ENV_FILE = self.pars(self.env_file.text(), file, self.stripRoot_env.isChecked())
            nml.MCM.MCM_FILE = self.pars(self.mcm_file.text(), file, self.stripRoot_mcm.isChecked())
            nml.ENV.LOSSES_FILE = self.pars(self.losses_file.text(), file, False)
            nml.PARTICLE.DMPS_FILE = self.pars(self.dmps_file.text(), file, self.stripRoot_par.isChecked())
            nml.PARTICLE.EXTRA_PARTICLES = self.pars(self.extra_particles.text(), file, self.stripRoot_xtr.isChecked())
            # nml.VAP.VAP_PROPS = self.pars(self.vap_props.text(), file, False)

            self.print_values(file=file,mode='silent',nobatch=False)

        if self.createBashFile.isChecked():
            bf.close()
            try:
                chmod(bashfile, 0o755)
            except:
                pass
        return


    def taghandler(self, tag):
        return (batch.tagparser(tag.group(0), self.index_for_parser))


    def pars(self, address, file=None, stripRoot=True):
        if address != '':
            if '>' in address and '<' in address:
                address = sub('<[iydm]*>', self.taghandler, address, flags=IGNORECASE)
            # else:
            if stripRoot:
                return (file[:file.rfind('/')+1] + address[address.rfind('/')+1:])
            return address
        else:
            return ''


    def resetSlider(self, slider, pos):
        slider.setProperty("value", pos)


    def scales(self):
        xf=1
        rt = self.runtime.value()
        scf = 24
        wScale = scf/2/200.0*xf
        pScale = scf*1.1905/200.0*xf
        aScale = 0.02*xf
        phScale = scf/0.4/200.0*xf
        ampScale = 1/20.0*xf
        return wScale,pScale,aScale,phScale,ampScale,rt


    def updteGraph(self, label='None', first=False):
        wScale,pScale,aScale,phScale,ampScale,rt = self.scales()

        x = linspace(0,rt,500)
        yscale = self.radio(self.fLin, self.fLog)

        dummy.sig = self.fWidth.value()*wScale
        if abs(dummy.sig)<0.001: dummy.sig = 0.001
        try:
            dummy.min = float(self.fMin.text())
        except:
            dummy.min = 0
        try:
            dummy.max = float(self.fMax.text())
        except:
            dummy.max=0
        dummy.mju = self.fPeak.value()*pScale
        dummy.fv = self.fFreq.value()*aScale
        dummy.ph = self.fPhase.value()*phScale
        dummy.am = self.fAmp.value()*ampScale
        dummy.sliderVls[0] = self.fWidth.value()
        dummy.sliderVls[1] = self.fPeak.value()
        dummy.sliderVls[2] = self.fFreq.value()
        dummy.sliderVls[3] = self.fPhase.value()
        dummy.sliderVls[4] = self.fAmp.value()

        if yscale == 'lin':
            dummy.mode = 1
        else:
            dummy.mode = 2

        self.monW.setValue(dummy.sig)
        self.monP.setValue(dummy.mju)
        self.monA.setValue(dummy.fv)
        self.monPh.setValue(dummy.ph)
        self.monAm.setValue(dummy.am)

        norm = self.gauss(dummy,yscale,rt)
        if float('.'.join(pg.__version__.split('.')[0:2]))<0.13:
            if (first or self.second) and self.legend in self.skene.items():
                self.skene.removeItem(self.legend)
                self.legend = self.PLOT.addLegend()
                if self.legend not in self.skene.items():
                    self.skene.addItem(self.legend)

        if first:
            self.currentPIVar = label
            self.editableselfPI = self.PLOT.plot(x,norm,pen=pg.mkPen('r', width=4), clear=True, name=label)
            self.second = True
        elif self.second:
            if not ' (mod)' in self.currentPIVar: self.currentPIVar = self.currentPIVar+' (mod)'
            self.editableselfPI = self.PLOT.plot(x,norm,pen=pg.mkPen('r', width=4), clear=True, name=self.currentPIVar)
        else:
            self.editableselfPI.setData(x,norm)

        if self.show_extra_plots != '':
            y=self.gauss(vars.mods[self.show_extra_plots], yscale,rt)
            if first or self.second:
                self.editableselfPI2 = self.PLOT.plot(x,y,pen=pg.mkPen(color='k', width=3,style=QtCore.Qt.DotLine), name=self.show_extra_plots)
            else:
                self.editableselfPI2.setData(x,y)
        if not first and self.second: self.second = False


    def saveParamValues(self):
        compound = self.names_sel_2.selectedItems()
        if compound != []:
            target = compound[0].text()
            for i in range(len(slMxs)):
                vars.mods[target].sl_x[i] = self.scaleFs[i].value()
            vars.mods[target].min  = dummy.min
            vars.mods[target].max  = dummy.max
            vars.mods[target].sig  = dummy.sig
            vars.mods[target].mju  = dummy.mju
            vars.mods[target].fv   = dummy.fv
            vars.mods[target].ph   = dummy.ph
            vars.mods[target].am   = dummy.am
            vars.mods[target].mode = dummy.mode
            for jj in range(len(vars.mods[target].sliderVls)):
                vars.mods[target].sliderVls[jj] = dummy.sliderVls[jj]
            vars.mods[target].pmInUse = 'Yes'
            confirmText = 'Saved to ' + target
            self.confirm.setText(confirmText)
            self.confirm.show()
            QtCore.QTimer.singleShot(2000, lambda : self.confirm.hide())
            self.updteGraph(label=target, first=True)
            for i in range(self.selected_vars.rowCount()):
                if self.selected_vars.item(i,0).text() == target:
                    self.selected_vars.cellWidget(i, 4).setCurrentIndex(1)
                    break


    def loadParamValues(self):
        compound = self.names_sel_2.selectedItems()
        if compound != []:
            target = compound[0].text()
            for i in range(len(slMxs)):
                self.scaleFs[i].setValue(vars.mods[target].sl_x[i])
            self.fMin.setText(str(vars.mods[target].min))
            self.fMax.setText(str(vars.mods[target].max))
            # self.gain.setProperty("value", vars.mods[target].gain)
            self.fWidth.setProperty("value", vars.mods[target].sliderVls[0])
            self.fPeak.setProperty("value", vars.mods[target].sliderVls[1])
            self.fFreq.setProperty("value", vars.mods[target].sliderVls[2])
            self.fPhase.setProperty("value", vars.mods[target].sliderVls[3])
            self.fAmp.setProperty("value", vars.mods[target].sliderVls[4])
            if vars.mods[target].mode < 2:
                self.fLin.setChecked(True)
                self.fLog.setChecked(False)
            else:
                self.fLin.setChecked(False)
                self.fLog.setChecked(True)
            self.updteGraph(label = target, first=True)
            confirmText = 'Loaded parameters from ' + target
            self.confirm.setText(confirmText)
            self.confirm.show()
            QtCore.QTimer.singleShot(2000, lambda : self.confirm.hide())


    def gauss(self,comp,ysc,rt):
        x = linspace(0,rt,500)
        mini = comp.min
        maxi = comp.max
        sigma = comp.sig
        peak = comp.mju
        freq = comp.fv
        phase = comp.ph
        amp = comp.am

        f = 1/sqrt(2*pi*sigma**2)
        D = peak + sin((x-peak)*freq)*amp + phase
        norm = 1/sqrt(2*pi*sigma**2)*exp(-(x-D)**2/(2*sigma**2))

        if ysc == 'lin':
            f = (maxi-mini)/f
            norm = norm*f + mini
        else:
            f = (log10(maxi-mini+1))/f
            norm = 10**(norm*f)-1 + mini
        return norm
        # D = peak + sin((x-peak)*freq)*amp + phase
        # norm = exp(-(x-D)**2/(2*sigma**2))
        #
        # if ysc == 'lin':
        #     h = (maxi-mini)
        #     norm = (norm-norm.min())/max(norm-norm.min())*h+mini
        # else:
        #     h = (log10(maxi-mini+1))
        #     norm = 10**((norm-norm.min())/max(norm-norm.min())*h)-1 + mini
        #
        # return norm


    def radio(self,*buts, action=True):
        for but in buts:
            if but.isChecked():
                if but.text() == 'Linear':
                    if action: self.PLOT.setLogMode(y=False)
                    return 'lin'
                if but.text() == 'Logarithmic':
                    if action: self.PLOT.setLogMode(y=True)
                    return 'log'

    def createCaseFolders(self,mode=0):
        cd = ''
        created = ''
        self.updatePath()
        relpath = self.currentAddressTb.text().split('/')
        for i,p in enumerate(relpath):
            if 'Common root does not exist' in p:
                self.popup('Error', 'Common root does not exist. Change "Common out" or create the necessary path: "'+self.inout_dir.text()+'"', icon=2)
                return
            if i>0:
                cd = osjoin(*relpath[:i])
                if cd != '' and not exists(cd):
                    mkdir(cd)
                    created = created  + cd +'\n'
        if created == '': created = 'No new directories created (all necessary directories existed already).'
        if mode == 0: self.popup('Created directories', created, icon=1)
        return

    def browse_path(self, target, mode, ftype=None, plWind=0, cmd = ''):
        """Browse for file or folder (depending on 'mode' and write the outcome to 'target')"""
        if mode == 'export':
            if 'Common root does not exist' in self.currentAddressTb.text():
                self.popup('Error', 'Common root does not exist. Change "Common out" or create the necessary path: "'+self.inout_dir.text()+'"', icon=2)
                return

        dialog = QtWidgets.QFileDialog()
        if ftype != None:
            dialog.setNameFilter(ftype)

        # if 'dir' in mode:
            # dialog.setFileMode(QtWidgets.QFileDialog.Directory)
            # dialog.setOption(QtWidgets.QFileDialog.ShowDirsOnly)
        options = dialog.Options()
        options |= dialog.DontUseNativeDialog
        if mode == 'dir':
            dialog.setFileMode(QtWidgets.QFileDialog.Directory)
            dialog.setOption(QtWidgets.QFileDialog.ShowDirsOnly)

            # path = dialog.getExistingDirectory(self, 'Choose Directory', options=options)
            dialog.setOption(QtWidgets.QFileDialog.DontUseNativeDialog)
            dialog.setWindowTitle('Choose directory')
            if dialog.exec() == 1:
                path = dialog.selectedFiles()[0]
            else: path=''

        elif mode == 'export':
            path = dialog.getSaveFileName(self, 'Save INITFILE', options=options)[0]
        else:
            dialog.setOption(QtWidgets.QFileDialog.DontUseNativeDialog)
            dialog.setWindowTitle('Choose File')
            if dialog.exec() == 1:
                path = dialog.selectedFiles()[0]
            else: path=''
        if path != '':
            if osrelpath(path, currentdir)[0] != '.':
            # if path[:currentdir_l] == currentdir and path[currentdir_l] == '/':
                path = osrelpath(path, currentdir)
            if mode == 'load':
                self.load_initfile(path)
                self.show_currentInit(path)
            elif mode == 'fixed':
                self.loadFixedFile(path)
            elif mode == 'fixed_cc':
                return path
            elif mode == 'plot':
                self.showOutput(path)
            elif mode == 'addplot':
                self.linePlotMulti(path)
            elif mode == 'addplot_more':
                self.linePlotMulti(path, new=False)
            elif mode == 'plot_mass':
                self.showMass(path)
            elif mode == 'plot_mass_add':
                self.showMass(path, add=True)
            elif mode == 'plotPar':
                self.showParOutput(path,plWind)
            elif mode == 'batchList':
                self.ListbatchCaller(path)
            # elif mode == 'dironly':
            #     path = path[path.rfind('/')+1:]
            #     target.clear()
            #     target.insert(path)
            elif mode == 'export':
                self.exportCurrentCase(path)
            elif mode == 'append':
                target.appendPlainText(path)
            else:
                target.clear()
                target.insert(path.replace('\\','/'))


    def save_file(self, file=None, mode=None,nobatch=True, changeTitle=True):
        if nobatch:
            status = self.update_nml()
            if status != 0: return
        if file==None:
            dialog = QtWidgets.QFileDialog()
            options = dialog.Options()
            options |= dialog.DontUseNativeDialog
            file = dialog.getSaveFileName(self, 'Save INITFILE', options=options)[0]
        if file != '':
            if file[:currentdir_l] == currentdir:
                file = file[currentdir_l+1:]
            self.print_values(file, mode)
            if file != defaults_file_path and changeTitle:
                self.show_currentInit(file)


    def markReverseSelection(self, op):
        """marks all variables or inverts selection"""
        for i in range(self.selected_vars.rowCount()):
            if i>1:
                if op == 'inv':
                    status = self.selected_vars.cellWidget(i,6).isChecked()
                    self.selected_vars.cellWidget(i,6).setChecked(not status)
                if op == 'all':
                    self.selected_vars.cellWidget(i,6).setChecked(True)


    def remv_item(self):
        """removes items from variable table"""
        # self.selected_vars.setSortingEnabled(False)
        for i in reversed(range(self.selected_vars.rowCount())):
            if self.selected_vars.cellWidget(i,6).isChecked():
                name = self.selected_vars.item(i,0).text()
                item = self.namesdat.item(namesPyInds[name])
                item.setFlags(item.flags() | QtCore.Qt.ItemIsEnabled | QtCore.Qt.ItemIsSelectable)
                item.setSelected(False)
                self.selected_vars.removeRow(i)
                vars.mods.pop(name)
                namesPyInds['PRESSURE']
        self.selected_vars.sortItems(7, QtCore.Qt.AscendingOrder)
        # self.selected_vars.setSortingEnabled(False)
        self.updateOtherTabs()


    def updateOtherTabs(self):
        """updates variable lists in other tabs"""
        self.names_sel_2.clear()
        for i in range(self.selected_vars.rowCount()):
            self.names_sel_2.addItem(self.selected_vars.item(i,0).text())


    def get_available_chemistry(self, checkonly=False):
        with open('makefile','r') as mk:
            for line in mk:
                if 'CHMDIR' in line and not '$' in line:
                    a,dir = line.replace('=','').split()
        if checkonly:
            return dir
        self.chemistryModules.clear()
        for (dirpath, dirnames, filenames) in walk('src/chemistry'):
            break
        dirnames = sorted(dirnames,key=str.lower)
        ii = -1
        for i,d in enumerate(dirnames):
            self.chemistryModules.addItem(d)
            if d == dir:
                if checkonly:
                    return dir
                else:
                    ii=i
        self.chemistryModules.setCurrentIndex(ii)
        if ii==-1:
            self.popup('Possible problems with makefile', '''
Chemistry scheme in "makefile" does not match with any
scheme names available in source. Please check your
makefile for variable "CHMDIR". You can also assign
a chemistry module in tab "Chemistry"''', icon=2)


    def editMakefile(self,mod):
        """Reads in the makefile, searches for line that contains Chemistry module directory name
        and replaces it with mod (which comes from current chemistry module drop down menu), replaces
        then the makefile.
        """
        replacement = 'CHMDIR = '+mod
        f = open('makefile','r')
        data = f.read()
        f.close()
        pattern = r'(CHMDIR)( )*(=)( )*(\S+)\n'
        data = sub(pattern,replacement, data)
        f = open('makefile','w')
        f.write(data)
        f.close()


    def loadFixedFile(self, path, cc=False):
        indef = False
        count = 0
        parsethis = False
        with open(path, 'r') as f:
            for line in f:
                if parsethis:
                    comps = line.upper().strip(' \n').replace('\t','').replace(' ','').split(',')
                    if cc: return comps
                    for comp in comps:
                        if comp not in vars.mods:
                            self.namesdat.item(namesPyInds[comp]).setSelected(True)
                            count = count +1
                    break
                if 'MOLECULAR WEIGHTS FOR SPECIES PRESENT IN THE SUBSET' in line.upper():
                    parsethis = True
        self.popup('File parsed', 'Selected %d variables'%count, icon=1)

    def vapours(self):
        """Envoke script to create new vapor file."""
        self.Wwin = VpressWin()
        response = self.Wwin.exec()
        if response == 0:
            return


    def createCC(self):
        """Envoke script to create new chemistry."""
        self.ccwin = CCWin()
        response = self.ccwin.exec()
        if response == 0:
            return


    def createVAR(self):
        """Envoke script to create Variations."""
        self.varWin = Variation()
        response = self.varWin.exec()
        if response == 0:
            return

    def createAb(self):
        """Envoke script to show About ARCA."""
        self.abwin = About()
        response = self.abwin.exec()
        if response == 0:
            return

    def editTxtFile(self, file):
        """Envoke script to edit some plain text file."""
        self.edwin = Editor(file=file,cursorToStart=True)
        response = self.edwin.exec()
        if response == 0:
            return
        else:
            t = self.edwin.editor.editedText.toPlainText()
        with open(file, 'w') as f:
            f.write(t)

    def loadFixedFromChemistry(self):
        '''Opens the chemistry fortran file and searches for fixed variables and selects them from the available vars'''
        chemistry = self.chemistryModules.currentText()
        try:
            count = 0
            with open(osjoin('src/chemistry/',chemistry,'second_Parameters.f90'),'r') as f:
                for line in f:
                    ind = line.upper().find('INDF_')
                    if ind > 0 and not '!' in line:
                        ind2 = line.find('=')
                        comp = line[ind+5:ind2].strip()
                        if comp not in vars.mods:
                            self.namesdat.item(namesPyInds[comp]).setSelected(True)
                            count += 1
                self.popup('File parsed', 'Selected %d variables'%count, icon=1)
        except:
            self.popup('Error','No \'second_Parameters.f90\' file\nin current chemistry.',icon=2)
        pass


    def highlightModifications(self, i):
        """For Multiply and Shift (in input variables tab), set font to bold if default values are changed"""
        row = i.row()
        c = i.column()
        if c==2 or c==3:
            try:
                if float(self.selected_vars.item(row,c).text())+c == 3.:
                    self.selected_vars.item(row,c).setFont(roman)
                else:
                    self.selected_vars.item(row,c).setFont(bold)
            except:
                self.selected_vars.item(row,c).setFont(bold)
        if c==1:
            try:
                float(self.selected_vars.item(row,c).text())
                self.selected_vars.item(row,c).setFont(roman)
                vars.mods[self.selected_vars.item(row,0).text()].tied = ''
            except:
                xx=self.selected_vars.item(row,1).text()
                vars.mods[self.selected_vars.item(row,0).text()].tied = xx
                self.selected_vars.item(row,c).setFont(bold)

    def add_new_line(self, name, unit_ind, cols=[],createNew=True, unt=0):
        """adds items to variable table"""
        if createNew:
            vars.mods[name] = Comp()
            vars.mods[name].Find = namesFoInds[name]
            vars.mods[name].name = name # Human readable name for modified variable
        # self.selected_vars.setSortingEnabled(False);
        row = self.selected_vars.rowCount()
        self.selected_vars.insertRow(row)
        item = self.namesdat.item(namesPyInds[name])
        item.setFlags(item.flags() & ~QtCore.Qt.ItemIsEnabled & ~QtCore.Qt.ItemIsSelectable)

        pmInUse = QtWidgets.QComboBox()
        pmInUse.addItems(['No','Yes'])
        unit = QtWidgets.QComboBox()
        unit.addItems(units.get(grepunit(name),units['REST']))
        unit.setCurrentIndex(unt)
        markBut = QtWidgets.QPushButton()
        markBut.setCheckable(True)
        if name == 'TEMPK' or name == 'PRESSURE' :
            unit.setCurrentIndex(1)
            markBut.setEnabled(False)
            markBut.setToolTip("Temperature and pressure are always required and cannot be marked for removal")
        else:
            markBut.setToolTip("Mark variable for removal")

        markBut.setText('mark')
        if cols==[]:
            cols = [name, '-1','1.0', '0.0',0]
        self.selected_vars.horizontalHeader().setStretchLastSection(True)
        self.selected_vars.horizontalHeader().setSectionResizeMode(0, QtWidgets.QHeaderView.Stretch)

        for i in range(4):
            tag = QtWidgets.QTableWidgetItem(cols[i])
            if i==0: tag.setFlags(tag.flags() & ~QtCore.Qt.ItemIsEditable)
            self.selected_vars.setItem(row, i, tag)
            if namesPyInds[name]<divider_i:
                self.selected_vars.item(row, i).setBackground(QtGui.QColor(*env_no))
            elif namesPyInds[name]>divider_xtr_i:
                self.selected_vars.item(row, i).setBackground(QtGui.QColor(*xtr_no))
            else:
                self.selected_vars.item(row, i).setBackground(QtGui.QColor(*org_no))
        if name == 'PRESSURE' :
            self.selected_vars.setItem(row, 3, QtWidgets.QTableWidgetItem('1e3'))
        self.selected_vars.setCellWidget(row, i+1, pmInUse )
        pmInUse.currentIndexChanged.connect(lambda: self.toggleColorPre(name))
        self.selected_vars.itemChanged.connect(self.highlightModifications)
        pmInUse.setCurrentIndex(cols[4])
        self.selected_vars.setCellWidget(row, i+2, unit )
        self.selected_vars.setCellWidget(row, i+3, markBut )

        self.selected_vars.setItem(row, i+4, QtWidgets.QTableWidgetItem('%03d'%(namesFoInds[name])))
        self.selected_vars.sortItems(7, QtCore.Qt.AscendingOrder)
        # self.selected_vars.setSortingEnabled(False)
        self.updateOtherTabs()


    def toggleColorPre(self, n):
        if namesPyInds[n]<divider_i:
            c = (env_yes,env_no)
        elif namesPyInds[n]>divider_xtr_i:
            c = (xtr_yes,xtr_no)
        else:
            c = (org_yes,org_no)
        for i in range(self.selected_vars.rowCount()):
            if self.selected_vars.item(i,0).text() == n: break
        self.toggleColor(i,c)


    def toggleColor(self,r,c):
        if self.selected_vars.cellWidget(r,4).currentText() == 'Yes':
            for z in range(4):
                self.selected_vars.item(r, z).setBackground(QtGui.QColor(*c[0]))
        else:
            for z in range(4):
                self.selected_vars.item(r, z).setBackground(QtGui.QColor(*c[1]))


    def toggle_frame(self, frame):
        if frame.isEnabled() == True:
            frame.setEnabled(False)
        else:
            frame.setEnabled(True)


    def toggle_gray(self, guard, group):
        index = group.count()
        for i in range(index):
            grayWidget = group.itemAt(i).widget()
            if guard.isChecked() == True:
                grayWidget.setEnabled(True)
            else:
                grayWidget.setEnabled(False)


    def grayIfNotChecked(self, guard, frame, hide=-1):
        if hide>=0:
            if guard.isChecked() == True:
                self.tabWidget.setTabEnabled(hide,True)
            else:
                self.tabWidget.setTabEnabled(hide,False)
        else:
            if guard.isChecked() == True:
                frame.setEnabled(True)
            else:
                frame.setEnabled(False)


    def grayIfChecked(self, guard, frame):
        if 'bool' in str(type(guard)):
            tester=guard
        else:
            tester=guard.isChecked()
        if tester == True:
            frame.setChecked(0)
            frame.setEnabled(False)
        else:
            frame.setEnabled(True)


    def toggle_printtime(self):
        if self.fsave_division.value() != 0 :
            self.fsave_interval.setEnabled(False)
        else:
            self.fsave_interval.setEnabled(True)


    def print_values(self, file=None, mode=None, nobatch=True,status=0):
        if nobatch:
            status = self.update_nml()
            if status != 0 and file: return status
        self.prints += 1
        if (file):
            f = open(file, 'w')
            f.write('#'+('-')*50+'\n')
            f.write('#      ARCA box setting file: %s\n'%(file))
            f.write('#         Created at: '+( time.strftime("%b %d %Y, %H:%M:%S", time.localtime()))+'\n')
            f.write('#'+('-')*50+'\n\n')
            nml.printall(vars.mods, target='f',f=f)
            f.close()
            if not mode=='silent':
                if mode == 'noteOnTer':
                    print('Settings were saved to file "'+file+'"')
                    self.popup('Saved','Settings were saved to file "'+file+'"',0)
                    return
                if file == defaults_file_path:
                    self.popup('Saved', 'Defaults saved', icon=0)
                elif file != tempfile:
                    self.popup('Saved settings', 'Settings were saved to '+file, icon=0)

        else:
            print(('\n')*10, )
            print('#',('-')*50)
            print('#              ARCA box setting file #%d'%(self.prints))
            print('#         Created at:', ( time.strftime("%b %d %Y, %H:%M:%S", time.localtime())))
            print('#',('-')*50, '\n')
            nml.printall(vars.mods, target='p')

        return status

    def select_compounds(self):
        compounds = self.namesdat.selectedItems()
        for c in compounds:
            self.add_new_line(c.text(), 2)


    def select_compounds_for_plot(self):
        if self.plotTo.isChecked() == True:
            compound = self.names_sel_2.selectedItems()
            if compound != []:
                self.show_extra_plots = compound[0].text()
                self.updteGraph(first=True, label=self.currentPIVar)
            else:
                self.plotTo.setChecked(False)
        if self.plotTo.isChecked() == False:
            self.show_extra_plots = ''
            self.updteGraph(first=True, label=self.currentPIVar)


    def stopBox(self):
        self.Timer.stop()
        # self.TimerPlot.stop()
        self.pollTimer.stop()
        self.liveUpdate.setChecked(False)
        self.liveUpdate.setEnabled(False)
        self.boxProcess.kill()
        tout = self.boxProcess.wait(timeout=10)
        self.boxProcess.poll()
        self.toggle_frame(self.frameStop)
        self.toggle_frame(self.frameStart)
        self.MonitorWindow.setTextInteractionFlags(QtCore.Qt.TextSelectableByMouse)
        # try: self.tabWidget.currentChanged.disconnect(self.activeTab)
        # except: pass

        if exists(self.saveCurrentOutputDir):
            f = open(osjoin(self.saveCurrentOutputDir,'runReport.txt'), 'w')
            f.write(self.MonitorWindow.toPlainText())
            f.close()
        else:
            return

    def softStop(self):
        if exists(self.saveCurrentOutputDir):
            f = open(osjoin(self.saveCurrentOutputDir,'ENDNOW.INIT'), 'w')
            f.write('STOP')
            f.close()
            self.popup('Asked the model to finish', "Depending on the state of the model this could take a while. "
            +"If you don't need the NetCDF output you can also Force Stop.",1)


    def showParOutput(self, file, windowInd):
        if windowInd == 0:
            window = self.surfacePlotWindow_0
            titleLoc = self.parPlotTitle_0
        if windowInd == 1:
            window = self.surfacePlotWindow_1
            titleLoc = self.parPlotTitle_1
        levels=(self.lowlev.value(),self.highlev.value())
        if self.firstParPlot[windowInd] == 0 and scipyIs:
            self.gauss_x.valueChanged.connect(lambda: self.drawSurf(window))
            self.gauss_y.valueChanged.connect(lambda: self.drawSurf(window))
            self.Filter_0.clicked.connect(lambda: self.drawSurf(window))
            self.Filter_1.clicked.connect(lambda: self.drawSurf(window))
            self.lowlev.valueChanged.connect(lambda: self.drawSurf(window))
            self.highlev.valueChanged.connect(lambda: self.drawSurf(window))
            self.cmJet.triggered.connect(lambda: self.drawSurf(window))
            self.firstParPlot[windowInd] = 1
        if '.nc' in file[-4:] and not netcdf:
            self.popup(*netcdfMissinnMes)
            return

        if file == 'load current':
            file = self.pars(self.dmps_file.text(), file=self.indir, stripRoot=self.stripRoot_par.isChecked())
            if not exists(file):
                self.popup('Not found','File %s not found'%file, icon=2)
                return
        titleLoc.setText(file)


        window.clear()
        ret = par.loadNC(file, self.actionAssume_lognormal_input.isChecked())
        if 'str' in str(type(ret)):
            self.popup('Error', 'File was not completely ok, maybe an interrupted run? Error message: \n'+ret, icon=2)
            return
        time,diam,n = ret

        if windowInd==0:
            self.z0 = n
            if self.X_axis_in_time.isChecked():
                self.time0 = (time[0],time[-1])
            else:
                self.time0 = (0,n.shape[0])
            if self.Y_axis_in_nm.isChecked():
                self.diam0 = (log10(diam[0]*1e9), log10(diam[-1]*1e9), log10(diam[-2]*1e9))
                self.log0 = True
            else:
                self.diam0 = (0,n.shape[1], n.shape[1]-1)
                self.log0 = False

        if windowInd==1:
            self.z1 = n
            if self.X_axis_in_time.isChecked():
                self.time1 = (time[0],time[-1])
            else:
                self.time1 = (0,n.shape[0])
            if self.Y_axis_in_nm.isChecked():
                self.diam1 = (log10(diam[0]*1e9), log10(diam[-1]*1e9), log10(diam[-2]*1e9))
                self.log1 = True
            else:
                self.diam1 = (0,n.shape[1], n.shape[1]-1)
                self.log1 = False
        self.drawSurf(window, new=1)


    def parse_ACDC_systems(self, num=0):

        def parsefile(d):
            f = open(osjoin('src/ACDC',d,'acdc_system_0x%d.f90'%int(d[-1])))
            t = f.read()
            nickname = ''.join(re.findall('! System name: .+',t)).replace('! System name: ','')
            if num==0:
                self.nickname.append(nickname)
            else:
                self.nickname[num-1] = nickname
            x = re.findall('integer, parameter :: .+ neutral_monomers.+= ?\(/( ?.+)/\)',t)
            y = re.findall(' ?clust\((\w+)\)\(\:\) ?\= ?(.+)',t)
            names = {}
            for i in range(len(y)):
                names[y[i][0]] = y[i][1]
            monomers=[]
            for j in x[0].split():
                monomers.append(names[j.strip(',')].replace('1','').strip("'"))
            return monomers

        def linker(n,new=True):
            for im,m in enumerate(monomers):
                if self.ACDC_linker.topLevelItem(n).child(im) != None:
                #     new = True
                # else:
                    csystm = self.ACDC_linker.topLevelItem(n)
                    children = [csystm.child(nnn) for nnn in reversed(range(csystm.childCount()))]
                    for child in children:
                        csystm.removeChild(child)
                item_1 = QtWidgets.QTreeWidgetItem(self.ACDC_linker.topLevelItem(n))
                self.ACDC_linker.topLevelItem(n).child(im).setText(0, m)
                linkto = '-'
                self.ACDC_linker.topLevelItem(n).child(im).setText(1, linkto)
                self.ACDC_linker.topLevelItem(n).child(im).setFlags(QtCore.Qt.ItemIsEnabled | QtCore.Qt.ItemIsEditable)



        if num==0: cl_all = True
        if num>0:  cl_all = False
        if cl_all:
            self.ACDC_linker.clear()
            for (dirpath, dirnames, filenames) in walk('src/ACDC'):
                break
            n = 0
            dirnames.sort()
        if num>0:
            monomers = parsefile('ACDC_%02d'%num)
            linker(num-1,new=False)
            self.set_ACDC_links()
        else:
            for i,d in enumerate(dirnames):
                if not 'ACDC_0' in d:
                    continue
                else:
                    item_0 = QtWidgets.QTreeWidgetItem(self.ACDC_linker)
                    item_0.setFlags(QtCore.Qt.ItemIsEnabled | QtCore.Qt.ItemIsUserCheckable )
                    item_0.setCheckState(0, QtCore.Qt.Unchecked)
                    self.ACDC_linker.topLevelItem(n).setText(0, d)
                    monomers = parsefile(d)
                    # self.ACDC_linker.topLevelItem(n).setText(2, self.nickname[n])

                    linker(n)
                    n +=1

            self.ACDC_n_systems = n
        self.ACDC_current_links = self.get_ACDC_links()
        self.ACDC_available_compounds = self.get_ACDC_links()


    def get_ACDC_links(self):
        root = self.ACDC_linker.invisibleRootItem()
        child_count = root.childCount()
        acdc_links = []
        for i in range(child_count):
            item = root.child(i)
            if item.checkState(0): self.acdc_systems_flags[i] = 1
            else: self.acdc_systems_flags[i] = 0
            acdc_links.append([item.checkState(0)])
            child_count1 = item.childCount()
            for j in range(child_count1):
                if j==0: acdc_links[i].append({})
                item2 = item.child(j)
                acdc_links[i][1][item2.text(0)] = item2.text(1)
        return acdc_links


    def set_ACDC_links(self,loose=False):
        root = self.ACDC_linker.invisibleRootItem()

        child_count = root.childCount()
        for i in range(child_count):
            item = root.child(i)
            item.setText(2,self.nickname[i])
            # print(self.nickname[i])
            if self.acdc_systems_flags[i]>0:
                item.setCheckState(0, QtCore.Qt.Checked)
                item.setExpanded(True)
            else:
                item.setCheckState(0, QtCore.Qt.Unchecked)
            child_count1 = item.childCount()
            for j in range(child_count1):
                item2 = item.child(j)
                if item2.text(0) in self.ACDC_current_links[i][1] or (loose and j==1):
                    item2.setText(1,self.ACDC_current_links[i][1][item2.text(0)])
                else:
                    item2.setText(1,'-')


    def drawSurf(self,window, new=0):
        use_filter = False
        if window==self.surfacePlotWindow_0:
            n_levelled = self.z0
            time = self.time0
            diam = self.diam0
            l_scale = self.log0
            if self.Filter_0.isChecked():
                use_filter = True
        else:
            n_levelled = self.z1
            time = self.time1
            diam = self.diam1
            l_scale = self.log1
            if self.Filter_1.isChecked():
                use_filter = True
        filter_MA = False
        levels=(self.lowlev.value(),self.highlev.value())
        if scipyIs and use_filter:
            if filter_MA:
                w = (int(shape(n_levelled)[1]/10))//2*2+1
                for i in range(shape(n_levelled)[0]):
                    for j in range(5):
                        n_levelled[i,:] = savgol_filter(n_levelled[i,:],w,2)
                n_levelled = where(n_levelled<0,0,n_levelled)
            else:
                n_levelled = gaussian_filter(n_levelled,(self.gauss_x.value(),self.gauss_y.value()),mode='constant')
        n_levelled = where(n_levelled>=levels[1],levels[1]*0.98,n_levelled)
        hm = pg.ImageItem(n_levelled)
        # if self.Y_axis_in_nm.isChecked():
        #     bx = QtCore.QRectF(time[0], log10(diam[0]*1e9), time[1], log10(diam[1]*1e9))
        # else:
        bx = QtCore.QRectF(time[0], diam[0], time[1], diam[1]-diam[0] + diam[1]-diam[2])
        # (cb[-1,0]-cb[0,0] + cb[-1,0]-cb[-2,0])
        hm.setRect(bx)

        cb = ndarray((20,1))
        cb[:,0] = linspace(self.lowlev.value(), self.highlev.value(), 20)

        ss = pg.ImageItem(cb.T)
        ss.setRect(QtCore.QRectF(0,self.lowlev.value(),1,(cb[-1,0]-cb[0,0] + cb[-1,0]-cb[-2,0])))

        try:
            # If matplotlib is installed, we get colours
            from matplotlib import cm
            if self.cmJet.isChecked():
                colormap = cm.get_cmap("jet")
            else:
                colormap = cm.get_cmap("viridis")
            colormap._init()
            lut = (colormap._lut * 255).view(ndarray)  # Convert matplotlib colormap from 0-1 to 0 -255 for Qt
            # Apply the colormap
            hm.setLookupTable(lut)
            ss.setLookupTable(lut)
        except:
            pass

        hm.setLevels(levels)
        ss.setLevels(levels)
        window.setMenuEnabled(False)
        window.addItem(hm)
        self.cbWindow.clear()
        self.cbWindow.enableAutoRange()
        self.cbWindow.addItem(ss)
        self.cbWindow.hideAxis('bottom')
        self.cbWindow.setLogMode(False, True)
        if l_scale:
            window.setLogMode(False, True)
        else:
            window.setLogMode(False, False)

    def plotChanges(self):
        import modules.optich as optich
        ret = optich.plot(self.currentAddressTb.text())
        if ret!='': self.popup('Plot changes:',ret,0)

    def StartboxShortcut(self):
        if self.frameStop.isEnabled():
            self.stopBox()
        self.tabWidget.setCurrentIndex(6)
        self.startBox()

    def QuitShortcut(self):
        if self.frameStop.isEnabled():
            self.stopBox()
        self.tabWidget.setCurrentIndex(6)

    def QuitGracefullyShortcut(self):
        if self.frameStop.isEnabled():
            self.softStop()
        self.tabWidget.setCurrentIndex(6)

    def startBox(self):
        self.closenetcdf()
        self.closenetcdf_mass()
        self.createCaseFolders(mode=1)
        self.MonitorWindow.setTextInteractionFlags(QtCore.Qt.NoTextInteraction)
        self.pauseScroll.setChecked(False)

        currentWait = self.wait_for
        self.wait_for = 0
        if not exists(gui_path+'tmp'):
            mkdir(gui_path+'tmp')
        st = self.print_values(tempfile)
        if st != 0: return
        self.wait_for = currentWait

        try:
            self.boxProcess = Popen(["./"+exe_name, "%s"%tempfile, '--gui'], stdout=PIPE,stderr=STDOUT,stdin=None)
            self.saveCurrentOutputDir = self.currentAddressTb.text()
            self.MonitorWindow.clear()
            self.Timer.start(5)
            # self.TimerPlot.start(1000)
            self.pollTimer.start(1)
            QtCore.QTimer.singleShot(2000, lambda : self.MonitorWindow.verticalScrollBar().setSliderPosition(self.MonitorWindow.verticalScrollBar().maximum()))
            self.liveUpdate.setEnabled(True)
            self.toggle_frame(self.frameStart)
            self.toggle_frame(self.frameStop)
        except:
            self.MonitorWindow.appendPlainText('\n    Could not start model executable, is it compiled?')


    def pollMonitor(self):
        self.fulltext = self.boxProcess.stdout.readline().decode("utf-8")
        if self.fulltext != '.\r\n' and self.fulltext != '.\n':
            self.MonitorWindow.appendPlainText(self.fulltext.strip('\n'))
        self.monStatus = self.boxProcess.poll()
        if self.monStatus != None and self.fulltext == '':
            self.stopBox()


    def updateOutput(self):
        self.fulltext = self.boxProcess.stdout.readline().decode("utf-8")
        if self.fulltext != '.\r\n' and self.fulltext != '.\n' and self.fulltext != '':
            self.MonitorWindow.appendPlainText(self.fulltext.strip('\n'))
        # if self.pauseScroll.isChecked() == False:
        #     self.MonitorWindow.verticalScrollBar().setSliderPosition(self.MonitorWindow.verticalScrollBar().maximum())
        # if 'SIMULATION HAS ENDED' in str(self.fulltext)[-50:]:
        #     # self.MonitorWindow.setPlainText(self.MonitorWindow.toPlainText())
        #     self.MonitorWindow.verticalScrollBar().setSliderPosition(self.MonitorWindow.verticalScrollBar().maximum())


    def checkboxToFOR(self, widget):
        if widget.isChecked() == True:
            return '.TRUE.'
        else:
            return '.FALSE.'


    def checkboxToINT(self, widget):
        if widget.isChecked() == True:
            return 1
        else:
            return 0


    def update_nml(self):
        status = 0
        self.ACDC_current_links = self.get_ACDC_links()
        # class _SETTINGS:
        nml.SETTINGS.BATCH = '%s:%s %s:%s %s:%s %s:%s %s:%s %s:%s %s:%s %s:%s %s:%s' %(
                                'batchRangeDayBegin', self.batchRangeDayBegin.text(),
                                'batchRangeDayEnd', self.batchRangeDayEnd.text(),
                                'batchRangeIndBegin', self.batchRangeIndBegin.text(),
                                'batchRangeIndEnd', self.batchRangeIndEnd.text(),
                                'indexRadioIndex', self.checkboxToFOR(self.indexRadioIndex),
                                'indexRadioDate', self.checkboxToFOR(self.indexRadioDate),
                                'createBashFile', self.checkboxToFOR(self.createBashFile),
                                'batchRangeDay', self.checkboxToFOR(self.batchRangeDay),
                                'batchRangeInd', self.checkboxToFOR(self.batchRangeInd),
        )
        nml.SETTINGS.INPUT = '%s:%s %s:%s %s:%s %s:%s %s:%s %s:%s %s:%s %s:%s %s:%s %s:%s' %(
                                'env_file', self.env_file.text(),
                                'mcm_file', self.mcm_file.text(),
                                'dmps_file', self.dmps_file.text(),
                                'extra_particles', self.extra_particles.text(),
                                'losses_file', self.losses_file.text(),
                                'spectralFunctions', self.spectralFunctions.text(),
                                'stripRoot_env', self.checkboxToFOR(self.stripRoot_env),
                                'stripRoot_mcm', self.checkboxToFOR(self.stripRoot_mcm),
                                'stripRoot_par', self.checkboxToFOR(self.stripRoot_par),
                                'stripRoot_xtr', self.checkboxToFOR(self.stripRoot_xtr),
        )

        # class NAMES:
        nml.NAMES.NAMESDAT=path_to_names
        # class _PATH:
        nml.PATH.INOUT_DIR=self.inout_dir.text()
        if nml.PATH.INOUT_DIR == '': nml.PATH.INOUT_DIR = default_inout
        nml.PATH.CASE_NAME=self.case_name.text().upper()
        if nml.PATH.CASE_NAME == '': nml.PATH.CASE_NAME = default_case
        nml.PATH.RUN_NAME = self.run_name.text().upper()
        if nml.PATH.RUN_NAME == '': nml.PATH.RUN_NAME = default_run

        # class _FLAG:
        nml.FLAG.CHEMISTRY_FLAG=self.checkboxToFOR(self.checkBox_che)
        nml.FLAG.AEROSOL_FLAG=self.checkboxToFOR(self.checkBox_aer)
        nml.FLAG.ACDC_SOLVE_SS=self.checkboxToFOR(self.acdc_solve_ss)
        nml.FLAG.ACDC=self.checkboxToFOR(self.checkBox_acd)
        nml.FLAG.CONDENSATION=self.checkboxToFOR(self.condensation)
        nml.FLAG.COAGULATION=self.checkboxToFOR(self.coagulation)
        nml.FLAG.DEPOSITION=self.checkboxToFOR(self.deposition)
        nml.FLAG.CHEM_DEPOSITION=self.checkboxToFOR(self.chemDeposition)
        # nml.FLAG.MODEL_H2SO4=self.checkboxToFOR(self.model_h2so4)
        nml.FLAG.ORG_NUCL=self.checkboxToFOR(self.Org_nucl)
        nml.FLAG.PRINT_ACDC=self.checkboxToFOR(self.print_acdc)
        nml.FLAG.OPTIMIZE_DT=self.checkboxToFOR(self.useSpeed)
        nml.FLAG.AFTER_CHEM_ON=self.checkboxToFOR(self.after_chem_on)
        nml.FLAG.AFTER_NUCL_ON=self.checkboxToFOR(self.after_nucl_on)
        nml.FLAG.FILE_TIME_UNIT=self.fileTimeUnit_a.currentText()
        nml.FLAG.LOSSFILE_TIME_UNIT=self.lossfileTimeUnit.currentText()

        # class _TIME:
        nml.TIME.RUNTIME=self.runtime.value()
        nml.TIME.DT=self.dt.value()
        nml.TIME.FSAVE_INTERVAL=self.fsave_interval.value()
        nml.TIME.PRINT_INTERVAL=self.print_interval.value()
        nml.TIME.FSAVE_DIVISION=self.fsave_division.value()
        if self.indexRadioDate.isChecked():
            self.index_for_parser = self.dateEdit.text()
            nml.TIME.DATE=self.dateEdit.text()
            nml.TIME.NUMBER=''
        if self.indexRadioIndex.isChecked():
            self.index_for_parser = self.indexEdit.value()
            nml.TIME.DATE=''
            nml.TIME.NUMBER='%04d'%self.indexEdit.value()

        # class _PARTICLE:
        nml.PARTICLE.PSD_MODE=self.psd_mode.currentIndex()
        nml.PARTICLE.N_BINS_PAR=self.n_bins_particle.value()
        try:
            float(self.min_particle_diam.text().lower().replace('d','e'))
            float(self.max_particle_diam.text().lower().replace('d','e'))
        except:
            self.popup('Error', 'Minimum and maximum particle diameters must be numbers', icon=3)
            status = 1
        nml.PARTICLE.MIN_PARTICLE_DIAM=self.min_particle_diam.text()
        nml.PARTICLE.MAX_PARTICLE_DIAM=self.max_particle_diam.text()
        nml.PARTICLE.MMODAL_INPUT_INUSE=self.checkboxToINT(self.multiModalBox)
        nml.PARTICLE.N_MODAL=self.n_modal.text()
        nml.PARTICLE.DMPS_FILE=self.pars(self.dmps_file.text(), file=self.indir, stripRoot=self.stripRoot_par.isChecked())
        nml.PARTICLE.EXTRA_PARTICLES=self.pars(self.extra_particles.text(), file=self.indir, stripRoot=self.stripRoot_xtr.isChecked())
        nml.PARTICLE.MMODAL_INPUT=self.mmodal_input.text()
        nml.PARTICLE.DMPS_READ_IN_TIME=self.dmps_read_in_time.value()
        nml.PARTICLE.DMPS_HIGHBAND_LOWER_LIMIT=self.dmps_highband_lower_limit.text()
        nml.PARTICLE.DMPS_LOWBAND_UPPER_LIMIT=self.dmps_lowband_upper_limit.text()
        nml.PARTICLE.USE_DMPS=self.checkboxToFOR(self.use_dmps)
        nml.PARTICLE.USE_DMPS_PARTIAL=self.checkboxToFOR(self.use_dmps_partial)
        nml.PARTICLE.DMPS_INTERVAL=self.dmpsIntvl.value()

        # class _ENV:
        nml.ENV.ENV_FILE=self.pars(self.env_file.text(), file=self.indir, stripRoot=self.stripRoot_env.isChecked())
        nml.ENV.SPECTRUMFILE=self.pars(self.spectralFunctions.text(), file=self.indir, stripRoot=False)
        nml.ENV.LOSSES_FILE=self.pars(self.losses_file.text(), file=self.indir, stripRoot=False)
        nml.ENV.CHAMBER_FLOOR_AREA=self.floorArea.value()
        nml.ENV.SWR_IS_ACTINICFLUX=self.checkboxToFOR(self.SW_is_AF)
        # nml.ENV.CHAMBER_CIRCUMFENCE=self.chamberCircumfence.value()
        nml.ENV.CHAMBER_HEIGHT=self.chamberHeight.value()
        nml.ENV.EDDYK=self.eddyK.text()
        nml.ENV.USTAR=self.ustar.text()
        nml.ENV.ALPHAWALL=self.alphaWall.text()
        nml.ENV.CW_EQV=self.Cw_eqv.text()
        nml.ENV.SWR_IN_LOWER=self.swr_in_lower.value()
        nml.ENV.SWR_IN_UPPER=self.swr_in_upper.value()

        # class _MCM:
        nml.MCM.MCM_FILE=self.pars(self.mcm_file.text(), file=self.indir, stripRoot=self.stripRoot_mcm.isChecked())

        # class _MISC:
        nml.MISC.LAT=self.lat.value()
        nml.MISC.LON=self.lon.value()
        nml.MISC.WAIT_FOR=self.wait_for
        nml.MISC.DESCRIPTION=self.description.toPlainText().replace('\n','<br>')
        nml.MISC.CH_ALBEDO=self.ch_albedo.value()
        nml.MISC.GR_SIZES=self.GR_sizes.text()

        # class _VAP:
        nml.VAP.VAP_NAMES=self.vap_names.text()
        nml.VAP.USE_ATOMS=self.checkboxToFOR(self.use_atoms)
        nml.VAP.VAP_ATOMS=self.vap_atoms.text()

        # class _PRECISION:
        nml.PRECISION.DDIAM_RANGE="%f,%f"%(self.prec_low_1.value(), self.prec_high_1.value())
        nml.PRECISION.DPNUM_RANGE="%f,%f"%(self.prec_low_2.value(), self.prec_high_2.value())
        nml.PRECISION.DVAPO_RANGE="%f,%f"%(self.prec_low_3.value(), self.prec_high_3.value())

        for i in range(self.selected_vars.rowCount()):
            name = self.selected_vars.item(i,0).text()
            if vars.mods[name].tied != '':
                vars.mods[name].col = (self.selected_vars.item(i,1).text())
            else:
                vars.mods[name].col = int(self.selected_vars.item(i,1).text())
            vars.mods[name].multi = float(self.selected_vars.item(i,2).text())
            vars.mods[name].shift = float(self.selected_vars.item(i,3).text())
            vars.mods[name].pmInUse = self.selected_vars.cellWidget(i,4).currentText()
            vars.mods[name].unit = self.selected_vars.cellWidget(i,5).currentText()

        nml.RAW.RAW = self.rawEdit.toPlainText()
        nml.CUSTOM.CUSTOMS = []
        for i in range(1,31):
            key = 'customKey_%d'%i
            value = 'customVal_%d'%i
            exec('nml.CUSTOM.CUSTOMS.append([self.%s.text(),self.%s.text()])'%(key,value))

        nml.ACDC.ACDC_SYSTEMS = ','.join([str(i) for i in self.acdc_systems_flags])
        nml.ACDC.LINKS = [' '] * self.ACDC_n_systems
        for i in range(self.ACDC_n_systems):
            nml.ACDC.LINKS[i] = ' '.join(['%s %s'%(j,self.ACDC_current_links[i][1][j]) for j in self.ACDC_current_links[i][1].keys()])

        return status

    def resolveHelper(self):
        text = self.fill_formation_with.currentText()
        if text == 'Fixed ratio':
            return ''
        elif 'NH3' in text:
            return 'NH3'
        elif 'DMA' in text:
            return 'DMA'


    def load_initfile(self,file):
        INC_COMP = False
        timeunits = {'day':0,'hrs':1,'min':2,'sec':3}
        self.fileTimeUnit_a.setCurrentIndex( timeunits['day'] )
        self.fileTimeUnit_b.setCurrentIndex( timeunits['day'] )
        self.lossfileTimeUnit.setCurrentIndex( timeunits['day'] )

        if not exists(file):
            if file != defaults_file_path:
                self.popup('Ooops', 'File "%s" not found'%file, icon=3)
                return 0
        self.fileLoadOngoing = True
        self.markReverseSelection('all')
        self.remv_item()
        if self.plotTo.isChecked() == True:
            self.plotTo.setChecked(False)
            self.show_extra_plots = ''
            self.updteGraph(first=True)

        # def solve_for_parser(query):
        #     if query.upper() == 'NH3': return 2
        #     elif query.upper() == 'DMA': return 1
        #     else: return 0

        def parse_date(strg):
            str = strg.strip()
            if len(str)==10:
                y = int(str[0:4])
                m = int(str[5:7])
                d = int(str[8:10])
                return QtCore.QDate(y, m, d)
            else:
                return QtCore.QDate(2000, 1, 1)


        try:
            f = open(file, 'r')
            for line in f:
                break
            f.close()
        except:
            self.popup('Not a valid file', 'You are trying to open a file which does not appear to be a valid ARCA configuration file')
            return

        f = open(file, 'r')
        in_custom = False

        for line in f:
            i = line.find('=')
            x = line.find('!')
            if i<x or (x==-1 and i!=-1):
                key = line[:i].upper().strip()
                if key == '# RAW_INPUT':
                    rawline = line[i+1:].lstrip().rstrip()
                # remove comma and excess whitespace
                strng = line[i+1:x].strip(' ,')
                if len(strng) == 0:
                    continue

                # remove apostrophes from ends if they exist
                if (strng[0] == "\'" and strng[-1] == "\'") or (strng[0] == "\"" and strng[-1] == "\""):
                    strng = strng[1:-1]
                    strng = strng.strip()

                # change fortran boolean to python boolean
                if not in_custom:
                    if strng == "T" or strng.upper() == ".TRUE.": strng = True
                    elif strng == "F" or strng.upper() == ".FALSE.": strng = False
                # Check if srtng is a number
                try:
                    float(strng)
                    isFl = True
                except: isFl = False

                if key[:5]=='MODS(':
                    y = key.find(')')
                    index = int(key[5:y])
                    if index>len(NAMES): continue
                    if index-1 == divider_i or index-1 == divider_xtr_i:
                        continue
                    name = NAMES[index-1]
                    if not name in vars.mods:
                        vars.mods[name] = Comp()
                        vars.mods[name].Find = namesFoInds[name]
                        vars.mods[name].name = name

                    if len(key)>y+1:
                        prop = key[y+2:].lower()
                        cmd = 'vars.mods[\'%s\'].%s'%(NAMES[index-1],prop)
                        if isFl:
                            exec("%s = %s"%(cmd, strng))
                        else:
                            exec("%s = \'%s\'"%(cmd, strng))
                    else:
                        props = strng.replace('d','e').split()
                        i = 0
                        n=len(props)
                        if i< n: vars.mods[name].mode   = int(props[i])   ; i=i+1
                        if i< n: vars.mods[name].col    = int(props[i])   ; i=i+1
                        if i< n: vars.mods[name].multi  = float(props[i]) ; i=i+1
                        if i< n: vars.mods[name].shift  = float(props[i]) ; i=i+1
                        if i< n: vars.mods[name].min    = float(props[i]) ; i=i+1
                        if i< n: vars.mods[name].max    = float(props[i]) ; i=i+1
                        if i< n: vars.mods[name].sig    = float(props[i]) ; i=i+1
                        if i< n: vars.mods[name].mju    = float(props[i]) ; i=i+1
                        if i< n: vars.mods[name].fv     = float(props[i]) ; i=i+1
                        if i< n: vars.mods[name].ph     = float(props[i]) ; i=i+1
                        if i< n: vars.mods[name].am     = float(props[i]) ; i=i+1
                        if i< n:
                            unt = props[i].strip('\'\"').replace('#','#/cm3')
                            vars.mods[name].unit = unt
                            i=i+1
                        if i< n:
                            vars.mods[name].tied = (props[i].strip('\'\"'))
                            if vars.mods[name].tied in NAMES:
                                vars.mods[name].col  = vars.mods[name].tied

            else:
                if line.strip() == '/' and in_custom:
                    in_custom = False
                    n = len(nml.CUSTOM.CUSTOMS)
                    ii = 0
                    for i in range(1,n+1):
                        ii += 1
                        keyW = 'customKey_%d'%ii
                        valueW = 'customVal_%d'%ii
                        exec("self.%s.setText(\'%s\')"%(keyW,nml.CUSTOM.CUSTOMS[ii-1][0]))
                        exec("self.%s.setText(\'%s\')"%(valueW,nml.CUSTOM.CUSTOMS[ii-1][1]))
                    for j in range(ii+1,31):
                        keyW = 'customKey_%d'%j
                        valueW = 'customVal_%d'%j
                        exec("self.%s.clear()"%(keyW))
                        exec("self.%s.clear()"%(valueW))

                elif '&NML_CUSTOM' in line:
                    in_custom = True
                    nml.CUSTOM.CUSTOMS = []
                continue

            if 'INOUT_DIR' == key: self.inout_dir.setText(strng)
            elif 'CASE_NAME' == key: self.case_name.setText(strng)
            elif 'RUN_NAME' == key: self.run_name.setText(strng)
            elif 'CHEMISTRY_FLAG' == key: self.checkBox_che.setChecked(strng)
            elif 'AEROSOL_FLAG' == key: self.checkBox_aer.setChecked(strng)
            elif 'ACDC_SOLVE_SS' == key: self.acdc_solve_ss.setChecked(strng)
            elif 'ACDC' == key: self.checkBox_acd.setChecked(strng)
            elif 'CONDENSATION' == key: self.condensation.setChecked(strng)
            elif 'COAGULATION' == key: self.coagulation.setChecked(strng)
            elif 'DEPOSITION' == key: self.deposition.setChecked(strng)
            elif 'CHEM_DEPOSITION' == key: self.chemDeposition.setChecked(strng)
            # elif 'MODEL_H2SO4' == key: self.model_h2so4.setChecked(strng)
            elif 'ORG_NUCL' == key: self.Org_nucl.setChecked(strng)
            elif 'RUNTIME' == key and isFl: self.runtime.setValue(float(strng))
            elif 'DT' == key and isFl: self.dt.setValue(float(strng))
            elif 'FILE_TIME_UNIT' == key:
                iii = timeunits.get(strng, 0)
                self.fileTimeUnit_a.setCurrentIndex( iii )
                self.fileTimeUnit_b.setCurrentIndex( iii )
            elif 'LOSSFILE_TIME_UNIT' == key:
                iii = timeunits.get(strng, 0)
                self.lossfileTimeUnit.setCurrentIndex( iii )
            elif 'PRINT_ACDC' == key: self.print_acdc.setChecked(strng)
            elif 'OPTIMIZE_DT' == key: self.useSpeed.setChecked(strng)
            elif 'AFTER_CHEM_ON' == key: self.after_chem_on.setChecked(strng)
            elif 'AFTER_NUCL_ON' == key: self.after_nucl_on.setChecked(strng)
            elif 'FSAVE_INTERVAL' == key and isFl: self.fsave_interval.setValue(int(float(strng)))
            elif 'PRINT_INTERVAL' == key and isFl: self.print_interval.setValue(int(float(strng)))
            elif 'FSAVE_DIVISION' == key and isFl: self.fsave_division.setValue(int(float(strng)))
            elif 'DATE' == key: self.dateEdit.setDate(parse_date(strng))
            elif 'INDEX' == key and isFl: self.indexEdit.setValue(int(float(strng)))
            elif 'NUMBER' == key and isFl: self.indexEdit.setValue(int(float(strng)))
            elif 'PSD_MODE' == key and isFl: self.psd_mode.setCurrentIndex(int(float(strng)))
            elif 'N_BINS_PAR' == key and isFl: self.n_bins_particle.setValue(int(float(strng)))
            elif 'MIN_PARTICLE_DIAM' == key: self.min_particle_diam.setText(strng)#   1.0000000000000001E-009,
            elif 'MAX_PARTICLE_DIAM' == key: self.max_particle_diam.setText(strng)#   9.9999999999999995E-007,
            elif 'N_MODAL' == key: self.n_modal.setText(strng)#   9.9999999999999995E-007,
            elif 'MMODAL_INPUT_INUSE' == key: self.multiModalBox.setChecked(int(float(strng)))#   9.9999999999999995E-007,
            elif 'DMPS_FILE' == key: self.dmps_file.setText(strng)
            elif 'EXTRA_PARTICLES' == key: self.extra_particles.setText(strng)
            elif 'MMODAL_INPUT' == key: self.mmodal_input.setText(strng)
            elif 'DMPS_READ_IN_TIME' == key and isFl: self.dmps_read_in_time.setValue(float(strng))
            elif 'DMPS_HIGHBAND_LOWER_LIMIT' == key: self.dmps_highband_lower_limit.setText(strng)
            elif 'DMPS_LOWBAND_UPPER_LIMIT' == key: self.dmps_lowband_upper_limit.setText(strng)
            elif 'USE_DMPS' == key: self.use_dmps.setChecked(strng)
            elif 'USE_DMPS_PARTIAL' == key: self.use_dmps_partial.setChecked(strng)
            elif 'DMPS_INTERVAL' == key: self.dmpsIntvl.setValue(float(strng))
            elif 'ENV_FILE' == key: self.env_file.setText(strng)
            elif 'SPECTRUMFILE' == key: self.spectralFunctions.setText(strng)
            elif 'SWR_IS_ACTINICFLUX' == key: self.SW_is_AF.setChecked(strng)
            elif 'CHAMBER_FLOOR_AREA' == key and isFl: self.floorArea.setValue(float(strng))#  0.20000000000000001     ,
            # elif 'CHAMBER_CIRCUMFENCE' == key and isFl: self.chamberCircumfence.setValue(float(strng))#  0.20000000000000001     ,
            elif 'CHAMBER_HEIGHT' == key and isFl: self.chamberHeight.setValue(float(strng))#  0.20000000000000001     ,
            elif 'EDDYK' == key and isFl: self.eddyK.setText(strng)#  0.05000000000000001     ,
            elif 'USTAR' == key and isFl: self.ustar.setText(strng)#  0.050000000000000001     ,
            elif 'ALPHAWALL' == key and isFl: self.alphaWall.setText(strng)#  0.050000000000000001     ,
            elif 'CW_EQV' == key and isFl: self.Cw_eqv.setText(strng)#  0.050000000000000001     ,
            elif 'MCM_FILE' == key: self.mcm_file.setText(strng)# "
            elif 'LOSSES_FILE' == key: self.losses_file.setText(strng)# "
            elif 'LAT' == key and isFl: self.lat.setValue(float(strng))
            elif 'LON' == key and isFl: self.lon.setValue(float(strng))
            elif 'WAIT_FOR' == key and isFl: self.wait_for = (int(float(strng)))
            elif 'DESCRIPTION' == key: self.description.setPlainText(strng.replace('<br>','\n'))# "Just some keying
            elif 'CH_ALBEDO' == key and isFl: self.ch_albedo.setValue(float(strng))#  0.20000000000000001     ,
            elif 'SWR_IN_LOWER' == key and isFl: self.swr_in_lower.setValue(int(float(strng)))#  0.20000000000000001     ,
            elif 'SWR_IN_UPPER' == key and isFl: self.swr_in_upper.setValue(int(float(strng)))#  0.20000000000000001     ,
            elif 'USE_ATOMS' == key: self.use_atoms.setChecked(strng)
            elif 'VAP_NAMES' == key: self.vap_names.setText(strng)
            elif 'VAP_ATOMS' == key: self.vap_atoms.setText(strng)
            elif 'GR_SIZES' == key: self.GR_sizes.setText(strng)
            elif 'NAMESDAT' == key:
                if strng != path_to_names:
                    self.popup('Hazard', """These settings need different set of input compounds, you should restart the GUI \n
                    and use -o flag in the bash command, or --names flag in the python command to define the correct path \n"""
                    +strng, 3)

            elif 'DDIAM_RANGE' == key:
                self.prec_low_1.setValue(float(strng.split(',')[0]))
                self.prec_high_1.setValue(float(strng.split(',')[1]))
            elif 'DPNUM_RANGE' == key:
                self.prec_low_2.setValue(float(strng.split(',')[0]))
                self.prec_high_2.setValue(float(strng.split(',')[1]))
            elif 'DVAPO_RANGE' == key:
                self.prec_low_3.setValue(float(strng.split(',')[0]))
                self.prec_high_3.setValue(float(strng.split(',')[1]))

            elif 'ACDC_SYSTEMS' == key:
                if ',' in strng:
                    self.acdc_systems_flags = [int(z) for z in strng.split(',')]
                elif '*' in strng:
                    self.acdc_systems_flags = [int(strng[2:])]*int(strng[0])
                lll = len(self.acdc_systems_flags)
                if lll==0: self.acdc_systems_flags = [1,1,0,0,0]
                if lll<self.ACDC_n_systems: self.acdc_systems_flags = self.acdc_systems_flags + [0]*(self.ACDC_n_systems-lll)

            elif len(re.findall('ACDC_LINKS\((\d)\)',key))>0:
                num = int(re.findall('ACDC_LINKS\((\d)\)',key)[0])
                if num<=self.ACDC_n_systems:
                    self.ACDC_current_links[num-1][1] = {a:b for a,b in zip(strng.split()[0::2],strng.split()[1::2])}
                    if self.ACDC_available_compounds[num-1][1].keys() != self.ACDC_current_links[num-1][1].keys():
                        INC_COMP = True
                        msg = """The existing ACDC system # %d does not contain the same compounds that were defined in the setting file.:
In the ACDC system: %s
In the loaded settings: %s""" %(num, ' '.join(self.ACDC_available_compounds[num-1][1].keys()),', '.join(self.ACDC_current_links[num-1][1].keys()))
                        self.popup('Incompatible definitions', msg)


            elif in_custom:
                if strng == 'T': strng = '.TRUE.'
                if strng == 'F': strng = '.FALSE.'
                if key=='DMPS_TRES_MIN':
                    print('Moving obsoleted DMPS_TRES_MIN to dmpsIntvl')
                    self.dmpsIntvl.setValue(float(strng))
                else:
                    nml.CUSTOM.CUSTOMS.append([key, strng])

            elif '# RAW_INPUT' == key:
                self.rawEdit.clear()
                self.rawEdit.insertPlainText(rawline.replace('<br>', '\n'))
            elif '# INPUT_SETTINGS' == key:
                sets = strng.split()
                for kv in sets:
                    sp_i = kv.find(':')
                    kk, val = kv[:sp_i], kv[sp_i+1:]
                    if 'FALSE' in val or 'TRUE' in val:
                        if 'TRUE' in val:
                            exec('self.'+kk+'.setChecked(True)')
                        else:
                            exec('self.'+kk+'.setChecked(False)')
                    else:
                        exec('self.'+kk+'.setText(\''+val+'\')')

            elif '# BATCH_SETTINGS' == key:
                sets = strng.split()
                for kv in sets:
                    sp_i = kv.find(':')
                    kk, val = kv[:sp_i], kv[sp_i+1:]
                    if 'FALSE' in val or 'TRUE' in val:
                        if 'TRUE' in val:
                            exec('self.'+kk+'.setChecked(True)')
                        else:
                            exec('self.'+kk+'.setChecked(False)')
                    elif 'batchRangeInd' in kk:
                        exec('self.'+kk+'.setValue('+val+')')
                    elif 'batchRangeDay' in kk:
                        exec('self.'+kk+'.setDate(parse_date(\''+val+'\'))')
            else:
                continue

        f.close()
        # manage different units
        for key in vars.mods:
            pmInUse = 0
            if vars.mods[key].mode > 0:
                pmInUse = 1
            unts = units.get(grepunit(key),units['REST'])
            unitIndex = 0
            for unit in unts:
                if unit.upper() == vars.mods[key].unit.upper(): break
                else: unitIndex += 1
            if unitIndex==len(unts): unitIndex=0

            if key == 'TEMPK' or key == 'PRESSURE':
                if key == 'TEMPK': row = 0
                if key == 'PRESSURE': row = 1
                self.selected_vars.setItem(row, 1, QtWidgets.QTableWidgetItem(str(vars.mods[key].col)))
                self.selected_vars.setItem(row, 2, QtWidgets.QTableWidgetItem(str(vars.mods[key].multi)))
                self.selected_vars.setItem(row, 3, QtWidgets.QTableWidgetItem(str(vars.mods[key].shift)))
                self.selected_vars.cellWidget(row,4).setCurrentIndex(pmInUse)
                self.selected_vars.cellWidget(row,5).setCurrentIndex(unitIndex)
                for z in range(4):
                    if namesPyInds[key]<divider_i:
                        self.selected_vars.item(row, z).setBackground(QtGui.QColor(*env_no))
                    elif namesPyInds[key]>divider_xtr_i:
                        self.selected_vars.item(row, z).setBackground(QtGui.QColor(*xtr_no))
                    else:
                        self.selected_vars.item(row, z).setBackground(QtGui.QColor(*org_no))

            else:
                cols = [key,str(vars.mods[key].col),str(vars.mods[key].multi),str(vars.mods[key].shift),pmInUse]
                self.add_new_line(key,2,cols=cols,createNew=False, unt=unitIndex)

            wScale,pScale,aScale,phScale,ampScale,_ = self.scales()
            vars.mods[key].sliderVls[0] = int(round(vars.mods[key].sig /wScale))
            vars.mods[key].sliderVls[1] = int(round(vars.mods[key].mju /pScale))
            vars.mods[key].sliderVls[2] = int(round(vars.mods[key].fv  /aScale))
            vars.mods[key].sliderVls[3] = int(round(vars.mods[key].ph  /phScale))
            vars.mods[key].sliderVls[4] = int(round(vars.mods[key].am  /ampScale))
            for i in range(5):
                trgt = vars.mods[key].sliderVls[i]
                for j in range(0,11):
                    if max(j, 0.5)*slMxs[i] > trgt:
                        vars.mods[key].sl_x[i] = j
                        break

        self.fileLoadOngoing = False
        if INC_COMP:
            INC_COMP = False
            self.popup('ACDC incompatibility','Check the linking of ACDC compounds.')
        self.set_ACDC_links()
        self.updatePath()
        self.updateEnvPath()



    def updateMass(self):
        if self.MPD != []:
            self.showMass(first=False, target='mass')


    def updateNumbers(self):
        if self.MPD != []:
            self.showMass(first=False, target='numb')

    def testNC(self,file):
        try:
            self.ncs = netCDF4.Dataset(file, 'r')
            self.ncs.close()
            return 0
        except:
            return 1

    def showMass(self, file=None, first=True, target=None, add=False):
        symbols = ['x','+','t1','t','star']

        if float('.'.join(pg.__version__.split('.')[0:2]))<0.13:
            self.unnecessaryLegendTrick()

        if first:
            if file==None: return
            if add:
                self.plotResultWindow_2.clear()
                self.plotResultWindow_3.clear()
                self.msspltsLines = []
            else:
                self.closenetcdf_mass()
            self.loadNetcdf_massAdd.setEnabled(True)
            if exists(file):
                if self.testNC(file)>0:
                    self.popup('Bummer...', 'Not a valid output file',icon=3)
                    if len(self.MPD)==0:
                        return
                    else:
                        self.showMass(first=False)
                        return
            if add:
                if len(self.MPD)==4:
                    self.popup('No more files', 'You can have maximum 4 files open in this tool')
                else:
                    self.MPD.append(NcPlot(file))
                    if len(self.MPD[0].diameter) != len(self.MPD[-1].diameter):
                        self.MPD.pop(len(self.MPD)-1)
                        self.popup('Cannot load new file','PSD must match with the first file.')
                    elif npsum(npround(log(self.MPD[0].diameter),2) - npround(log(self.MPD[-1].diameter), 2))>0:
                        self.MPD.pop(len(self.MPD)-1)
                        self.popup('Cannot load new file','PSD must match with the first file.')
                    elif int(self.MPD[0].timeint*10) != int(self.MPD[-1].timeint*10):
                        self.MPD.pop(len(self.MPD)-1)
                        self.popup('Cannot load new file','Time interval must match with the first file.')
                    elif self.MPD[0].time[0] != self.MPD[-1].time[0]:
                        self.MPD.pop(len(self.MPD)-1)
                        self.popup('Cannot load new file','Start time must match with the first file.')
                    if round(self.MPD[-1].time[-1],3)>round(self.mp_time[-1],3):
                        self.mp_time = self.MPD[-1].time
                        self.times.itemSelectionChanged.disconnect(self.updateNumbers)
                        self.times.clear()
                        self.times.itemSelectionChanged.connect(self.updateNumbers)
                        self.times.addItems(['%7.2f'%(i) for i in self.mp_time])
            else:
                self.MPD = [NcPlot(file)]
            if not add:
                self.diams.addItems(['%7.2f'%(1e9*i) for i in self.MPD[0].diameter])
                self.diams.setSelectionMode(QtWidgets.QAbstractItemView.ExtendedSelection)
                self.mp_time = self.MPD[0].time
                indstime = [0]
                self.plotResultWindow_3.setLogMode(x=True)
                self.plotResultWindow_2.setLabel('bottom', 'Time', units='h')
                self.plotResultWindow_2.setLabel('left', 'Mass', units='g')
                self.plotResultWindow_3.setLabel('bottom', 'Diameter', units='m')
                self.plotResultWindow_3.setLabel('left', '# normalized')
            # NUMBER_CONCENTRATION = self.ncs_mass.variables['NUMBER_CONCENTRATION'][:]
            # self.mass_in_bin     = self.ncs_mass.variables['MASS'][:]*NUMBER_CONCENTRATION*1e3 # from kg to g
            # self.lognormconc     = NUMBER_CONCENTRATION/log10(self.DIAMETER[0,1]/self.DIAMETER[0,0])
                self.times.addItems(['%7.2f'%(i) for i in self.mp_time])
                self.times.setSelectionMode(QtWidgets.QAbstractItemView.ExtendedSelection)
                # self.massPlotTitle.setText(self.MPD[0].legend.replace(': Particles.nc',''))
                self.times.item(0).setSelected(True)
                self.diams.selectAll()
                self.diams.itemSelectionChanged.connect(self.updateMass)
                self.times.itemSelectionChanged.connect(self.updateNumbers)
            # try:
            #     DMPS_CONCENTRATION = self.ncs_mass.variables['INPUT_CONCENTRATION'][:]
            #     self.massdmps = self.ncs_mass.variables['MASS'][:]*DMPS_CONCENTRATION
            #     self.lognormdmps = DMPS_CONCENTRATION/log10(self.DIAMETER[0,1]/self.DIAMETER[0,0])
            #     self.measdmps = True
            # except:
            #     print('File does not contain measured PSD')

            # self.ncs_mass.close()

        # self.NC_lines[j].setData(TT[j],Y)
        # Update figures, start with cleaning
        # else:
        #     if target=='mass':
        self.plotResultWindow_2.clear()
            # if target=='numb':
        self.plotResultWindow_3.clear()

        indstime = [i.row() for i in self.times.selectedIndexes()]
        inds = [i.row() for i in self.diams.selectedIndexes()]
        # y is the array that gets plotted
        miny, maxy = 0,0
        for j,m in enumerate(self.MPD):
            y  = npsum(self.MPD[j].totalmass_ng_m[:,inds],axis=1)*1e-9
            y[y<1e-22] = 0e0
            miny, maxy = min(miny,y.min()),max(maxy,y.max())

        if self.MPD[0].measdmps:
            y2 = npsum(self.MPD[0].massdmps[:,inds],axis=1)*1e-9
            miny, maxy = min(miny,y2.min()),max(maxy,self.showAlsoMeasInMassConc.isChecked()*y2.max())
        if maxy>0:
            if abs(1-miny/maxy) <1e-12:
                maxy = miny*100000
        if self.MPD[0].measdmps and self.showAlsoMeasInMassConc.isChecked():
            self.outplot_mass = self.plotResultWindow_2.plot(self.mp_time,
                                y2,
                                pen={'color':'r','width': 2.0,'style': QtCore.Qt.DotLine},
                                symbol='x',
                                symbolPen='b',
                                symbolBrush='b',
                                symbolSize=6,
                                name='Model init'
                                )

        for j,mpd_file in enumerate(self.MPD):
            self.outplot_mass = self.plotResultWindow_2.plot(mpd_file.time,
                                npsum(mpd_file.totalmass_ng_m[:,inds],axis=1)*1e-9,
                                pen={'color':colors[j],'width': 2.0},
                                symbol=symbols[j+1],
                                symbolPen=colors[j],
                                symbolBrush=colors[j],
                                symbolSize=8,
                                name=self.MPD[j].legend.replace(': Particles.nc','')
                                )

        if target == 'mass' or first: self.plotResultWindow_2.setRange(yRange=[miny*0.95,maxy*1.05])
        nmax = 0.01
        for i,ii in enumerate(indstime):
            if self.MPD[0].measdmps and self.showAlsoMeasInMassConc.isChecked():
                nmax = max(nmax,max(self.MPD[0].lognormdmps[ii,:]))
                self.outplot_numb2 = self.plotResultWindow_3.plot(self.MPD[0].diameter,
                                self.MPD[0].lognormdmps[ii,:],
                                pen={'color':colors[i%10],'width': 2.0,'style': QtCore.Qt.DotLine},
                                symbol='x',
                                symbolPen=colors[i%10],
                                symbolBrush=colors[i%10],
                                symbolSize=6,
                                name='Model init'
                                )
            for j,mpd_file in enumerate(self.MPD):
                if ii >= len(mpd_file.time): continue
                nmax = max(nmax,max(mpd_file.lognorm_nc_cm3[ii,:]))
                self.outplot_numb = self.plotResultWindow_3.plot(self.MPD[0].diameter,
                                mpd_file.lognorm_nc_cm3[ii,:],
                                pen={'color':colors[i%10],'width': 2.0},
                                symbol=symbols[j+1],
                                symbolPen=colors[i%10],
                                symbolBrush=colors[i%10],
                                symbolSize=8,
                                )

        # yr=self.plotResultWindow_3.getViewBox().state['viewRange'][1]
        self.plotResultWindow_3.setRange(yRange=[-0.02*nmax,nmax*1.05])

#
    def toggleppm(self,what):
        if what == 'off':
            self.n_conc.setChecked(True)
            self.ppm.setEnabled(False)
            self.ppb.setEnabled(False)
            self.ppt.setEnabled(False)

        if what == 'on':
            self.ppm.setEnabled(True)
            self.ppb.setEnabled(True)
            self.ppt.setEnabled(True)

    def addAnotherNC(self):
        if self.LPD == []:
            return
        else:
            ftype = self.LPD[-1].masterfile
        self.browse_path(None, 'addplot_more', ftype=ftype)



    def linePlotMulti(self, file, new=True):
        self.findComp.clear()
        if new or self.fLin_2.isChecked():
            putBackLog = False
        if not new and self.fLog_2.isChecked():
            putBackLog = True
            self.plotResultWindow.setLogMode(y=False)
        self.fLin_2.setChecked(True)
        self.fLog_2.setChecked(False)
        if new:
            self.plotResultWindow.setLogMode(y=False)
            self.findComp.clear()
            # Close all previus netcdf-files and clear plot
            self.closenetcdf()
            self.addSimilar.setEnabled(True)
            self.actionShow_variable_attributes.setEnabled(True)
            self.CloseLinePlotsButton.setEnabled(True)
            self.findComp.setEnabled(True)

        else:
            comp = self.availableVars.currentItem().text()
        # Try to open netCDF-file
        if exists(file):
            if self.testNC(file)>0:
                self.popup('Bummer...', 'Not a valid output file',icon=3)
                return
        else: return
        # the file checks, move on
        if new:
            if 'General.nc' in file:
                mdim = True
                self.LPD = [NcPlot(file,mdim=mdim)]
            else:
                mdim = False
                self.LPD = [NcPlot(file)]
            if self.LPD[-1].masterfile == 'Particles.nc': self.ShowPPC.setEnabled(True)
        else:
            if 'General.nc' in file:
                mdim = True
                self.LPD.append(NcPlot(file,mdim=mdim))
            else:
                mdim = False
                self.LPD.append(NcPlot(file,mdim=mdim))
            if not self.LPD[0].convars.keys() <= self.LPD[-1].convars.keys():
                self.LPD.pop(len(self.LPD)-1)
                self.popup('No can do', 'The file you are trying to add has less compounds than the original, cannot add the new file. Loading the file which has less compounds first might work.')
                return
        if new:
            Y,unit = self.LPD[-1].getloc(0, return_unit=True)
            if unit == '[#/cm3]': self.toggleppm('on')
        else:
            Y = self.LPD[-1].getconc(comp)

        self.NC_lines.append( self.plotResultWindow.plot(
        self.LPD[-1].time,
        Y,
        pen={'color':colors[(len(self.LPD)-1)],'width': 2.0},
        name=self.LPD[-1].legend
        ) )
        if new:
            self.availableVars.clear()
            self.availableVars.addItems(self.LPD[-1].names)
            self.availableVars.item(0).setSelected(True)
            self.availableVars.setCurrentItem(self.availableVars.item(0))
            self.availableVars.itemSelectionChanged.connect(self.showOutputUpdate)
            self.ShowPPC.toggled.connect(lambda: self.showOutputUpdate(info=False))
        # If fails, give information and return
        # except:
        #     self.popup('Bummer...', "Output file does not contain any plottable data",icon=3)
        #     return
        if new:
            # All's well, finish plot
            self.plotResultWindow.setLabel('bottom', 'Time', units='h')
            self.plotTitle = file + ': ' + list(self.LPD[-1].convars.keys())[0]
            self.plotResultTitle.setText(self.plotTitle+' '+unit)
        else:
            self.fLin_2.setChecked(~putBackLog)
            self.fLog_2.setChecked(putBackLog)
            self.showOutputUpdate()


    def selectionMode(self):
        if self.sumSelection.isChecked():
            self.availableVars.setSelectionMode(QtWidgets.QAbstractItemView.ExtendedSelection)
        else:
            self.availableVars.setSelectionMode(QtWidgets.QAbstractItemView.SingleSelection)
            if self.availableVars.currentItem() != None:
                self.availableVars.currentItem().setSelected(True)


    def loadsecondFile(self):
        newfile = self.browse_path(target=None,mode='fixed_cc', ftype=self.plotMasterFile)
        self.listOfplottedFiles.append(newfile)
        print(self.listOfplottedFiles)


    def showOutputUpdate(self,info=False):
        """This function is called when lin/log radio button or any variable in the list is changed"""
        # find out which y-scale should be used
        if self.ShowPPC.isChecked():
            self.toggleppm('off')
            PPconc = True
            self.plotResultWindow.setBackground('#F3F0EF')

        else:
            PPconc = False
            self.plotResultWindow.setBackground('w')
        scale = self.radio(self.fLin_2, self.fLog_2, action=False)
        if scale == 'log':loga = True
        else: loga = False

        # find out which variable should be plotted
        if self.availableVars.currentItem() is None:
            return
        else:
            comp = self.availableVars.currentItem().text()
            if self.LPD[-1].masterfile == 'Particles.nc':
                csat = f' c*:{self.LPD[-1].csat[comp]:7.1e} {self.LPD[-1].csat_unit}'
            else:
                csat = ''
        if not PPconc:
            if self.sumSelection.isChecked():
                if any(c.text() in units for c in self.availableVars.selectedItems()):
                    self.n_conc.setChecked(True)
                    self.toggleppm('off')
                else:
                    self.toggleppm('on')
            else:
                if comp not in units:
                    self.toggleppm('on')
                else:
                    self.n_conc.setChecked(True)
                    self.toggleppm('off')
        if info:
            self.LPD[-1].getinfo(comp)
            return
        # Exctract that variable from netCDF-dataset and save to Y
        YY = []
        TT = []
        if not self.sumSelection.isChecked():
            if PPconc:
                for z in self.LPD:
                    YY.append(z.getcom(comp))
                    TT.append(z.time)
            else:
                for z in self.LPD:
                    YY.append(z.getconc(comp))
                    TT.append(z.time)

            self.plotTitle = self.plotTitle[:self.plotTitle.rfind(':')+2] + comp+' '
        else:
            if self.availableVars.selectedItems() != []:
                if PPconc:
                    for z in self.LPD:
                        YY.append(z.getcomsum([c.text() for c in self.availableVars.selectedItems()]))
                        TT.append(z.time)
                else:
                    for z in self.LPD:
                        YY.append(z.getconcsum([c.text() for c in self.availableVars.selectedItems()]))
                        TT.append(z.time)
                self.plotTitle = self.plotTitle[:self.plotTitle.rfind(':')+2] + self.availableVars.selectedItems()[0].text()+' etc. '
            else:
                if PPconc:
                    for z in self.LPD:
                        YY.append(z.getcom(comp))
                        TT.append(z.time)
                else:
                    for z in self.LPD:
                        YY.append(z.getconc(comp))
                        TT.append(z.time)
                self.plotTitle = self.plotTitle[:self.plotTitle.rfind(':')+2] + comp+' '

        for j,Y in enumerate(YY):
            if max(Y)!=0:
                # if the variable is PRACTICALLY nonvariant, set the value to constant
                if (max(Y)-min(Y))/max(Y)<1e-5: YY[j][:]=mean(Y)

        have_aircc = False
        if all([z.have_aircc for z in self.LPD]): have_aircc = True

        un = ''
        if have_aircc:
            for j,Y in enumerate(YY):
                if self.ppm.isChecked():
                    YY[j] = Y/self.LPD[j].aircc * 1e6
                    un = '[ppm]'
                elif self.ppb.isChecked():
                    YY[j] = Y/self.LPD[j].aircc * 1e9
                    un = '[ppb]'
                elif self.ppt.isChecked():
                    YY[j] = Y/self.LPD[j].aircc * 1e12
                    un = '[ppt]'
                elif PPconc:
                    un = '[ng/m3]'
                elif self.sumSelection.isChecked():
                    if any(c.text() in units for c in self.availableVars.selectedItems()):
                        if len(self.availableVars.selectedItems())==1:
                            un = '['+units.get(grepunit(self.availableVars.selectedItems()[0].text()),units['REST'])[0]+']'
                        else:
                            un = '[-]'
                    else:
                        un = '[#/cm3]'
                else : un = '['+units.get(grepunit(comp),units['REST'])[0]+']'


        # Are the values non-negative?
        positive = all([all(Y>=0) for Y in YY])
        # Are all the values zeros?
        zeros = all([all(Y==0) for Y in YY])
        thickness = 2.0
        if PPconc: thickness = 2.5
        # if non-negative and not all zeros, and log-scale is possible without any fixes
        if not zeros and not positive and scale=='log':
            for j,Y in enumerate(YY):
                YY[j] = where(Y<=0, 1e-20, Y)
            # Mark plot with red to warn that negative values have been deleted
                self.NC_lines[j].setPen({'style':QtCore.Qt.DashLine,'color':colors[j],'width': thickness})
            positive = True
        else:
            for j,p in enumerate(self.NC_lines):
                p.setPen({'style':QtCore.Qt.SolidLine,'color':colors[j],'width': thickness})

        if scale=='log' and positive and not zeros:
            # To avoid very small exponentials, zeros are changed to nearest small number
            if not all([all(Y>0) for Y in YY]):
                for j,Y in enumerate(YY):
                    uniqs = unique(Y)
                    if len(uniqs)>1:
                        YY[j] = where(Y==0,uniqs[1],Y)
                    # if the above fails, use linear scale
                    else:
                        loga = False
            # if everything ok, use log scale
            loga = True

        # if linear scale was chosen, or negative values, or all zeros, use lin-scale
        else:
            self.plotResultWindow.setLogMode(y=False)
            loga = False
        # update data
        for j,Y in enumerate(YY):
        # time = self.ncs.variables[self.hnames[0]][:]/3600
        # if ma.is_masked(time):
            self.NC_lines[j].setData(TT[j],Y)
        # else:
        #     self.outplot.setData(time,Y)
        self.plotResultWindow.setLogMode(y=loga)
        # update title
        self.plotResultTitle.setText(self.plotTitle+un+csat)



    ## Popup message function -icon sets the icon:----------------------------------------------------------------------
    # QMessageBox::NoIcon      0   the message box does not have any icon.
    #              Information 1   an icon indicating that the message is nothing out of the ordinary.
    #              Warning     2   an icon indicating that the message is a warning, but can be dealt with.
    #              Critical    3   an icon indicating that the message represents a critical problem.
    #              Question	   4   an icon indicating that the message is asking a question.
    def popup(self,title,message,icon=2):
        """Handle for giving popup messages"""
        msg = QtWidgets.QMessageBox()
        msg.setIcon(icon)
        msg.setWindowTitle(title)
        msg.setText(message)
        retval = msg.exec_()


    def closenetcdf(self):
        # try: self.ncs.close()
        # except: pass
        self.findComp.clear()
        for z in self.LPD:
            z.closenc()
        self.LPD = []
        self.NC_lines = []
        self.addSimilar.setEnabled(False)
        self.actionShow_variable_attributes.setEnabled(False)
        self.CloseLinePlotsButton.setEnabled(False)
        self.findComp.setEnabled(False)

        if float('.'.join(pg.__version__.split('.')[0:2]))<0.13:
            if self.ncleg in self.ncleg_skene.items():
                self.ncleg_skene.removeItem(self.ncleg)
                self.ncleg = self.plotResultWindow.addLegend()
                if self.ncleg not in self.ncleg_skene.items():
                    self.ncleg_skene.addItem(self.ncleg)


        try:
            while True:
                try: self.availableVars.itemSelectionChanged.disconnect(self.showOutputUpdate)
                except TypeError: break
            self.availableVars.clear()
            self.plotResultWindow.clear()
            self.ShowPPC.toggled.disconnect(self.showOutputUpdate)
            # self.ShowPPC.toggled.disconnect(self.toggleTimeInavailableVars)

        except:
            pass
        self.ShowPPC.setChecked(False)
        self.ShowPPC.setEnabled(False)
        self.toggleppm('off')
        self.listOfplottedFiles = []

    def unnecessaryLegendTrick(self):
        if self.massleg in self.massleg_skene.items():
            self.massleg_skene.removeItem(self.massleg)
            self.massleg = self.plotResultWindow_2.addLegend()
            if self.massleg not in self.massleg_skene.items():
                self.massleg_skene.addItem(self.massleg)

    def closenetcdf_mass(self):
        self.MPD = []
        self.loadNetcdf_massAdd.setEnabled(False)
        if float('.'.join(pg.__version__.split('.')[0:2]))<0.13:
            self.unnecessaryLegendTrick()
        try: self.ncs_mass.close()
        except: pass
        try:
            while True:
                try:
                    self.diams.itemSelectionChanged.disconnect(self.updateMass)
                    self.times.itemSelectionChanged.disconnect(self.updateNumbers)
                except TypeError: break
            self.diams.clear()
            self.times.clear()
            self.plotResultWindow_2.clear()
            self.plotResultWindow_3.clear()
        except:
            pass


    def progress(self):
        i = self.compileProgressBar.value()
        if i==100 or i==0:
            if self.inv == -1:
                self.compileProgressBar.setInvertedAppearance(True)
                self.inv = 1
            else:
                self.compileProgressBar.setInvertedAppearance(False)
                self.inv = -1
        self.compileProgressBar.setValue(i-self.inv)
        self.running = self.compile.poll()

        if self.running != None:
            self.TimerCompile.stop()
            self.recompile.setEnabled(True)
            self.compileProgressBar.hide()
            self.compileProgressBar.setInvertedAppearance(False)
            self.compileProgressBar.setValue(0)
            if self.running == 0:
                self.popup('Compiled', 'Compiling done', icon=0)
            else:
                self.popup('Compile error', 'Error in compiling, see output from terminal', icon=2)
            if self.ReplChem.isChecked(): self.ReplChem.setChecked(False)
            if self.makeClean.isChecked(): self.makeClean.setChecked(False)


    def remake(self, syst=0):
        if self.running != None:
            if self.ReplChem.isChecked():
                self.editMakefile(mod=self.chemistryModules.currentText())
                writeRATESdat(self.chemistryModules.currentText())
                if self.makeClean.isChecked():
                    target = "clean_current_chemistry"
                else:
                    target = "clean"
                self.compile = Popen(["make", target])#, stdout=subprocess.PIPE,stderr=subprocess.STDOUT,stdin=None)
                while True:
                    self.running = self.compile.poll()
                    if self.running != None: break
            if syst!=0:
                self.compile = Popen(["make", "cleanoneacdc","ACDCTARG=%s"%syst])#, stdout=subprocess.PIPE,stderr=subprocess.STDOUT,stdin=None)
            self.compile = Popen(["make"])#, stdout=subprocess.PIPE,stderr=subprocess.STDOUT,stdin=None)
            self.recompile.setEnabled(False)
            self.compileProgressBar.show()
            self.inv = 1
            self.compileProgressBar.setValue(0)
            self.TimerCompile.start(10)
            self.chemLabel.setText('Current chemistry scheme in makefile: '+self.get_available_chemistry(checkonly=True))
            self.chemistryModules.setCurrentIndex(self.chemistryModules.findText(self.get_available_chemistry(checkonly=True)))
            self.currentAddressTb.deselect()

    def filterListOfComp(self):
        text = self.findComp.text().upper()
        strict = False
        nonzeros = False
        sortByPsat = False
        if text != '':
            if text[-1] == '.':
                text = text[:-1]
                strict = True
            if '+' in text:
                nonzeros = True
                text = text.replace('+','')
            if '^' in text:
                text = text.replace('^','')
                sortByPsat = True

        self.availableVars.clear()
        if sortByPsat and self.LPD[0].masterfile == 'Particles.nc':
            names_ = array(self.LPD[0].names)[self.LPD[0].sorter]
        else:
            names_ = self.LPD[0].names
        if nonzeros:
            names_ = filter(lambda x: self.LPD[0].getconc(x).sum()>0, names_)

        if text == '':
            self.availableVars.addItems(names_)
            self.availableVars.item(0).setSelected(True)
            self.availableVars.setCurrentItem(self.availableVars.item(0))

        else:
            for c in names_:
                if strict:
                    if text == c.upper():
                        self.availableVars.addItem(c)
                        break
                else:
                    if text in c.upper():
                        self.availableVars.addItem(c)


    def filterListOfInput(self):
        for z in self.namesdat.selectedItems():
            z.setSelected(False)
        text = self.findInput.text().upper()
        strict = False
        if text != '':
            if text[-1] == '.':
                text = text[:-1]
                strict = True
        if text == '':
            self.namesdat.scrollToItem(self.namesdat.item(0), QtWidgets.QAbstractItemView.PositionAtTop)
        else:
            for c in NAMES:
                if strict:
                    if text == c.upper():
                        if (not namesPyInds[c] == divider_i) and (not namesPyInds[c] == divider_xtr_i):
                            self.namesdat.item(namesPyInds[c]).setSelected(True)
                else:
                    if text in c.upper():
                        if(not namesPyInds[c] == divider_i) and (not namesPyInds[c] == divider_xtr_i):
                            item = self.namesdat.item(namesPyInds[c])
                            item.setSelected(True)
                            self.namesdat.scrollToItem(item, QtWidgets.QAbstractItemView.PositionAtTop)

    def update_title(self):
        self.chamberLossBox.setTitle("Chamber properties (assuming square floor, volume: %0.2f m³)"
            %(self.floorArea.value()*self.chamberHeight.value()))


    def printHeaders(self):
        print('\n+----------------- Print input headers with column numbers -----------------+')
        for type, ff in zip(['ENV input', 'MCM input'], [nml.ENV.ENV_FILE, nml.MCM.MCM_FILE]):
            if exists(ff):
                f = open(ff)
            else:
                if ff=='':
                    print(type+' was not defined.')
                else:
                    print(type+' file "'+ff+'" was not found.')
                continue
            print()
            print('Header from '+ff+':')
            print()
            for line in f:
                for i,val in enumerate(line.split()):
                    if i==0 and '#' in val:
                        if val != '#':
                            print('   There appears to be some junk attached to the "#", maybe a column name? \n   Make sure there is space after the hashtag.')
                            break
                    elif i==0 and '#' not in val:
                        print('   < No header in file: '+ff+'>')
                        break
                    if i==0:
                        print('+-col-+-----column name-------')
                        continue
                    print(' %3d     %s'%(i, val))
                break
            f.close()
            print()
            print()

dummy = Comp()
defCompound = Comp()


def writeRATESdat(scheme):
    commandstring = [currentPythonVer, osjoin(gui_path, 'modules', 'parseRates.py'), osjoin('src','chemistry', scheme)]
    procs = Popen([*commandstring], stdout=PIPE,stderr=STDOUT,stdin=None)


def CustomCommandsCheatSheet():
    print("""
List of Custom commands; can be set in Run ARCA -> Custom model options
TYPE UNIT  NAME                       Options, (default), DESCRIPTION
----|-----|--------------------------|---------------------------------------
    """)
    with open('ModelLib/gui/conf/Custom_description.txt.csv', encoding="utf8", errors='ignore') as csvfile:
        f = csv.reader(csvfile, delimiter=',', quotechar='"')
        for k,line in enumerate(f):
            if k>0:
                print('\n%s'%(''.join(['_']*80)))
                print('')

            # if k>0: print('\n....|.....|..........................|.......................................')
            print('%-04s %-05s %-28s'%(line[0],line[1],line[2]), end='')
            fmt = '\n%38s %s'
            words = line[3].split()
            j=0
            for i in range(len(words)):
                j = j + len(words[i]) + 1
                if j<=40:
                    print(words[i],end=' ')
                else:
                    print('\n%38s %s'%(' ',words[i]),end=' ')
                    j=len(words[i]) + 1
    print('')

# with open("ModelLib/gui/conf/Custom_description.txt.csv") as f:
#     for line in f:
#         print(o,line.split(','))
        # print(' %04s  %04s %32s %s'%(line.split(',')[0],line.split(',')[1],line.split(',')[2],line.split(',')[3]))

if __name__ == '__main__':
    print(CurrentVersion+' started at:', ( time.strftime("%B %d %Y, %H:%M:%S", time.localtime())))
    app = QtWidgets.QApplication([sys.argv])
    qt_box = QtBoxGui()
    qt_box.setGeometry(30, 30, 900, 700)
    qt_box.show()
    styles = QtWidgets.QStyleFactory.keys()
    if "Fusion" in styles:
        app.setStyle(QtWidgets.QStyleFactory.create("Fusion"))
    else:
        print('Available styles: ',styles)
    app.exec_()
