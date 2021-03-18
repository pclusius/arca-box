#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
=============================================================================
Created By: Atmospheric modelling group AMG, Universities of Helsinki, Lund
and Saltzburg. To report bugs and/or make feature request contact
arca@helsinki.fi
=============================================================================
"""

from PyQt5 import QtCore, QtWidgets, QtGui, uic
import pyqtgraph as pg
from layouts import varWin,gui8,batchDialog1,batchDialog2,batchDialog3,vdialog,cc,about,input,t_editor
from modules import variations,vars,batch,GetVapourPressures as gvp
from subprocess import Popen, PIPE, STDOUT
from numpy import linspace,log10,sqrt,exp,pi,sin,shape,unique,array,ndarray,where,flip,zeros, sum as npsum
import numpy.ma as ma
from re import sub, finditer
from os import walk, mkdir, getcwd, chdir, chmod, environ, system, name as osname
from os.path import exists, dirname, getmtime, abspath, split as ossplit, join as osjoin, relpath as osrelpath
from shutil import copyfile as cpf
from re import sub,IGNORECASE, findall
import time
import pickle
from urllib.parse import urljoin

try:
    from scipy.ndimage.filters import gaussian_filter
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

QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling, True) # enable highdpi scaling
QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_UseHighDpiPixmaps, True)    # use highdpi icons
# See if scaling is necessary, currently only on Windows
if osname.upper() == 'NT' or operatingsystem == 'Windows': # this or is only till I sort out if platform usually works for people
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

## Some constants --------------------------------------------
# widths of the columns in "Input variables" tab
column_widths = [120,90,70,70,70,90,50,3]

# available units for variables, used to fill the tables and graphs with appropriate units
units = {
'TEMPK': ['K','C'],
'PRESSURE': ['Pa','hPa','mbar','kPa','bar','atm'],
'REL_HUMIDITY': ['%'],
'CONDENS_SINK':['1/s'],
'CON_SIN_NITR':['1/s'],
'SW_RADIATION':['W/m2'],
'ION_PROD_RATE':['ip/cm3 s'],
'NUC_RATE_IN':['1/cm3 s'],
'REST':['#/cm3','ppm','ppb','ppt','ppq']
}
# Name of the executable -------------------------------------------
exe_name = 'arcabox.exe'
# Path to variable names -------------------------------------------
path_to_names = 'ModelLib/required/NAMES.dat'
# GUI root
gui_path = 'ModelLib/gui/'
# GUI defaults are saved into this file. If it doesn't exist, it gets created in first start
defaults_file_path = gui_path+'conf/defaults.init'
minimal_settings_path = gui_path+'conf/minimal.init'
monitorfont_pickle = 'conf/monitorfont.pickle'
globalfont_pickle = 'conf/globalfont.pickle'

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
# BG colour for ENV vars
env_no = (215,238,244)
env_yes = (128,179,255)
# BG colour for ORG vars
org_no = (243,240,239)
org_yes = (172,147,147)

# icon
modellogo = gui_path+"/icons/ArcaLogo.png"
boxicon = gui_path+"/icons/thebox_ico.png"
CurrentVersion = "ARCA Box Model 0.9"
# Some messages
netcdfMissinnMes = ('Please note:',
'To open NetCDF-files you need netCDF4 for Python.\nYou can istall it with pip, package manager (or perhaps: python3 -m pip install --user netCDF4.')

# get current directory (to render relative paths) ----------
currentdir   = getcwd()
currentdir   = currentdir.replace('/ModelLib/gui', '')
currentdir_l = len(currentdir)
chdir(currentdir)

helpd = {}
with open(osjoin(gui_path,'conf','helplinks.txt'), 'r') as b:
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
AmmSystemFile = "src/ACDC/ACDC_module_ions_2018_08_31/Perl_input/input_ANnarrow_neutral_neg_pos.inp"
Amm_EnergyFile = "src/ACDC/ACDC_module_ions_2018_08_31/Perl_input/HS298.15K_426clusters2016Apr25.txt"
Amm_DipoleFile = "src/ACDC/ACDC_module_ions_2018_08_31/Perl_input/dip_pol_298.15K_426clusters2016Apr25.txt"
DMASystemFile = "src/ACDC/ACDC_module_2016_09_23/Perl_input/input_AD.inp"
DMA_EnergyFile = "src/ACDC/ACDC_module_2016_09_23/Perl_input/dH_dS.txt"
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
with open(path_to_names) as f:
    for line in f:
        name = line[:-1]
        if '#' in line:
            divider_i=i
            name = 'MCM compounds start here'
        NAMES.append(name)
        namesPyInds[name] = i
        namesFoInds[name] = i+1
        i += 1

## -----------------------------------------------------------
nml = vars.INITFILE(NAMES)
mmc = {'True': 'k', 'False': 0.95}
# The popup window for batch file preview
class batchW(QtGui.QDialog):
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
        self.ui.bDialogbuttonBox.helpRequested.connect(lambda: qt_box.helplink(helpd['batch']))

    def settext(self,a):
        """Setter for window text"""
        c = 1
        for i in range(3):
            if a[0][i]>0:
                exec('self.ui.label_%d.setText(a[1][%d])'%(c,i))
                exec('self.ui.tb_%d.appendPlainText(\'\'.join(a[2][%d]))'%(c,i))
                c +=1

# The popup window for About ARCA
class About(QtGui.QDialog):
    def __init__(self, parent = None):
        super(About, self).__init__(parent)
        self.ab = about.Ui_Dialog()
        self.ab.setupUi(self)
        self.ab.okgreat.clicked.connect(self.reject)
        self.ab.logo.setPixmap(QtGui.QPixmap(modellogo.replace('.png', 'HR.png')))

# The popup window for Create KPP files
class CCWin(QtGui.QDialog):
    def __init__(self, parent = None):
        super(CCWin, self).__init__(parent)
        self.ccw = cc.Ui_Dialog()
        self.ccw.setupUi(self)
        self.ccw.ccClose.clicked.connect(self.reject)
        self.ccw.createKPPsettings.clicked.connect(self.kpp)
        self.ccw.openOutput.clicked.connect(lambda: qt_box.openOutputDir(None, self.ccw.outDir.text()))
        self.ccw.browseOut.clicked.connect(lambda: qt_box.browse_path(self.ccw.outDir, 'dir'))
        self.ccw.browseSourceFile.clicked.connect(lambda: qt_box.browse_path(self.ccw.sourceFile, 'file'))
        self.ccw.browseIncludes.clicked.connect(lambda: qt_box.browse_path(self.ccw.includedFiles, 'append'))
        self.ccw.mainFrame.setFont(qt_box.font)
        self.ccw.manualCC.clicked.connect(lambda: qt_box.helplink(helpd['ccmanual']))
        self.ccw.pickDeffix.clicked.connect(self.ListFixed)

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
            includes.append(osjoin(ccloc,'PRAM_v21.txt'))
        out = osjoin(self.ccw.outDir.text(),'second.def')
        log = osjoin(self.ccw.outDir.text(),'second.log')
        commandstring = [currentPythonVer,ccloc+'/create_chemistry.py',cmds,'-o',out,'-l',log]
        if len(fixed)>0:
            commandstring.append('-d')
            commandstring += fixed
        if len(includes)>0:
            commandstring.append('-f')
            commandstring += includes

        print( 'Calling chemistry script with:\n'+' '.join(commandstring))
        self.kppProcess = Popen([*commandstring], stdout=PIPE,stderr=STDOUT,stdin=None)
        lines = True
        error = False
        warnings = False
        output = ['Chemistry definitions were created.\n']
        boilerplate = '\n1) Run KPP in the output directory: "kpp second.kpp"\n2) Recompile ARCA in tab "Chemistry".'
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

        if error: return
        cpf(osjoin(ccloc,'mcm_module.f90'),osjoin(self.ccw.outDir.text(),'mcm_module.f90'))
        cpf(osjoin(ccloc,'second.kpp'),osjoin(self.ccw.outDir.text(),'second.kpp'))
        if warnings: output.append('Read the .log and resolve the problems, then:')
        qt_box.popup('Chemistry created', '\n'.join(output)+boilerplate,0)


# The popup window for variations
class Variation(QtGui.QDialog):
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
        self.vary.manualVar.clicked.connect(lambda: qt_box.helplink(helpd['variations']))
        self.vary.mainFrame.setFont(qt_box.font)

    def vars(self):
        p=self.vary.lineEdit.text()
        if self.vary.table.rowCount() == 0: return
        ops = zeros((self.vary.table.rowCount(), 7))
        for i in range(self.vary.table.rowCount()):
            for j in range(7):
                if j==0:
                    try:
                        ops[i,j] = float(self.vary.table.item(i,j).text())
                    except:
                        if self.vary.table.item(i,j).text().upper() in self.indices:
                            ops[i,j] = self.indices[self.vary.table.item(i,j).text().upper()]
                        else:
                            print('Check the input, compound '+self.vary.table.item(i,j).text().upper()+' was not found.')
                else:
                    ops[i,j] = float(self.vary.table.item(i,j).text())
        print(variations.zzzz(p, ossplit(p)[0], ops, dryrun=False, nopause=True))

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
class VpressWin(QtGui.QDialog):
    def __init__(self, parent = None):
        super(VpressWin, self).__init__(parent)
        self.vp = vdialog.Ui_Dialog()
        self.vp.setupUi(self)
        self.vp.VapourClose.clicked.connect(self.reject)
        self.vp.UmanFrame.setEnabled(False)
        self.vp.useUMan.toggled.connect(lambda: qt_box.grayIfNotChecked(self.vp.useUMan,self.vp.UmanFrame))
        self.vp.massSmilesButton.clicked.connect(lambda: qt_box.browse_path(self.vp.lineEdit, 'file'))
        self.vp.PramButton.clicked.connect(lambda: qt_box.browse_path(self.vp.pramFile, 'file'))
        self.vp.browseVapourPath.clicked.connect(self.filename)
        self.vp.createVapourFileButton.clicked.connect(self.saveVapours)
        self.vp.UmanWWW.clicked.connect(lambda: qt_box.helplink(helpd['umanweb']))
        self.vp.manualVap.clicked.connect(lambda: qt_box.helplink(helpd['CreateVapourFile']))
        self.vp.mainFrame.setFont(qt_box.font)

    def filename(self):
        dialog = QtWidgets.QFileDialog()
        options = dialog.Options()
        options |= dialog.DontUseNativeDialog
        file = dialog.getSaveFileName(self, 'Save Vapours', options=options)[0]
        if file != '': self.vp.VapourPath.setText(file)

    def saveVapours(self):
        """This tool creates the Vapour and Elements files from user input or from the AMG server."""
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

        message = gvp.getVaps(args={
        'vp_method':vp_method[self.vp.vpCombo.currentIndex()],
        'bp_method':vp_method[self.vp.bpCombo.currentIndex()],
        'server':source,
        'smilesfile':self.vp.lineEdit.text(),
        'pram':self.vp.usePRAM.isChecked(),
        'pramfile':self.vp.pramFile.text(),
        'psat_lim':plim,
        'saveto':self.vp.VapourPath.text(),
        'filter':filterlist,
        'save_atoms':self.vp.saveElements.isChecked()
        })
        if len(message)==1: qt_box.popup('Oops...', 'No dice: '+message[0],3)
        if len(message)==2: qt_box.popup(message[0],message[1],0)

# The popup window for simple input
class Input(QtGui.QDialog):
    def __init__(self, parent = None, default = 0):
        super(Input, self).__init__(parent)
        self.inp = input.Ui_Dialog()
        self.inp.setupUi(self)
        self.inp.input.setValue(default)

# The popup window for editing simple text files
class Editor(QtGui.QDialog):
    def __init__(self, parent = None, file = None):
        super(Editor, self).__init__(parent)
        self.editor = t_editor.Ui_Dialog()
        self.editor.setupUi(self)
        self.setWindowTitle("Editing "+file)
        self.editor.editedText.setFont(qt_box.mfont)
        f = open(file, 'r')
        t = f.read()
        self.editor.editedText.appendPlainText(t)

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

class QtBoxGui(gui8.Ui_MainWindow,QtWidgets.QMainWindow):
    """Main program window."""
    def __init__(self):
        super(QtBoxGui,self).__init__()
        self.setupUi(self)

    # -----------------------
    # Common stuff
    # -----------------------
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
        self.actionARCA_webpage.triggered.connect(lambda: self.helplink(helpd['arcaweb']))
        self.actionOnline_manual.triggered.connect(lambda: self.helplink(helpd['manual']))
        self.actionFileHelp.triggered.connect(lambda: self.helplink(helpd['filehelp']))
        self.saveDefaults.clicked.connect(lambda: self.save_file(file=defaults_file_path))
        self.actionSave_as_defaults.triggered.connect(lambda: self.save_file(file=defaults_file_path))
        self.label_10.setPixmap(QtGui.QPixmap(modellogo))
        self.actionPrint_input_headers.triggered.connect(self.printHeaders)
        self.actionOpen_output_directory.triggered.connect(lambda: self.openOutputDir(None, self.currentAddressTb.text()))
        self.actionStopCurrentRunClean.triggered.connect(self.softStop)
        self.actionPrint_Custom_commands_cheat_sheet.triggered.connect(lambda: print(CustomCommandsCheatSheet))
        self.actionPlt_changes_from_current_dir.triggered.connect(self.plotChanges)

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
            elif i==divider_i:
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
        self.helpGroupName.clicked.connect(lambda: self.helplink(helpd['groupName']))
        self.helpRunName.clicked.connect(lambda: self.helplink(helpd['runName']))
        self.browseCommonIn.clicked.connect(lambda: self.browse_path(self.inout_dir, 'dir'))
        self.browseEnv.clicked.connect(lambda: self.browse_path(self.env_file, 'file'))
        self.browseMcm.clicked.connect(lambda: self.browse_path(self.mcm_file, 'file'))
        self.browsePar.clicked.connect(lambda: self.browse_path(self.dmps_file, 'file'))
        self.browseXtr.clicked.connect(lambda: self.browse_path(self.extra_particles, 'file'))
        # self.checkBox_aer.stateChanged.connect(lambda: self.grayIfNotChecked(self.checkBox_aer,self.groupBox_8))
        self.checkBox_aer.toggled.connect(lambda: self.grayIfNotChecked(self.checkBox_aer,None, 4))
        self.checkBox_che.stateChanged.connect(lambda: self.grayIfNotChecked(self.checkBox_che,None, 2))
        self.fsave_division.valueChanged.connect(self.toggle_printtime)
        self.checkBox_acd.stateChanged.connect(lambda: self.grayIfNotChecked(self.checkBox_acd,None, hide=3))
        # self.simpleton.stateChanged.connect(lambda: self.grayIfNotChecked(self.simpleton,None, hide=6))
        # self.simpleton.stateChanged.connect(lambda: self.grayIfNotChecked(self.simpleton,None, hide=5))
        # self.use_dmps.stateChanged.connect(lambda: self.grayIfNotChecked(self.use_dmps,self.dmps_read_in_time))
        self.dateEdit.dateChanged.connect(self.updatePath)
        self.indexEdit.valueChanged.connect(self.updatePath)
        self.case_name.textChanged.connect(self.updatePath)
        self.case_name.editingFinished.connect(lambda: self.allcaps(self.case_name))
        self.run_name.textChanged.connect(self.updatePath)
        self.run_name.editingFinished.connect(lambda: self.allcaps(self.run_name))
        self.inout_dir.textChanged.connect(self.updatePath)
        self.indexRadioDate.toggled.connect(self.updatePath)
        # self.useSpeed.stateChanged.connect(lambda: self.grayIfNotChecked(self.useSpeed,self.precLimits))
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

        self.batchHelp.clicked.connect(lambda: self.helplink(helpd['batch']))
        self.PrecHelp.clicked.connect(lambda: self.helplink(helpd['precision']))

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
    # tab Losses
    # -----------------------
        self.browseLosses.clicked.connect(lambda: self.browse_path(self.losses_file, 'file'))

    # -----------------------
    # tab Chemistry
    # -----------------------

        self.kppTool.clicked.connect(self.createCC)
        self.vapTool.clicked.connect(self.vapours)
        self.recompile.clicked.connect(self.remake)
        self.TimerCompile = QtCore.QTimer(self);
        self.TimerCompile.timeout.connect(self.progress)
        self.compileProgressBar.hide()
        self.running = 0
        self.get_available_chemistry(checkonly=False)
        self.chemLabel.setText('Current chemistry scheme in makefile: '+self.get_available_chemistry(checkonly=True))

    # -----------------------
    # tab Aerosols
    # -----------------------
        # self.frameBase.setEnabled(False)
        self.butVapourNames.clicked.connect(lambda: self.browse_path(self.vap_names, 'file'))
        self.butVapourAtoms.clicked.connect(lambda: self.browse_path(self.vap_atoms, 'file'))
        self.use_atoms.stateChanged.connect(lambda: self.grayIfNotChecked(self.use_atoms,self.vap_atoms))
        self.use_atoms.stateChanged.connect(lambda: self.grayIfNotChecked(self.use_atoms,self.butVapourAtoms))
        self.Org_nucl.toggled.connect(lambda: self.grayIfChecked(self.Org_nucl,self.resolve_base))
        self.resolve_base.toggled.connect(lambda: self.grayIfChecked(self.resolve_base,self.groupBox_16))
        # self.resolve_base.stateChanged.connect(lambda: self.grayIfNotChecked(self.resolve_base,self.frameBase))
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

        self.EditAmmSystem.clicked.connect(lambda: self.editTxtFile(AmmSystemFile))
        self.EditAmm_E.clicked.connect(lambda: self.editTxtFile(Amm_EnergyFile))
        self.EditAmmDipoles.clicked.connect(lambda: self.editTxtFile(Amm_DipoleFile))
        self.EditDMASystem.clicked.connect(lambda: self.editTxtFile(DMASystemFile))
        self.EditDMA_E.clicked.connect(lambda: self.editTxtFile(DMA_EnergyFile))

    # -----------------------
    # tab Process Monitor
    # -----------------------
        self.currentEndLine = 0
        self.fulltext = ''
        self.tabWidget.currentChanged.connect(self.activeTab)
        self.tabWidget_4.currentChanged.connect(self.activeSubTab)
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
        self.pollTimer.timeout.connect(self.updateOutput)
        self.pauseScroll.clicked.connect(lambda: self.MonitorWindow.verticalScrollBar().setSliderPosition(self.MonitorWindow.verticalScrollBar().maximum()))
        if exists(osjoin(gui_path,monitorfont_pickle)):
            try:
                sf = pickle.load(open(osjoin(gui_path,monitorfont_pickle), "rb"))
                self.mfont = self.MonitorWindow.font()
                self.mfont.setFamily(sf[0])
                self.mfont.setPointSize(sf[1])
                self.mfont.setBold(sf[2])
                self.mfont.setItalic(sf[3])
                self.MonitorWindow.setFont(self.mfont)
            except:
                pass
        else:
            self.mfont = self.MonitorWindow.font()
            savefont = [self.mfont.family(),self.mfont.pointSize(),self.mfont.bold(),self.mfont.italic()]
            pickle.dump(savefont, open(osjoin(gui_path,monitorfont_pickle), 'wb'))

        if exists(osjoin(gui_path,globalfont_pickle)):
            try:
                sf = pickle.load(open(osjoin(gui_path,globalfont_pickle), "rb"))
                self.font = self.centralwidget.font()
                self.font.setFamily(sf[0])
                self.font.setPointSize(sf[1])
                self.font.setBold(sf[2])
                self.font.setItalic(sf[3])
                self.centralwidget.setFont(self.font)
            except:
                pass
        else:
            self.font = self.centralwidget.font()
            savefont = [self.font.family(),self.font.pointSize(),self.font.bold(),self.font.italic()]
            pickle.dump(savefont, open(osjoin(gui_path,globalfont_pickle), 'wb'))
        self.viewPrintNML.clicked.connect(lambda: self.editTxtFile(SCREENPRINT_NML))
    # -----------------------
    # tab Output Graph
    # -----------------------
        if netcdf:
            self.show_netcdf.hide()
            self.show_netcdf_2.hide()
            self.fLog_2.clicked.connect(self.showOutputUpdate)
            self.fLin_2.clicked.connect(self.showOutputUpdate)
            self.findComp.textChanged.connect(self.filterListOfComp)
            self.loadNetcdf.clicked.connect(lambda: self.browse_path(None, 'plot', ftype="NetCDF (*.nc)"))
            self.loadNetcdf_mass.clicked.connect(lambda: self.browse_path(None, 'plot_mass', ftype="ARCA particle file (Particles.nc)"))
            self.loadNetcdfPar.clicked.connect(lambda: self.browse_path(None, 'plotPar', ftype="NetCDF, sum (*.nc *.sum *.dat)",plWind=0))
            self.loadSumPar.clicked.connect(lambda: self.browse_path(None, 'plotPar', ftype="NetCDF, sum (*.nc *.sum *.dat)",plWind=1))
            self.CloseLinePlotsButton.clicked.connect(self.closenetcdf)
        else:
            self.sumSelection.setEnabled(False)
            self.show_netcdf.show()
            self.fLin_2.setEnabled(False)
            self.fLog_2.setEnabled(False)
            self.findComp.setEnabled(False)
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
    # Load preferences, or create preferences if not found
    # -----------------------
        try:
            # if defaults exist, use them
            self.load_initfile(defaults_file_path)
        except:
            try:
                # If defaults did not exist, load minimal working settings
                self.load_initfile(minimal_settings_path)
            except:
                # If they also were missing, use "factory settings"
                pass
            # Save the obtained settings as default
            self.save_file(file=defaults_file_path, mode='silent')
        self.updateEnvPath()

    # -----------------------
    # Class methods
    # -----------------------
    def openOutputDir(self, a, dir):
        import os
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
            QtGui.QDesktopServices.openUrl(QtCore.QUrl(linkStr))

    def activeTab(self,i):
        pass
        # if i == 6:
        #     # self.MonitorWindow.verticalScrollBar().setSliderPosition(self.currentEndLine)
        #     self.MonitorWindow.verticalScrollBar().setSliderPosition(self.MonitorWindow.verticalScrollBar().maximum()-40)
        #     self.MonitorWindow.verticalScrollBar().setSliderPosition(self.MonitorWindow.verticalScrollBar().maximum())

    def activeSubTab(self,i):
        pass
        # if i == 0: self.MonitorWindow.verticalScrollBar().setSliderPosition(self.currentEndLine)

    def guiSetFont(self, wdgt, name, reset=False):
        if not reset:
            dialog = QtWidgets.QFontDialog()
            font, ok = dialog.getFont(wdgt.font(), parent=self)
            if ok:
                wdgt.setFont(font)
                savefont = [font.family(),font.pointSize(),font.bold(),font.italic()]
                pickle.dump(savefont, open(osjoin(gui_path,globalfont_pickle.replace('global', name)), 'wb'))
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
            savefont = [font.family(),font.pointSize(),font.bold(),font.italic()]
            pickle.dump(savefont, open(osjoin(gui_path,globalfont_pickle.replace('global', name)), 'wb'))

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
        x = 10**linspace(log10(x0),log10(x1),nb)
        acl = zeros(len(x))
        k = log10(x[1]/x[0])
        for i in range(len(luvut)//3):
            if abs(luvut[3*i+1])>0:
                acl = acl + luvut[3*i+2]*(gaussian(log10(x), log10(luvut[3*i+0]),luvut[3*i+1]))
        if sum( acl )>0:
            Z = N * acl / (sum( acl ))
            ndel = -min(nb-1, 8)
            Z[ndel:] = where(Z[ndel:]>1e-12, 0,Z[ndel:])
            self.HPLotter.plot(x,Z/k,pen=pg.mkPen('w', width=4), clear=True, name='PSD')
            self.HPLotter.setLogMode(x=True)
            self.HPLotter.showGrid(x=True,y=True)


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
        self.saveCurrentButton.setStyleSheet("background-image: url('%s'); background-repeat: no-repeat" %(gui_path+"icons/save.png"))
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
        self.update_nml()
        self.env_file.setToolTip('Location: "'+nml.ENV.ENV_FILE+'"')
        self.mcm_file.setToolTip('Location: "'+nml.MCM.MCM_FILE+'"')
        self.dmps_file.setToolTip('Location: "'+nml.PARTICLE.DMPS_FILE+'"')
        self.extra_particles.setToolTip('Location: "'+nml.PARTICLE.EXTRA_PARTICLES+'"')


    def get_case_kwargs(self,r):
        return {'begin':r[0],'end':r[1],'case':nml.PATH.CASE_NAME,'run':nml.PATH.RUN_NAME, 'common_root':nml.PATH.INOUT_DIR}


    def updatePath(self):
        if self.fileLoadOngoing:
            return
        self.update_nml()

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
        self.update_nml()
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
                bf.write('./'+exe_name+' '+file+' |tee '+dirname(dirname(file)[:-1])+'/'+nml.PATH.RUN_NAME+'/RunReport.txt'+'\n' )
            if self.batchRangeDay.isChecked():
                nml.TIME.DATE='%s'%(date)
                nml.TIME.INDEX=''
            else:
                nml.TIME.DATE=''
                nml.TIME.INDEX='%s'%(date)

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


    def radio(self,*buts):
        for but in buts:
            if but.isChecked():
                if but.text() == 'Linear':
                    self.PLOT.setLogMode(y=False)
                    return 'lin'
                if but.text() == 'Logarithmic':
                    self.PLOT.setLogMode(y=True)
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

        if 'dir' in mode:
            dialog.setFileMode(QtWidgets.QFileDialog.Directory)
            dialog.setOption(QtWidgets.QFileDialog.ShowDirsOnly)
        options = dialog.Options()
        options |= dialog.DontUseNativeDialog
        if mode == 'dir':
            path = dialog.getExistingDirectory(self, 'Choose Directory', options=options)
        elif mode == 'export':
            path = dialog.getSaveFileName(self, 'Save INITFILE', options=options)[0]
        else:
            dialog.setOption(QtWidgets.QFileDialog.DontUseNativeDialog)
            dialog.setWindowTitle('Choose File')
            if dialog.exec() == 1:
                path = dialog.selectedFiles()[0]
            else: path=''
        if path != '':
            if path[:currentdir_l] == currentdir and path[currentdir_l] == '/':
                path = path[currentdir_l+1:]
            if mode == 'load':
                self.load_initfile(path)
                self.show_currentInit(path)
            elif mode == 'fixed':
                self.loadFixedFile(path)
            elif mode == 'fixed_cc':
                return path
            elif mode == 'plot':
                self.showOutput(path)
            elif mode == 'plot_mass':
                self.showMass(path)
            elif mode == 'plotPar':
                self.showParOutput(path,plWind)
            elif mode == 'batchList':
                self.ListbatchCaller(path)
            elif mode == 'dironly':
                path = path[path.rfind('/')+1:]
                target.clear()
                target.insert(path)
            elif mode == 'export':
                self.exportCurrentCase(path)
            elif mode == 'append':
                target.appendPlainText(path)
            else:
                target.clear()
                target.insert(path)


    def save_file(self, file=None, mode=None,nobatch=True, changeTitle=True):
        if nobatch:
            self.update_nml()
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
        pattern = '(CHMDIR)( )*(=)( )*[a-z|A-Z|0-9|-|_]*'
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
        self.edwin = Editor(file=file)
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
        unit.addItems(units.get(name,units['REST']))
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
        if guard.isChecked() == True:
            frame.setChecked(0)
            frame.setEnabled(False)
        else:
            frame.setEnabled(True)


    def toggle_printtime(self):
        if self.fsave_division.value() != 0 :
            self.fsave_interval.setEnabled(False)
        else:
            self.fsave_interval.setEnabled(True)


    def print_values(self, file=None, mode=None, nobatch=True):
        if nobatch:
            self.update_nml()
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
        try: self.tabWidget.currentChanged.disconnect(self.activeTab)
        except: pass

        if exists(self.saveCurrentOutputDir):
            f = open(osjoin(self.saveCurrentOutputDir,'runReport.txt'), 'w')
            f.write(self.MonitorWindow.toPlainText())
            f.close()
        else:
            return

    def softStop(self,a):
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
                self.scaling0 = (time[-1]/n.shape[0],log10(diam[-1]/diam[0])/(len(diam)-1), log10(diam[0]*1e9))
            else:
                self.scaling0 = (1,log10(diam[-1]/diam[0])/(len(diam)-1), log10(diam[0]*1e9))
        if windowInd==1:
            self.z1 = n
            if self.X_axis_in_time.isChecked():
                self.scaling1 = (time[-1]/n.shape[0],log10(diam[-1]/diam[0])/(len(diam)-1), log10(diam[0]*1e9))
            else:
                self.scaling1 = (1,log10(diam[-1]/diam[0])/(len(diam)-1), log10(diam[0]*1e9))
        self.drawSurf(window, new=1)


    def drawSurf(self,window, new=0):
        use_filter = False
        if window==self.surfacePlotWindow_0:
            n_levelled = self.z0
            scale = self.scaling0
            if self.Filter_0.isChecked():
                use_filter = True
        else:
            n_levelled = self.z1
            scale = self.scaling1
            if self.Filter_1.isChecked():
                use_filter = True

        levels=(self.lowlev.value(),self.highlev.value())
        if scipyIs and use_filter: n_levelled = gaussian_filter(n_levelled,(self.gauss_x.value(),self.gauss_y.value()),mode='constant')
        n_levelled = where(n_levelled>=levels[1],levels[1]*0.98,n_levelled)

        hm = pg.ImageItem(n_levelled)
        cb = ndarray((20,1))
        cb[:,0] = linspace(self.lowlev.value(), self.highlev.value(), 20)
        ss = pg.ImageItem(cb.T)
        ss.translate(0,self.lowlev.value())
        ss.scale(1,(self.highlev.value()-self.lowlev.value())/19)
        if self.Y_axis_in_nm.isChecked():
            hm.translate(0,scale[2])
            hm.scale(scale[0],scale[1])

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
        if self.Y_axis_in_nm.isChecked():
            window.setLogMode(False, True)
        else:
            window.setLogMode(False, False)

    def plotChanges(self):
        import modules.optich as optich
        ret = optich.plot(self.currentAddressTb.text())
        if ret!='': self.popup('Plot changes:',ret,0)

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
        self.print_values(tempfile)
        self.wait_for = currentWait

        try:
            self.boxProcess = Popen(["./"+exe_name, "%s"%tempfile, '--gui'], stdout=PIPE,stderr=STDOUT,stdin=None)
            self.saveCurrentOutputDir = self.currentAddressTb.text()
            self.MonitorWindow.clear()
            self.Timer.start(5)
            # self.TimerPlot.start(1000)
            self.pollTimer.start(2000)
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
        nml.SETTINGS.INPUT = '%s:%s %s:%s %s:%s %s:%s %s:%s %s:%s %s:%s %s:%s' %(
                                'env_file', self.env_file.text(),
                                'mcm_file', self.mcm_file.text(),
                                'dmps_file', self.dmps_file.text(),
                                'extra_particles', self.extra_particles.text(),

                                'stripRoot_env', self.checkboxToFOR(self.stripRoot_env),
                                'stripRoot_mcm', self.checkboxToFOR(self.stripRoot_mcm),
                                'stripRoot_par', self.checkboxToFOR(self.stripRoot_par),
                                'stripRoot_xtr', self.checkboxToFOR(self.stripRoot_xtr),
        )

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
        nml.FLAG.MODEL_H2SO4=self.checkboxToFOR(self.model_h2so4)
        nml.FLAG.ORG_NUCL=self.checkboxToFOR(self.Org_nucl)
        nml.FLAG.RESOLVE_BASE=self.checkboxToFOR(self.resolve_base)
        nml.FLAG.PRINT_ACDC=self.checkboxToFOR(self.print_acdc)
        nml.FLAG.USE_SPEED=self.checkboxToFOR(self.useSpeed)
        nml.FLAG.AFTER_CHEM_ON=self.checkboxToFOR(self.after_chem_on)
        nml.FLAG.AFTER_NUCL_ON=self.checkboxToFOR(self.after_nucl_on)

        # class _TIME:
        nml.TIME.RUNTIME=self.runtime.value()
        nml.TIME.DT=self.dt.value()
        nml.TIME.FSAVE_INTERVAL=self.fsave_interval.value()
        nml.TIME.PRINT_INTERVAL=self.print_interval.value()
        nml.TIME.FSAVE_DIVISION=self.fsave_division.value()
        if self.indexRadioDate.isChecked():
            self.index_for_parser = self.dateEdit.text()
            nml.TIME.DATE=self.dateEdit.text()
            nml.TIME.INDEX=''
        if self.indexRadioIndex.isChecked():
            self.index_for_parser = self.indexEdit.value()
            nml.TIME.DATE=''
            nml.TIME.INDEX='%04d'%self.indexEdit.value()

        # class _PARTICLE:
        nml.PARTICLE.PSD_MODE=self.psd_mode.currentIndex()
        nml.PARTICLE.N_BINS_PAR=self.n_bins_particle.value()
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

        # class _ENV:
        nml.ENV.ENV_FILE=self.pars(self.env_file.text(), file=self.indir, stripRoot=self.stripRoot_env.isChecked())
        nml.ENV.LOSSES_FILE=self.pars(self.losses_file.text(), file=self.indir, stripRoot=False)
        nml.ENV.CHAMBER_FLOOR_AREA=self.floorArea.value()
        nml.ENV.CHAMBER_CIRCUMFENCE=self.chamberCircumfence.value()
        nml.ENV.CHAMBER_HEIGHT=self.chamberHeight.value()
        nml.ENV.EDDYK=self.eddyK.value()
        nml.ENV.USTAR=self.ustar.value()
        nml.ENV.ALPHAWALL=self.alphaWall.value()

        # class _MCM:
        nml.MCM.MCM_FILE=self.pars(self.mcm_file.text(), file=self.indir, stripRoot=self.stripRoot_mcm.isChecked())

        # class _MISC:
        nml.MISC.LAT=self.lat.value()
        nml.MISC.LON=self.lon.value()
        nml.MISC.WAIT_FOR=self.wait_for
        nml.MISC.DESCRIPTION=self.description.toPlainText().replace('\n','<br>')
        nml.MISC.CH_ALBEDO=self.ch_albedo.value()
        nml.MISC.DMA_F=self.dma_f.value()
        nml.MISC.RESOLVE_BASE_PRECISION=self.resolve_base_precision.value()
        nml.MISC.FILL_FORMATION_WITH=self.resolveHelper()
        nml.MISC.SKIP_ACDC=self.checkboxToFOR(self.skip_acdc)
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
        return


    def resolveHelper(self):
        text = self.fill_formation_with.currentText()
        if text == 'Fixed ratio':
            return ''
        elif 'NH3' in text:
            return 'NH3'
        elif 'DMA' in text:
            return 'DMA'


    def load_initfile(self,file):
        if not exists(file):
            if file != defaults_file_path:
                self.popup('Ooops', 'File "%s" not found'%file, icon=3)
                return
        self.fileLoadOngoing = True
        self.markReverseSelection('all')
        self.remv_item()
        if self.plotTo.isChecked() == True:
            self.plotTo.setChecked(False)
            self.show_extra_plots = ''
            self.updteGraph(first=True)

        def solve_for_parser(query):
            if query.upper() == 'NH3': return 2
            elif query.upper() == 'DMA': return 1
            else: return 0

        def parse_date(str):
            if len(str)==10:
                y = int(str[0:4])
                m = int(str[5:7])
                d = int(str[8:10])
                return QtCore.QDate(y, m, d)
            else:
                return QtCore.QDate(2000, 1, 1)

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
                    if index-1 == divider_i:
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
            elif 'MODEL_H2SO4' == key: self.model_h2so4.setChecked(strng)
            elif 'ORG_NUCL' == key: self.Org_nucl.setChecked(strng)
            elif 'RESOLVE_BASE' == key: self.resolve_base.setChecked(strng)
            elif 'RUNTIME' == key and isFl: self.runtime.setValue(float(strng))
            elif 'DT' == key and isFl: self.dt.setValue(float(strng)),
            elif 'PRINT_ACDC' == key: self.print_acdc.setChecked(strng)
            elif 'USE_SPEED' == key: self.useSpeed.setChecked(strng)
            elif 'AFTER_CHEM_ON' == key: self.after_chem_on.setChecked(strng)
            elif 'AFTER_NUCL_ON' == key: self.after_nucl_on.setChecked(strng)
            elif 'FSAVE_INTERVAL' == key and isFl: self.fsave_interval.setValue(int(strng))
            elif 'PRINT_INTERVAL' == key and isFl: self.print_interval.setValue(int(strng))
            elif 'FSAVE_DIVISION' == key and isFl: self.fsave_division.setValue(int(strng))
            elif 'DATE' == key: self.dateEdit.setDate(parse_date(strng))
            elif 'INDEX' == key and isFl: self.indexEdit.setValue(int(strng))
            elif 'PSD_MODE' == key and isFl: self.psd_mode.setCurrentIndex(int(strng))
            elif 'N_BINS_PAR' == key and isFl: self.n_bins_particle.setValue(int(strng))
            elif 'MIN_PARTICLE_DIAM' == key: self.min_particle_diam.setText(strng)#   1.0000000000000001E-009,
            elif 'MAX_PARTICLE_DIAM' == key: self.max_particle_diam.setText(strng)#   9.9999999999999995E-007,
            elif 'N_MODAL' == key: self.n_modal.setText(strng)#   9.9999999999999995E-007,
            elif 'MMODAL_INPUT_INUSE' == key: self.multiModalBox.setChecked(int(strng))#   9.9999999999999995E-007,
            elif 'DMPS_FILE' == key: self.dmps_file.setText(strng)
            elif 'EXTRA_PARTICLES' == key: self.extra_particles.setText(strng)
            elif 'MMODAL_INPUT' == key: self.mmodal_input.setText(strng)
            elif 'DMPS_READ_IN_TIME' == key and isFl: self.dmps_read_in_time.setValue(float(strng))
            elif 'DMPS_HIGHBAND_LOWER_LIMIT' == key: self.dmps_highband_lower_limit.setText(strng)
            elif 'DMPS_LOWBAND_UPPER_LIMIT' == key: self.dmps_lowband_upper_limit.setText(strng)
            elif 'USE_DMPS' == key: self.use_dmps.setChecked(strng)
            elif 'USE_DMPS_PARTIAL' == key: self.use_dmps_partial.setChecked(strng)
            elif 'ENV_FILE' == key: self.env_file.setText(strng)
            elif 'CHAMBER_FLOOR_AREA' == key and isFl: self.floorArea.setValue(float(strng))#  0.20000000000000001     ,
            elif 'CHAMBER_CIRCUMFENCE' == key and isFl: self.chamberCircumfence.setValue(float(strng))#  0.20000000000000001     ,
            elif 'CHAMBER_HEIGHT' == key and isFl: self.chamberHeight.setValue(float(strng))#  0.20000000000000001     ,
            elif 'EDDYK' == key and isFl: self.eddyK.setValue(float(strng))#  0.05000000000000001     ,
            elif 'USTAR' == key and isFl: self.ustar.setValue(float(strng))#  0.050000000000000001     ,
            elif 'ALPHAWALL' == key and isFl: self.alphaWall.setValue(float(strng))#  0.050000000000000001     ,
            elif 'MCM_FILE' == key: self.mcm_file.setText(strng)# "
            elif 'LOSSES_FILE' == key: self.losses_file.setText(strng)# "
            elif 'LAT' == key and isFl: self.lat.setValue(float(strng))
            elif 'LON' == key and isFl: self.lon.setValue(float(strng))
            elif 'WAIT_FOR' == key and isFl: self.wait_for = (int(strng))
            elif 'DESCRIPTION' == key: self.description.setPlainText(strng.replace('<br>','\n'))# "Just some keying
            elif 'CH_ALBEDO' == key and isFl: self.ch_albedo.setValue(float(strng))#  0.20000000000000001     ,
            elif 'DMA_F' == key and isFl: self.dma_f.setValue(float(strng))
            elif 'RESOLVE_BASE_PRECISION' == key and isFl: self.resolve_base_precision.setValue(float(strng))
            elif 'FILL_FORMATION_WITH' == key: self.fill_formation_with.setCurrentIndex(solve_for_parser(strng))
            elif 'SKIP_ACDC' == key: self.skip_acdc.setChecked(strng)
            elif 'USE_ATOMS' == key: self.use_atoms.setChecked(strng)
            elif 'VAP_NAMES' == key: self.vap_names.setText(strng)
            elif 'VAP_ATOMS' == key: self.vap_atoms.setText(strng)
            elif 'GR_SIZES' == key: self.GR_sizes.setText(strng)

            elif 'DDIAM_RANGE' == key:
                self.prec_low_1.setValue(float(strng.split(',')[0]))
                self.prec_high_1.setValue(float(strng.split(',')[1]))
            elif 'DPNUM_RANGE' == key:
                self.prec_low_2.setValue(float(strng.split(',')[0]))
                self.prec_high_2.setValue(float(strng.split(',')[1]))
            elif 'DVAPO_RANGE' == key:
                self.prec_low_3.setValue(float(strng.split(',')[0]))
                self.prec_high_3.setValue(float(strng.split(',')[1]))

            elif in_custom:
                nml.CUSTOM.CUSTOMS.append([key, strng])
                pass

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
            unts = units.get(key,units['REST'])
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
        self.updatePath()
        self.updateEnvPath()


    def updateMass(self):
        if (self.ncs_mass != 0):
            self.showMass(first=False, target='mass')


    def updateNumbers(self):
        if (self.ncs_mass != 0):
            self.showMass(first=False, target='numb')


    def showMass(self, file=None, first=True, target=None):
        if first:
            if file==None: return
            self.closenetcdf_mass()
            if exists(file):
                try:
                    self.ncs_mass = netCDF4.Dataset(file, 'r')
                except:
                    self.popup('Bummer...', 'Not a valid output file',icon=3)
                    return
            DIAMETER = self.ncs_mass.variables['DIAMETER'][:]
            self.diams.addItems(['%7.2f'%(1e9*i) for i in DIAMETER[0,:]])
            self.diams.setSelectionMode(QtWidgets.QAbstractItemView.ExtendedSelection)
            self.diams.itemSelectionChanged.connect(self.updateMass)
            self.diams.selectAll()
            time = self.ncs_mass.variables['time_in_hrs'][:]
            self.times.itemSelectionChanged.connect(self.updateNumbers)
            self.times.addItems(['%7.2f'%(i) for i in time])
            self.times.setSelectionMode(QtWidgets.QAbstractItemView.ExtendedSelection)
            self.times.item(0).setSelected(True)
            indstime = [0]
            self.plotResultWindow_3.setLogMode(x=True)
            self.plotResultWindow_2.setLabel('bottom', 'Time', units='h')
            self.plotResultWindow_2.setLabel('left', 'Mass', units='g')

            self.plotResultWindow_3.setLabel('bottom', 'Diameter', units='m')
            self.plotResultWindow_3.setLabel('left', '# normalized', units=None)

        else:
            if target=='mass':
                self.plotResultWindow_2.clear()
            if target=='numb':
                self.plotResultWindow_3.clear()
            indstime = [i.row() for i in self.times.selectedIndexes()]
        inds = [i.row() for i in self.diams.selectedIndexes()]

        DIAMETER                = self.ncs_mass.variables['DIAMETER'][:]
        NUMBER_CONCENTRATION    = self.ncs_mass.variables['NUMBER_CONCENTRATION'][:]
        MASS_OF_SINGLE_PAR      = self.ncs_mass.variables['MASS'][:]
        mass_in_bin             = MASS_OF_SINGLE_PAR*NUMBER_CONCENTRATION
        lognormconc             = NUMBER_CONCENTRATION/log10(DIAMETER[0,1]/DIAMETER[0,0])
        self.measdmps = False
        if self.showAlsoMeasInMassConc.isChecked():
            try:
                DMPS_CONCENTRATION = self.ncs_mass.variables['INPUT_CONCENTRATION'][:]
                massdmps = MASS_OF_SINGLE_PAR*DMPS_CONCENTRATION
                lognormdmps = DMPS_CONCENTRATION/log10(DIAMETER[0,1]/DIAMETER[0,0])
                self.measdmps = True
            except:
                print('File did not contain measured PSD')

        time = self.ncs_mass.variables['time_in_hrs'][:]
        y  = npsum(mass_in_bin[:,inds],axis=1)*1e3
        miny, maxy = y.min(),y.max()
        if self.measdmps:
            y2 = npsum(massdmps[:,inds],axis=1)*1e3
            miny, maxy = min(miny,y2.min()),max(maxy,y2.max())
        if maxy>0:
            if abs(1-miny/maxy) <1e-12:
                maxy = miny*100000
        if self.measdmps and self.showAlsoMeasInMassConc.isChecked():
            self.outplot_mass = self.plotResultWindow_2.plot(time,
                                                            y2,
                                                            pen={'color':'r','width': 2.0,'style': QtCore.Qt.DotLine},
                                                            symbol='x',
                                                            symbolPen='r',
                                                            symbolSize=6
                                                            )
        self.outplot_mass = self.plotResultWindow_2.plot(time,
                                                        y,
                                                        pen={'color':'b','width': 2.0},
                                                        symbol='o',
                                                        symbolPen='b',
                                                        symbolSize=3
                                                        )

        if target == 'mass' or first: self.plotResultWindow_2.setRange(yRange=[miny*0.95,maxy*1.05])

        for i,ii in enumerate(indstime):
            if self.measdmps and self.showAlsoMeasInMassConc.isChecked():
                self.outplot_numb2 = self.plotResultWindow_3.plot(DIAMETER[0,:],
                                                                lognormdmps[ii,:],
                                                                pen={'color':colors[i%10],'width': 2.0,'style': QtCore.Qt.DotLine},
                                                                symbol='x',
                                                                symbolPen=colors[i%10],
                                                                symbolSize=6
                                                                )
            self.outplot_numb = self.plotResultWindow_3.plot(DIAMETER[0,:],
                                                            lognormconc[ii,:],
                                                            pen={'color':colors[i%10],'width': 2.0},
                                                            symbol='o',
                                                            symbolPen='k',
                                                            symbolSize=3
                                                            )
#

    def showOutput(self, file, add=False):
        # First set plot mode to linear in order to avoid errors with zero values
        self.fLin_2.setChecked(True)
        self.fLog_2.setChecked(False)
        self.plotResultWindow.setLogMode(y=False)
        self.findComp.clear()
        # Close all previus netcdf-files and clear plot
        self.closenetcdf()
        # Try to open netCDF-file
        if exists(file):
            try:
                self.ncs = netCDF4.Dataset(file, 'r')
            except:
                self.popup('Bummer...', 'Not a valid output file',icon=3)
                return

        # find out the time dimension, using unlimited dimension here
        for timedim in self.ncs.dimensions:
            if self.ncs.dimensions[timedim].isunlimited():
                break

        checker = lambda v,n: v.lower() in timedim and 'Shifter' not in n and 'Multipl' not in n

        # collect all variables and dimensions from netCDF-dataset to hnames
        cache = array([i.name for i in self.ncs.get_variables_by_attributes(ndim=1)])
        timevars = [checker(i.dimensions[0], i.name) for i in self.ncs.get_variables_by_attributes(ndim=1)]
        self.hnames = cache[timevars]

        # Now try to plot first the third line, which is the first non-time variable
        try:
            time = self.ncs.variables[self.hnames[0]][:]/3600
            vari = self.ncs.variables[self.hnames[2]][:]
            if ma.is_masked(time):
                self.outplot = self.plotResultWindow.plot(
                time[~time.mask],
                vari[~time.mask],
                pen={'color':'b','width': 2.0}
                )
            else:
                self.outplot = self.plotResultWindow.plot(
                time,
                vari,
                pen={'color':'b','width': 2.0}
                )
            self.availableVars.clear()
            self.availableVars.addItems(self.hnames)
            self.availableVars.item(2).setSelected(True)
            self.availableVars.itemSelectionChanged.connect(self.showOutputUpdate)

        # If fails, give information and return
        except:
            self.popup('Bummer...', "Output file does not contain any plottable data",icon=3)
            return

        # All's well, finish plot
        self.plotResultWindow.setLabel('bottom', 'Time', units='h')
        self.plotTitle = file + ': ' + self.hnames[2]+' ['+units.get(self.hnames[2],units['REST'])[0]+']'
        self.plotResultTitle.setText(self.plotTitle)


    def selectionMode(self):
        if self.sumSelection.isChecked():
            self.availableVars.setSelectionMode(QtWidgets.QAbstractItemView.ExtendedSelection)
        else:
            self.availableVars.setSelectionMode(QtWidgets.QAbstractItemView.SingleSelection)
            if self.availableVars.currentItem() != None:
                self.availableVars.currentItem().setSelected(True)


    def showOutputUpdate(self):
        """This function is inwoked when lin/log radio button or any variable in the list is changed"""
        # find out which y-scale should be used
        scale = self.radio(self.fLin_2, self.fLog_2)
        if scale == 'log':loga = True
        else: loga = False
        # find out which variable should be plotted
        try:
            comp = self.availableVars.currentItem().text()
        except:
            return
        # Exctract that variable from netCDF-dataset and save to Y
        self.plotTitle = self.plotTitle[:self.plotTitle.rfind(':')+2] + comp+' ['+units.get(comp,units['REST'])[0]+']'
        if not self.sumSelection.isChecked():
            Y = self.ncs.variables[comp][:]
        else:
            if self.availableVars.selectedItems() != []:
                Y = sum(self.ncs.variables[c.text()][:] for c in self.availableVars.selectedItems())
                self.plotTitle = self.plotTitle[:self.plotTitle.rfind(':')+2] + self.availableVars.selectedItems()[0].text()+' etc.'
            else:
                Y = self.ncs.variables[comp][:]

        # Are the values non-negative?
        positive = all(Y>=0)
        # Are all the values zeros?
        zeros = all(Y==0)
        # if non-negative and not all zeros, and log-scale is possible without any fixes
        if not zeros and not positive and scale=='log':
            Y[Y<=0] = 1e-20
            # Mark plot with red to warn that negative values have been deleted
            self.outplot.setPen({'color':'r','width': 2.0})
            positive = True
        else:
            self.outplot.setPen({'color':'b','width': 2.0})

        if scale=='log' and positive and not zeros:
            # To avoid very small exponentials, zeros are changed to nearest small number
            if not all(Y>0):
                uniqs = unique(Y)
                if len(uniqs)>1:
                    Y[Y==0] = uniqs[1]
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
        time = self.ncs.variables[self.hnames[0]][:]/3600
        if ma.is_masked(time):
            self.outplot.setData(time[~time.mask],Y[~time.mask])
        else:
            self.outplot.setData(time,Y)
        self.plotResultWindow.setLogMode(y=loga)
        # update title
        self.plotResultTitle.setText(self.plotTitle)


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
        try: self.ncs.close()
        except: pass
        try:
            while True:
                try: self.availableVars.itemSelectionChanged.disconnect(self.showOutputUpdate)
                except TypeError: break
            self.availableVars.clear()
            self.plotResultWindow.clear()
        except:
            pass


    def closenetcdf_mass(self):
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
                self.popup('Compiled', 'Compiled ARCA succesfully', icon=0)
            else:
                self.popup('Compile error', 'Error in compiling, see output from terminal', icon=2)
            if self.ReplChem.isChecked(): self.ReplChem.setChecked(False)
            if self.makeClean.isChecked(): self.makeClean.setChecked(False)


    def remake(self):
        if self.running != None:
            if self.ReplChem.isChecked():
                self.editMakefile(mod=self.chemistryModules.currentText())
                if self.makeClean.isChecked():
                    target = "clean_current_chemistry"
                else:
                    target = "clean"
                self.compile = Popen(["make", target])#, stdout=subprocess.PIPE,stderr=subprocess.STDOUT,stdin=None)
                while True:
                    self.running = self.compile.poll()
                    if self.running != None: break
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
        if text != '':
            if text[-1] == '.':
                text = text[:-1]
                strict = True
        self.availableVars.clear()
        if text == '':
            self.availableVars.addItems(self.hnames)
        else:
            for c in self.hnames:
                if strict:
                    if text == c.upper():
                        self.availableVars.addItem(c)
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
            self.namesdat.scrollToItem(self.namesdat.item(0), QtGui.QAbstractItemView.PositionAtTop)
        else:
            for c in NAMES:
                if strict:
                    if text == c.upper():
                        if not namesPyInds[c] == divider_i:
                            self.namesdat.item(namesPyInds[c]).setSelected(True)
                else:
                    if text in c.upper():
                        if not namesPyInds[c] == divider_i:
                            item = self.namesdat.item(namesPyInds[c])
                            item.setSelected(True)
                            self.namesdat.scrollToItem(item, QtGui.QAbstractItemView.PositionAtTop)


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

CustomCommandsCheatSheet = """
Options available in CUSTOM_NML (and defaults)
LOG: USE_RAOULT=T
LOG: VARIABLE_DENSITY=F
REA: DMPS_TRES_MIN=  10.0
REA: START_TIME_S=  0.0
REA: DMPS_MULTI=  1000000.0
CHR: INITIALIZE_WITH=""
INT: INITIALIZE_FROM=0
REA: VP_MULTI=  1.0000000000000000,
LOG: DONT_SAVE_CONDENSIBLES=F
INT: LIMIT_VAPOURS=999999
REA: END_DMPS_SPECIAL=1.0000000000000000E+100
LOG: NO2_IS_NOX=F
LOG: NO_NEGATIVE_CONCENTRATIONS=T
REA: FLOAT_CHEMISTRY_AFTER_HRS=1.0000000000000000E+100
LOG: USE_RH_CORRECTION=T
LOG: TEMP_DEP_SURFACE_TENSION=F
LOG: USE_DIFF_DIA_FROM_DIFF_VOL=F
REA: SPEED_DT_LIMIT=300d0, 300d0, 300d0
LOG: ENABLE_END_FROM_OUTSIDE=T
REA: LIMIT_FOR_EVAPORATION=0.0
REA: MIN_CONCTOT_CC_FOR_DVAP=1d3
"""


if __name__ == '__main__':
    print(CurrentVersion+' started at:', ( time.strftime("%B %d %Y, %H:%M:%S", time.localtime())))
    app = QtWidgets.QApplication([])
    qt_box = QtBoxGui()
    qt_box.show()
    styles = QtWidgets.QStyleFactory.keys()
    if "Fusion" in styles:
        app.setStyle(QtWidgets.QStyleFactory.create("Fusion"))
    else:
        print('Available styles: ',styles)
    app.exec_()
