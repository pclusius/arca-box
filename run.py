from PyQt5 import QtCore, QtWidgets, QtGui, uic
import pyqtgraph as pg
from gui import vars, gui5
import subprocess
from numpy import linspace,log10,sqrt,exp,pi,sin
import numpy as np
try:
    import netCDF4
    netcdf = True
except:
    netcdf = False

## Some constants --------------------------------------------
column_widths = [140,70,70,70,70,90,50,3]

# available units for variables
# used to fill the tables and graphs with appropriate units
units = {
'TEMPK': ['K','C'],
'PRESSURE': ['Pa','hPa','bar','kPa', 'mbar'],
'REL_HUMIDITY': ['%'],
'CONDENS_SINK':['1/s'],
'CON_SIN_NITR':['1/s'],
'SW_RADIATION':['W/m2'],
'ION_PROD_RATE':['ip/cm3 s'],
'NUC_RATE_IN':['1/cm3 s'],
'REST':['#/cm3','ppm','ppb','ppt','ppq']
}

## Read model names--------------------------------------------
path_to_names = 'src/NAMES.dat'
f = open(path_to_names)
NAMES = []
namesPyInds = {}
namesFoInds = {}
default_path = 'gui/defaults'
tempfile = 'input/tmp_b65d729f784bc8fcfb4beb009ac7a31d'
## -----------------------------------------------------------

## Create lists and dictionaries related to NAMES-------------
i = 0;k=1
for line in f:
    name = line[:-1]
    if '#' in line:
        NAMES.append('MCM compounds start here')
        divider_i=i
    else:
        NAMES.append(name)
    namesPyInds[name] = i
    i = i+1
    namesFoInds[name] = k
    k = k+1
f.close()
## -----------------------------------------------------------
nml = vars.INITFILE(NAMES)

class Comp:
    def __init__(self):
        self.index  = 0
        self.mode  = 0
        self.col  = -1
        self.multi = 1e0    # Multiplication factor in MODE0
        self.shift = 0e0    # Constant to be addded in MODE0
        self.min = 1e1      # Minimum value for the parametrized concentration OR constant value if max <= min
        self.max = 1e5      # Peak value
        self.sig = 2.34e0      # Standard deviation for the Gaussian=sig of the bell curve
        self.mju = 12e0     # Time of peak value
        self.fv  = 0e0      # Angular frequency [hours] of modifying sine function
        self.ph  = 0e0      # Angular frequency [hours] of modifying sine function
        self.am  = 1e0      # Amplitude of modificaion
        self.name = 'NONAME'# Human readable name for modified variable
        self.unit = '#/cm3'     # unit name
        self.Find = 1
        self.pmInUse = 'No'
        self.sl_1 = 39
        self.sl_2 = 84
        self.sl_3 = 0
        self.sl_4 = 0
        self.sl_5 = 20
        self.gain = 50

class QtBoxGui(gui5.Ui_MainWindow,QtWidgets.QMainWindow):
    """Main program window."""
    def __init__(self):
        super(QtBoxGui,self).__init__()
        self.setupUi(self)

    # -----------------------
    # Common stuff
    # -----------------------
        self.setWindowTitle('Superbox utility')
        self.prints = 1
        self.plots  = 0
        self.show_extra_plots = ''
        self.printButton.clicked.connect(lambda: self.print_values())
        self.saveButton.clicked.connect(lambda: self.save_file())
        self.loadButton.clicked.connect(lambda: self.browse_path(None, 'load'))
        self.actionSave_2.triggered.connect(self.save_file)
        self.actionPrint.triggered.connect(lambda: self.print_values())
        self.actionOpen.triggered.connect(lambda: self.browse_path(None, 'load'))
        self.actionQuit_Ctrl_Q.triggered.connect(self.close)
        self.saveDefaults.clicked.connect(lambda: self.save_file(file=default_path))

    # -----------------------
    # tab General options
    # -----------------------
        self.namesdat.clear()
        self.namesdat.addItems(NAMES)
        item = self.namesdat.item(divider_i)
        item.setFlags(item.flags() & ~QtCore.Qt.ItemIsEnabled & ~QtCore.Qt.ItemIsSelectable)

        self.runtime.valueChanged.connect(lambda: self.updteGraph())

        # Prepare the variable table
        for i in range(len(column_widths)):
            self.selected_vars.setColumnWidth(i, column_widths[i])

        # add minimum requirements
        self.add_new_line('TEMPK', 0)
        self.add_new_line('PRESSURE', 1)

        self.browseCommonIn.clicked.connect(lambda: self.browse_path(self.case_dir, 'dir'))
        self.browseCommonOut.clicked.connect(lambda: self.browse_path(self.lineEdit_13, 'dir'))
        self.browseEnv.clicked.connect(lambda: self.browse_path(self.env_file, 'file'))
        self.browseMcm.clicked.connect(lambda: self.browse_path(self.mcm_file, 'file'))
        self.browsePar.clicked.connect(lambda: self.browse_path(self.dmps_file, 'file'))
        self.browseXtr.clicked.connect(lambda: self.browse_path(self.extra_particles, 'file'))
        self.checkBox_aer.stateChanged.connect(lambda: self.grayIfNotChecked(self.checkBox_aer,self.groupBox_8))
        self.fsave_division.valueChanged.connect(self.toggle_printtime)

    # -----------------------
    # tab Input variables
    # -----------------------
        self.butMoveToSelVars.clicked.connect(self.select_compounds)
        self.markAll.clicked.connect(lambda: self.markReverseSelection('all'))
        self.invertMarks.clicked.connect(lambda: self.markReverseSelection('inv'))
        self.butRemoveSelVars.clicked.connect(self.remv_item)
        self.selected_vars.setColumnHidden(7, True)
        self.selected_vars.verticalHeader().setVisible(False);
        self.loadFixed.clicked.connect(lambda: self.browse_path(None, 'fixed', ftype="KPP def (*.def)"))

    # -----------------------
    # tab Function creator
    # -----------------------
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
        self.resetSlider(self.fWidth, default.sl_1)
        self.resetSlider(self.fPeak,  default.sl_2)
        self.resetSlider(self.fFreq,  default.sl_3)
        self.resetSlider(self.fPhase,default.sl_4)
        self.resetSlider(self.fAmp,  default.sl_5)
        self.resetSlider(self.gain,   default.gain)

        self.resW.clicked.connect(lambda: self.resetSlider(self.fWidth, default.sl_1))
        self.resP.clicked.connect(lambda: self.resetSlider(self.fPeak,  default.sl_2))
        self.resA.clicked.connect(lambda: self.resetSlider(self.fFreq,  default.sl_3))
        self.resPh.clicked.connect(lambda: self.resetSlider(self.fPhase,default.sl_4))
        self.resAm.clicked.connect(lambda: self.resetSlider(self.fAmp,  default.sl_5))
        self.resG.clicked.connect(lambda: self.resetSlider(self.gain,   default.gain))
        self.fPeak.valueChanged.connect(lambda: self.updteGraph())
        self.fFreq.valueChanged.connect(lambda: self.updteGraph())
        self.fPhase.valueChanged.connect(lambda: self.updteGraph())
        self.fAmp.valueChanged.connect(lambda: self.updteGraph())
        self.gain.valueChanged.connect(lambda: self.updteGraph())
        self.PLOT.showGrid(x=True,y=True)
        self.PLOT.showButtons()
        self.updteGraph()
    # -----------------------
    # tab Advanced
    # -----------------------
        self.frameBase.setEnabled(False)
        self.butVapourNames.clicked.connect(lambda: self.browse_path(self.vap_names, 'file'))
        self.butVapour.clicked.connect(lambda: self.browse_path(self.vap_props, 'file'))
        self.vap_logical.stateChanged.connect(lambda: self.grayIfNotChecked(self.vap_logical,self.frameVap))
        self.resolve_base.stateChanged.connect(lambda: self.grayIfNotChecked(self.resolve_base,self.frameBase))
        self.use_dmps_special.stateChanged.connect(lambda: self.toggle_gray(self.use_dmps_special,self.gridLayout_11))
    # -----------------------
    # tab Process Monitor
    # -----------------------
        fixedFont = QtGui.QFontDatabase.systemFont(QtGui.QFontDatabase.FixedFont)
        fixedFont.setPointSize(10)
        self.MonitorWindow.setFont(fixedFont)
        self.frameStop.setEnabled(False)
        self.startButton.clicked.connect(self.startBox)
        self.stopButton.clicked.connect(self.stopBox)
        self.boxProcess = 0 # Superbox run handle
        # self.r = 0 # Superbox output file handle
        self.Timer = QtCore.QTimer(self);
        self.pollTimer = QtCore.QTimer(self);
        self.Timer.timeout.connect(self.updateOutput)
        self.pollTimer.timeout.connect(self.pollMonitor)
        self.pauseScroll.clicked.connect(self.checkPauseScroll)
        self.pauseScrolling = False

    # -----------------------
    # tab Plot Monitor
    # -----------------------

        if netcdf:
            self.show_netcdf.hide()
            self.loadNetcdf.setEnabled(True)
            self.fLog_2.clicked.connect(self.showOutputUpdate)
            self.fLin_2.clicked.connect(self.showOutputUpdate)

        else:
            self.show_netcdf.show()
            self.loadNetcdf.setEnabled(False)
            self.fLin_2.setEnabled(False)
            self.fLog_2.setEnabled(False)

        self.plotResultWindow.showGrid(x=True,y=True)
        self.plotResultWindow.setBackground('w')
        pen = pg.mkPen(color=(0,0,0), width=1)
        self.plotResultWindow.getAxis('left').setPen(pen)
        self.plotResultWindow.getAxis('bottom').setPen(pen)
        self.loadNetcdf.clicked.connect(lambda: self.browse_path(None, 'plot', ftype="NetCDF (*.nc)"))
        try: self.load_initfile(default_path)
        except: self.save_file(file=default_path)


    # -----------------------
    # Class methods
    # -----------------------
    def checkPauseScroll(self):
        if self.pauseScroll.isChecked() == True:
            self.pauseScrolling = True
        else:
            self.pauseScrolling = False

    def resetSlider(self, slider, pos):
        slider.setProperty("value", pos)

    def scales(self):
        rt = self.runtime.value()
        # rt=24
        scf = 24
        wScale = scf/2/200.0
        pScale = scf*1.1905/200.0
        aScale = 0.1
        phScale = scf/0.2/200.0
        ampScale = 1/20.0
        return wScale,pScale,aScale,phScale,ampScale,rt

    def updteGraph(self, label='Current edit'):
        gain = 10**(self.gain.value()/50.-1)
        wScale,pScale,aScale,phScale,ampScale,rt = self.scales()
        wScale   = wScale   * gain
        pScale   = pScale   * gain
        aScale   = aScale   * gain
        phScale  = phScale  * gain
        ampScale = ampScale * gain

        x = linspace(0,rt,200)
        yscale = self.radio(self.fLin, self.fLog)

        dummy.sig = self.fWidth.value()*wScale
        if abs(dummy.sig)<0.01: dummy.sig = 0.01
        try:
            dummy.min = float(self.fMin.text())
        except:
            dummy.min = 0
        try:
            dummy.max = float(self.fMax.text())
        except:
            dummy.max=0
        dummy.gain = self.gain.value()
        dummy.mju = self.fPeak.value()*pScale
        dummy.fv = self.fFreq.value()*aScale
        dummy.ph = self.fPhase.value()*phScale
        dummy.am = self.fAmp.value()*ampScale
        dummy.sl_1 = self.fWidth.value()
        dummy.sl_2 = self.fPeak.value()
        dummy.sl_3 = self.fFreq.value()
        dummy.sl_4 = self.fPhase.value()
        dummy.sl_5 = self.fAmp.value()

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
        try: # delete old legend if it exists:
            if self.plots > 0:
                self.legend.scene().removeItem(self.legend)
            self.plots = 1
        except Exception as e:
            print(e)

        self.legend = self.PLOT.addLegend()
        self.PLOT.plot(x,norm,pen=pg.mkPen('r', width=4), clear=True, name=label)

        if self.show_extra_plots != '':
            y=self.gauss(vars.mods[self.show_extra_plots], yscale,rt)
            self.PLOT.plot(x,y,pen=pg.mkPen(color='k', width=3,style=QtCore.Qt.DotLine), name=self.show_extra_plots)

    def saveParamValues(self):
        compound = self.names_sel_2.selectedItems()
        if compound != []:
            target = compound[0].text()
            vars.mods[target].min  = dummy.min
            vars.mods[target].max  = dummy.max
            vars.mods[target].sig  = dummy.sig
            vars.mods[target].mju  = dummy.mju
            vars.mods[target].fv   = dummy.fv
            vars.mods[target].ph   = dummy.ph
            vars.mods[target].am   = dummy.am
            vars.mods[target].mode = dummy.mode
            vars.mods[target].gain = dummy.gain
            vars.mods[target].sl_1 = dummy.sl_1
            vars.mods[target].sl_2 = dummy.sl_2
            vars.mods[target].sl_3 = dummy.sl_3
            vars.mods[target].sl_4 = dummy.sl_4
            vars.mods[target].sl_5 = dummy.sl_5
            vars.mods[target].pmInUse = 'Yes'
            confirmText = 'Saved to ' + target
            self.confirm.setText(confirmText)
            self.confirm.show()
            QtCore.QTimer.singleShot(2000, lambda : self.confirm.hide())
            self.updteGraph(label=target)
            for i in range(self.selected_vars.rowCount()):
                if self.selected_vars.item(i,0).text() == target:
                    self.selected_vars.cellWidget(i, 4).setCurrentIndex(1)
                    break

    def loadParamValues(self):
        compound = self.names_sel_2.selectedItems()
        if compound != []:
            target = compound[0].text()
            self.fMin.setText(str(vars.mods[target].min))
            self.fMax.setText(str(vars.mods[target].max))
            self.gain.setProperty("value", vars.mods[target].gain)
            self.fWidth.setProperty("value", vars.mods[target].sl_1)
            self.fPeak.setProperty("value", vars.mods[target].sl_2)
            self.fFreq.setProperty("value", vars.mods[target].sl_3)
            self.fPhase.setProperty("value", vars.mods[target].sl_4)
            self.fAmp.setProperty("value", vars.mods[target].sl_5)
            if vars.mods[target].mode < 2:
                self.fLin.setChecked(True)
                self.fLog.setChecked(False)
            else:
                self.fLin.setChecked(False)
                self.fLog.setChecked(True)
            self.updteGraph(label = target)
            confirmText = 'Loaded parameters from ' + target
            self.confirm.setText(confirmText)
            self.confirm.show()
            QtCore.QTimer.singleShot(2000, lambda : self.confirm.hide())


    def gauss(self,comp,ysc,rt):
        x = linspace(0,rt,200)
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

    def radio(self,*buts):
        for but in buts:
            if but.isChecked():
                if but.text() == 'Linear':
                    self.PLOT.setLogMode(y=False)
                    return 'lin'
                if but.text() == 'Logarithmic':
                    self.PLOT.setLogMode(y=True)
                    return 'log'

    def browse_path(self, target, mode, ftype=None):
        """Browse for file or folder (depending on 'mode' and write the outcome to 'target')"""
        dialog = QtWidgets.QFileDialog()
        if ftype != None:
            dialog.setNameFilter(ftype)

        if mode == 'dir': dialog.setFileMode(QtWidgets.QFileDialog.Directory)
        if mode == 'dir': dialog.setOption(QtWidgets.QFileDialog.ShowDirsOnly)
        options = dialog.Options()
        options |= dialog.DontUseNativeDialog
        if mode == 'dir':
            path = dialog.getExistingDirectory(self, 'Choose Directory', options=options)
        # if mode == 'file' or mode == 'load' or mode == 'fixed': path = dialog.getOpenFileName(self, 'Choose File', options=options)[0]
        else:
            dialog.setOption(QtWidgets.QFileDialog.DontUseNativeDialog)
            dialog.setWindowTitle('Choose File')
            if dialog.exec() == 1:
                path = dialog.selectedFiles()[0]
            else: path=''
            # path = dialog.getOpenFileName(self, 'Choose File', options=options)[0]
        if path != '':
            if mode == 'load':
                self.load_initfile(path)
            elif mode == 'fixed':
                self.loadFixedFile(path)
            elif mode == 'plot':
                self.showOutput(path)
            else:
                target.clear()
                target.insert(path)

    def save_file(self, file=None):
        self.update_nml()
        if file==None:
            dialog = QtWidgets.QFileDialog()
            options = dialog.Options()
            options |= dialog.DontUseNativeDialog
            file = dialog.getSaveFileName(self, 'Save INITFILE', options=options)[0]
        if file != '':
            self.print_values(file)

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
        self.selected_vars.setSortingEnabled(False)
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
        self.selected_vars.setSortingEnabled(True)
        self.updateOtherTabs()

    def updateOtherTabs(self):
        """updates variable lists in other tabs"""
        self.names_sel.clear()
        self.names_sel_2.clear()
        for i in range(self.selected_vars.rowCount()):
            self.names_sel.addItem(self.selected_vars.item(i,0).text())
            self.names_sel_2.addItem(self.selected_vars.item(i,0).text())

    def loadFixedFile(self, path):
        f = open(path, 'r')
        indef = False
        count = 0
        for line in f:
            if '#DEFFIX' in line.upper():
                indef = True
                continue
            if indef and '#' in line and not '#DEFFIX' in line.upper():
                break
            if indef and '=' in line and '//' not in line:
                i = line.find('=')
                comp = line[:i].strip()
                if comp not in vars.mods:
                    self.namesdat.item(namesPyInds[comp]).setSelected(True)
                    count = count +1
        self.popup('File parsed', 'Selected %d variables'%count)

    def add_new_line(self, name, unit_ind, cols=[],createNew=True, unt=0):
        """adds items to variable table"""
        self.selected_vars.setSortingEnabled(False);
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
        pmInUse.setCurrentIndex(cols[4])
        self.selected_vars.horizontalHeader().setStretchLastSection(True)

        for i in range(4):
            self.selected_vars.setItem(row, i, QtWidgets.QTableWidgetItem(cols[i]))
        if name == 'PRESSURE' :
            self.selected_vars.setItem(row, 3, QtWidgets.QTableWidgetItem('1e5'))
        self.selected_vars.setCellWidget(row, i+1, pmInUse )
        self.selected_vars.setCellWidget(row, i+2, unit )
        self.selected_vars.setCellWidget(row, i+3, markBut )

        self.selected_vars.setItem(row, i+4, QtWidgets.QTableWidgetItem('%03d'%(namesFoInds[name])))

        self.selected_vars.sortItems(7, QtCore.Qt.AscendingOrder)
        self.selected_vars.setSortingEnabled(True)
        self.updateOtherTabs()
        if createNew:
            vars.mods[name] = Comp()
            vars.mods[name].Find = namesFoInds[name]
            vars.mods[name].name = name# Human readable name for modified variable

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

    def grayIfNotChecked(self, guard, frame):
        if guard.isChecked() == True:
            frame.setEnabled(True)
        else:
            frame.setEnabled(False)

    def toggle_printtime(self):
        if self.fsave_division.value() != 0:
            self.fsave_interval.setEnabled(False)
        else:
            self.fsave_interval.setEnabled(True)

    def print_values(self, file=None):
        self.update_nml()
        self.prints += 1
        if (file):
            f = open(file, 'w')
            f.write('#'+('-')*50+'\n')
            f.write('#      Superbox setting file: %s\n'%('file name'))
            f.write('#'+('-')*50+'\n\n')
            nml.printall(vars.mods, target='f',f=f)
            f.close()
            if file == default_path:
                self.popup('', 'Defaults saved')
            elif file != tempfile:
                self.popup('Saved settings to', file)
        else:
            print(('\n')*10, )
            print('#',('-')*50)
            print('#              Superbox setting file #%d'%(self.prints))
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
                self.updteGraph()
            else:
                self.plotTo.setChecked(False)
        if self.plotTo.isChecked() == False:
            self.show_extra_plots = ''
            self.updteGraph()

    def stopBox(self):
        self.Timer.stop()
        self.pollTimer.stop()
        self.boxProcess.kill()
        tout = self.boxProcess.wait(timeout=10)
        self.boxProcess.poll()
        self.toggle_frame(self.frameStop)
        self.toggle_frame(self.frameStart)

    def startBox(self):
        self.closenetcdf()

        currentWait = self.wait_for.value()
        if self.python.isChecked():
            currentPython = True
            self.python.setChecked(False)
        else:
            currentPython = False
        self.wait_for.setValue(-1)
        self.print_values(tempfile)
        self.python.setChecked(currentPython)
        self.wait_for.setValue(currentWait)

        try:
            self.boxProcess = subprocess.Popen(["./superbox.exe", " %s"%tempfile], stdout=subprocess.PIPE,stderr=subprocess.STDOUT,stdin=None)
            self.MonitorWindow.clear()
            self.Timer.start(10)
            self.pollTimer.start(1000)
            self.toggle_frame(self.frameStart)
            self.toggle_frame(self.frameStop)
        except:
            self.MonitorWindow.appendPlainText('\n    Could not start model executable, is it compiled?')


    def pollMonitor(self):
        status = self.boxProcess.poll()
        if status != None:
            self.stopBox()

    def updateOutput(self):
        fulltext = self.boxProcess.stdout.readline().decode("utf-8")
        self.MonitorWindow.insertPlainText(fulltext)
        if self.pauseScrolling == False:
            self.MonitorWindow.verticalScrollBar().setSliderPosition(self.MonitorWindow.verticalScrollBar().maximum());
        if 'SIMULATION HAS ENDED' in str(fulltext)[-50:]:
            self.MonitorWindow.setPlainText(self.MonitorWindow.toPlainText())
            self.MonitorWindow.verticalScrollBar().setSliderPosition(self.MonitorWindow.verticalScrollBar().maximum());
            self.stopBox()


    def checkbox(self, widget):
        if widget.isChecked() == True:
            return '.TRUE.'
        else:
            return '.FALSE.'

    def update_nml(self):
        # class _PATH:
        nml.PATH.CASE_DIR=self.case_dir.text()
        nml.PATH.CASE_NAME=self.case_name.text()
        nml.PATH.RUN_NAME = self.run_name.text()
        # class _FLAG:
        nml.FLAG.CHEMISTRY_FLAG=self.checkbox(self.checkBox_che)
        nml.FLAG.AEROSOL_FLAG=self.checkbox(self.checkBox_aer)
        nml.FLAG.ACDC_SOLVE_SS=self.checkbox(self.acdc_solve_ss)
        nml.FLAG.NUCLEATION='T'
        nml.FLAG.ACDC=self.checkbox(self.checkBox_acd)
        nml.FLAG.EXTRA_DATA='F'
        nml.FLAG.CURRENT_CASE='F'
        nml.FLAG.CONDENSATION=self.checkbox(self.condensation)
        nml.FLAG.COAGULATION=self.checkbox(self.coagulation)
        nml.FLAG.MODEL_H2SO4=self.checkbox(self.model_h2so4)
        nml.FLAG.RESOLVE_BASE=self.checkbox(self.resolve_base)
        # class _TIME:
        nml.TIME.RUNTIME=self.runtime.value()
        nml.TIME.DT=self.dt.value()
        nml.TIME.FSAVE_INTERVAL=self.fsave_interval.value()
        nml.TIME.PRINT_INTERVAL=self.print_interval.value()
        nml.TIME.FSAVE_DIVISION=self.fsave_division.value()
        nml.TIME.DATE=self.dateEdit.text()
        # class _PARTICLE:
        nml.PARTICLE.PSD_MODE=self.psd_mode.currentIndex()+1
        nml.PARTICLE.N_BINS_PARTICLE=self.n_bins_particle.value()
        nml.PARTICLE.MIN_PARTICLE_DIAM=self.min_particle_diam.text()
        nml.PARTICLE.MAX_PARTICLE_DIAM=self.max_particle_diam.text()
        # nml.PARTICLE.DMPS_DIR=self.dmps_dir.text()
        nml.PARTICLE.EXTRA_P_DIR=''
        nml.PARTICLE.DMPS_FILE=self.dmps_file.text()
        nml.PARTICLE.EXTRA_PARTICLES=self.extra_particles.text()
        nml.PARTICLE.DMPS_READ_IN_TIME=self.dmps_read_in_time.value()
        nml.PARTICLE.DMPS_HIGHBAND_LOWER_LIMIT=self.dmps_highband_lower_limit.text()
        nml.PARTICLE.DMPS_LOWBAND_UPPER_LIMIT=self.dmps_lowband_upper_limit.text()
        nml.PARTICLE.USE_DMPS=self.checkbox(self.use_dmps)
        nml.PARTICLE.USE_DMPS_SPECIAL=self.checkbox(self.use_dmps_special)
        # class _ENV:
        nml.ENV.ENV_PATH=self.case_dir.text()
        nml.ENV.ENV_FILE=self.env_file.text()
        # class _MCM:
        nml.MCM.MCM_PATH=self.case_dir.text()
        nml.MCM.MCM_FILE=self.mcm_file.text()
        # class _MISC:
        nml.MISC.LAT=self.lat.value()
        nml.MISC.LON=self.lon.value()
        nml.MISC.WAIT_FOR=self.wait_for.value()
        nml.MISC.PYTHON=self.checkbox(self.python)
        nml.MISC.DESCRIPTION=self.description.toPlainText()
        nml.MISC.CH_ALBEDO=self.ch_albedo.value()
        nml.MISC.DMA_F=self.dma_f.value()
        nml.MISC.RESOLVE_BASE_PRECISION=self.resolve_base_precision.value()
        nml.MISC.FILL_FORMATION_WITH=self.resolveHelper()
        # class _VAP:
        nml.VAP.VAP_LOGICAL=self.checkbox(self.vap_logical)
        nml.VAP.VAP_NAMES=self.vap_names.text()
        nml.VAP.VAP_PROPS=self.vap_props.text()

        for i in range(self.selected_vars.rowCount()):
            name = self.selected_vars.item(i,0).text()
            vars.mods[name].col = int(self.selected_vars.item(i,1).text())
            vars.mods[name].multi = float(self.selected_vars.item(i,2).text())
            vars.mods[name].shift = float(self.selected_vars.item(i,3).text())
            vars.mods[name].pmInUse = self.selected_vars.cellWidget(i,4).currentText()
            vars.mods[name].unit = self.selected_vars.cellWidget(i,5).currentText()
        nml.ENV.TEMPUNIT=vars.mods['TEMPK'].unit

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
        self.markReverseSelection('all')
        self.remv_item()
        if self.plotTo.isChecked() == True:
            self.plotTo.setChecked(False)
            self.show_extra_plots = ''
            self.updteGraph()


        def solve_for_parser(query):
            if query.upper() == 'NH3':
                return 1
            elif query.upper() == 'DMA':
                return 2
            else:
                return 0

        def parse_date(str):
            y = int(str[0:4] )
            m = int(str[5:7])
            d = int(str[8:10])
            return QtCore.QDate(y, m, d)


        f = open(file, 'r')

        for line in f:
            i = line.find('=')
            x = line.find('!')
            if i<x or (x==-1 and i!=-1):
                key = line[:i].upper().strip()
                strng = line[i+1:x].strip()
                if len(strng) > 0:
                    # remove end comma if exists
                    if strng[-1] == ",":
                        strng = strng[:-1]
                    # remove apostrophes if exists
                    if (strng[0] == "\'" and strng[-1] == "\'") or (strng[0] == "\"" and strng[-1] == "\""):
                        strng = strng[1:-1]
                    # change fortran boolean to python boolean
                    if strng == "T" or strng.upper() == ".TRUE." or strng == "F" or strng.upper() == ".FALSE.":
                        if strng == "T" or strng.upper() == ".TRUE.":
                            strng = True
                        elif strng == "F" or strng.upper() == ".FALSE.":
                            strng = False
                    # Check if srtng is a number
                    try:
                        float(strng)
                        isFl = True
                    except:
                        isFl = False

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
                            props = strng.split()
                            i = 0
                            n=len(props)
                            if i< n: vars.mods[name].mode   = int(props[i])  ;i=i+1
                            if i< n: vars.mods[name].col    = int(props[i])  ;i=i+1
                            if i< n: vars.mods[name].multi  = float(props[i].replace('d','e',1));i=i+1
                            if i< n: vars.mods[name].shift  = float(props[i].replace('d','e',1));i=i+1
                            if i< n: vars.mods[name].min    = float(props[i].replace('d','e',1));i=i+1
                            if i< n: vars.mods[name].max    = float(props[i].replace('d','e',1));i=i+1
                            if i< n: vars.mods[name].sig    = float(props[i].replace('d','e',1));i=i+1
                            if i< n: vars.mods[name].mju    = float(props[i].replace('d','e',1));i=i+1
                            if i< n: vars.mods[name].fv     = float(props[i].replace('d','e',1));i=i+1
                            if i< n: vars.mods[name].ph     = float(props[i].replace('d','e',1));i=i+1
                            if i< n: vars.mods[name].am     = float(props[i].replace('d','e',1));i=i+1
                            if i< n:
                                unt = props[i][1:-1]
                                if unt == '#': unt = '#/cm3'
                                vars.mods[name].unit = unt


            else:
                continue

            if   'WORK_DIR' == key: self.case_dir.setText(strng)
            elif 'CASE_DIR' == key: self.case_dir.setText(strng)
            elif 'CASE_NAME' == key: self.case_name.setText(strng)
            elif 'RUN_NAME' == key: self.run_name.setText(strng)
            elif 'CHEMISTRY_FLAG' == key: self.checkBox_che.setChecked(strng)
            elif 'AEROSOL_FLAG' == key: self.checkBox_aer.setChecked(strng)
            elif 'ACDC_SOLVE_SS' == key: self.acdc_solve_ss.setChecked(strng)
            # elif 'NUCLEATION' == key: self.      .setChecked(strng)
            elif 'ACDC' == key: self.checkBox_acd.setChecked(strng)
            # elif 'EXTRA_DATA' == key: self.      .setChecked(strng)
            # elif 'CURRENT_CASE' == key: self.      .setChecked(strng)
            elif 'CONDENSATION' == key: self.condensation.setChecked(strng)
            elif 'COAGULATION' == key: self.coagulation.setChecked(strng)
            elif 'MODEL_H2SO4' == key: self.model_h2so4.setChecked(strng)
            elif 'RESOLVE_BASE' == key: self.resolve_base.setChecked(strng)
            elif 'RUNTIME' == key and isFl: self.runtime.setValue(float(strng))
            # elif 'DT' == key: print(strng)#          -1,
            elif 'FSAVE_INTERVAL' == key and isFl: self.fsave_interval.setValue(float(strng))
            elif 'PRINT_INTERVAL' == key and isFl: self.print_interval.setValue(float(strng))
            elif 'FSAVE_DIVISION' == key and isFl: self.fsave_division.setValue(float(strng))
            elif 'DATE' == key: self.dateEdit.setDate(parse_date(strng))
            elif 'PSD_MODE' == key and isFl: self.psd_mode.setCurrentIndex(int(strng))
            elif 'N_BINS_PARTICLE' == key and isFl: self.n_bins_particle.setValue(int(strng))
            elif 'MIN_PARTICLE_DIAM' == key: self.min_particle_diam.setText(strng)#   1.0000000000000001E-009,
            elif 'MAX_PARTICLE_DIAM' == key: self.max_particle_diam.setText(strng)#   9.9999999999999995E-007,
            # elif 'DMPS_DIR' == key: print(strng)# "
            # elif 'EXTRA_P_DIR' == key: self.extra_p_dir.setText(strng)
            elif 'DMPS_FILE' == key: self.dmps_file.setText(strng)
            elif 'EXTRA_PARTICLES' == key: self.extra_particles.setText(strng)
            elif 'DMPS_READ_IN_TIME' == key and isFl: self.dmps_read_in_time.setValue(float(strng))
            elif 'DMPS_HIGHBAND_LOWER_LIMIT' == key: self.dmps_highband_lower_limit.setText(strng)
            elif 'DMPS_LOWBAND_UPPER_LIMIT' == key: self.dmps_lowband_upper_limit.setText(strng)
            elif 'USE_DMPS' == key: self.use_dmps.setChecked(strng)
            elif 'USE_DMPS_SPECIAL' == key: self.use_dmps_special.setChecked(strng)
            # elif 'ENV_PATH' == key: self.case_dir.setText(strng)
            elif 'ENV_FILE' == key: self.env_file.setText(strng)
            # elif 'TEMPUNIT' == key: print(strng)# "K
            # elif 'MCM_PATH' == key: print(strng)# "
            elif 'MCM_FILE' == key: self.mcm_file.setText(strng)# "
            # elif 'JD' == key: print(strng)#           2,
            elif 'LAT' == key and isFl: self.lat.setValue(float(strng))
            elif 'LON' == key and isFl: self.lon.setValue(float(strng))
            elif 'WAIT_FOR' == key and isFl: self.wait_for.setValue(int(strng))
            elif 'PYTHON' == key: self.python.setChecked(strng)
            elif 'DESCRIPTION' == key: self.description.setPlainText(strng)# "Just some keying
            # elif 'SOLVER' == key: print(strng)# "
            elif 'CH_ALBEDO' == key and isFl: self.ch_albedo.setValue(float(strng))#  0.20000000000000001     ,
            elif 'DMA_F' == key and isFl: self.dma_f.setValue(float(strng))
            elif 'RESOLVE_BASE_PRECISION' == key and isFl: self.resolve_base_precision.setValue(float(strng))
            elif 'FILL_FORMATION_WITH' == key: self.fill_formation_with.setCurrentIndex(solve_for_parser(strng))
            elif 'VAP_LOGICAL' == key: self.vap_logical.setChecked(strng)
            elif 'VAP_NAMES' == key: self.vap_names.setText(strng)
            elif 'VAP_PROPS' == key: self.vap_props.setText(strng)

        # manage different units
        for key in vars.mods:
            pm = 0
            if vars.mods[key].mode > 0:
                pm = 1
            unts = units.get(key,units['REST'])
            ui = 0
            for unit in unts:
                if unit.upper() == vars.mods[key].unit.upper():
                    break
                else:
                    ui=ui+1
            if ui==len(unts): ui=0

            if key == 'TEMPK' or key == 'PRESSURE':
                if key == 'TEMPK': row = 0
                if key == 'PRESSURE': row = 1
                self.selected_vars.setItem(row, 1, QtWidgets.QTableWidgetItem(str(vars.mods[key].col)))
                self.selected_vars.setItem(row, 2, QtWidgets.QTableWidgetItem(str(vars.mods[key].multi)))
                self.selected_vars.setItem(row, 3, QtWidgets.QTableWidgetItem(str(vars.mods[key].shift)))
                self.selected_vars.cellWidget(row,4).setCurrentIndex(pm)
                self.selected_vars.cellWidget(row,5).setCurrentIndex(ui)

            else:
                cols = [key,str(vars.mods[key].col),str(vars.mods[key].multi),str(vars.mods[key].shift),pm]
                self.add_new_line(key,2,cols=cols,createNew=False, unt=ui)

            wScale,pScale,aScale,phScale,ampScale,_ = self.scales()
            vars.mods[key].sl_1 = int(round(vars.mods[key].sig /wScale,0))
            vars.mods[key].sl_2 = int(round(vars.mods[key].mju /pScale,0))
            vars.mods[key].sl_3 = int(round(vars.mods[key].fv  /aScale,0))
            vars.mods[key].sl_4 = int(round(vars.mods[key].ph  /phScale,0))
            vars.mods[key].sl_5 = int(round(vars.mods[key].am  /ampScale,0))

    def showOutput(self, file):
        if netcdf:
            # First set plot mode to linear in order to avoid errors with zero values
            self.fLin_2.setChecked(True)
            self.fLog_2.setChecked(False)
            self.plotResultWindow.setLogMode(y=False)

            # Close all previus netcdf-files and clear plot
            self.closenetcdf()
            # Try to open netCDF-file
            try:
                self.ncs = netCDF4.Dataset(file, 'r')
            except:
                self.popup('Bummer...', 'Not a valid output file')
                return

            # collect all variables and dimensions from netCDF-dataset to hnames
            self.hnames = [i for i in self.ncs.variables]
            dims = [i for i in self.ncs.dimensions]

            # remove string variables (non-numbers and thus not plottable)
            cache = []
            for i,v in enumerate(self.hnames):
                try:
                    int(self.ncs.variables[v][0])
                    cache.append(v)
                except:
                    pass
            # remove constants
            self.hnames = []
            for i,v in enumerate(cache):
                if (len(np.shape(self.ncs.variables[v][:])) >1):
                    self.hnames.append(v)
                # Shifter ands multiplyers are left off for clarity
                if not 'Shifter' in v and not 'Multipl' in v:
                    self.hnames.append(v)
            # Now try plot first the third line, which is the first non-time variable
            try:
                self.outplot = self.plotResultWindow.plot(
                self.ncs.variables[self.hnames[0]][:]/3600,
                self.ncs.variables[self.hnames[2]][:],
                pen={'color':'b','width': 2.0}
                )
                self.availableVars.clear()
                self.availableVars.addItems(self.hnames)
                self.availableVars.item(2).setSelected(True)
                self.availableVars.itemSelectionChanged.connect(self.showOutputUpdate)

            # If fails, give information and return
            except:
                self.popup('Bummer...', "Output file does not contain any plottable data")
                return

            # All's well, decorate plot
            self.plotResultWindow.setLabel('bottom', 'Time', units='h')
            self.plotResultWindow.setTitle(self.hnames[2]+' '+units.get(self.hnames[2],units['REST'])[0])

        # if netCDF is not installed, return and do nothing
        else:
            return

    # This function is inwoked when lin/log radio button or any variable in the list is changed:
    def showOutputUpdate(self):
        # find out which y-scale should be used
        scale = self.radio(self.fLin_2, self.fLog_2)
        # find out which variable should be plotted
        i = self.availableVars.currentRow()
        # Exctract that variable from netCDF-dataset and save to Y
        Y = self.ncs.variables[self.hnames[i]][:]
        # Are the values non-negative?
        positive = all(Y>=0)
        # Are all the values zeros?
        zeros = all(Y==0)
        # if non-negative and not all zeros, and log-scale is possible
        if scale=='log' and positive and not zeros:
            # To avoid very small exponentials, zeros are changed to nearest small number
            if not all(Y>0):
                uniqs = np.unique(Y)
                if len(uniqs)>1:
                    Y[Y==0] = uniqs[1]
                # if the above fails, use linear scale
                else:
                    self.plotResultWindow.setLogMode(y=False)
            # if everything ok, use log scale
            self.plotResultWindow.setLogMode(y=True)
        # if linear scale was chosen, or negative values, or all zeros, use lin-scale
        else:
            self.plotResultWindow.setLogMode(y=False)
        # update data
        self.outplot.setData(self.ncs.variables[self.hnames[0]][:]/3600,Y)
        # update title
        self.plotResultWindow.setTitle(self.hnames[i]+' '+units.get(self.hnames[i],units['REST'])[0])


    def popup(self,title,message):
        msg = QtWidgets.QMessageBox()
        msg.setIcon(QtWidgets.QMessageBox.Warning)
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


dummy = Comp()
default = Comp()

if __name__ == '__main__':
    app = QtWidgets.QApplication([])
    qt_box = QtBoxGui()
    qt_box.show()
    app.exec_()
    # sys.exit(app.exec_())
