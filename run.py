from PyQt5 import QtCore, QtWidgets, QtGui, uic
import pyqtgraph as pg
from gui import vars, gui5
import subprocess
from numpy import linspace,log10,sqrt,exp,pi,sin

## Some constants --------------------------------------------
column_widths = [140,70,70,70,70,80,60,3]
all_units = [['C','K'],['Pa','hPa','bar','kPa', 'mbar'],['#','ppm','ppb','ppt','ppq'], ['as is']]
## Read model names--------------------------------------------
path_to_names = 'src/NAMES.dat'
f = open(path_to_names)
NAMES = []
namesPyInds = {}
namesFoInds = {}
## -----------------------------------------------------------

## Create lists and dictionaries related to NAMES-------------
i = 0;k=1
for line in f:
    name = line[:-1]
    if not '#' in line:
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
        self.unit = '#'     # unit name
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
        self.printButton.clicked.connect(self.print_values)
        self.saveButton.clicked.connect(self.save_file)
        self.loadButton.clicked.connect(self.load_file)
        self.actionSave_2.triggered.connect(self.save_file)
        self.actionPrint.triggered.connect(self.print_values)
        self.actionOpen.triggered.connect(self.browse_file)
        self.actionQuit_Ctrl_Q.triggered.connect(self.close)

    # -----------------------
    # tab General options
    # -----------------------
        self.namesdat.clear()
        self.namesdat.addItems(NAMES)
        self.runtime.valueChanged.connect(self.updteGraph)

        # Prepare the variable table
        for i in range(len(column_widths)):
            self.selected_vars.setColumnWidth(i, column_widths[i])

        # add minimum requirements
        self.add_new_line('TEMPK', 0)
        self.add_new_line('PRESSURE', 1)

        self.browseCommonIn.clicked.connect(lambda: self.browse_folder(self.case_dir))
        self.browseCommonOut.clicked.connect(lambda: self.browse_folder(self.lineEdit_13))
        self.browseEnv.clicked.connect(lambda: self.browse_file(self.env_file))
        self.browseMcm.clicked.connect(lambda: self.browse_file(self.mcm_file))
        self.browsePar.clicked.connect(lambda: self.browse_file(self.dmps_file))
        self.browseXtr.clicked.connect(lambda: self.browse_file(self.extra_particles))
        self.checkBox_aer.clicked.connect(lambda: self.grayIfNotChecked(self.checkBox_aer,self.groupBox_8))
        self.fsave_division.valueChanged.connect(self.toggle_printtime)

    # -----------------------
    # tab Input variables
    # -----------------------
        self.butMoveToSelVars.clicked.connect(self.select_compounds)
        self.markAll.clicked.connect(lambda: self.markReverseSelection('all'))
        self.invertMarks.clicked.connect(lambda: self.markReverseSelection('inv'))
        self.butRemoveSelVars.clicked.connect(self.remv_item)
        self.selected_vars.setColumnHidden(7, True)

    # -----------------------
    # tab Function creator
    # -----------------------
        self.confirm.hide()
        self.plotTo.clicked.connect(self.select_compounds_for_plot)
        self.saveParams.clicked.connect(self.saveParamValues)
        self.loadParams.clicked.connect(self.loadParamValues)
        self.fMin.editingFinished.connect(self.updteGraph)
        self.fMax.editingFinished.connect(self.updteGraph)
        self.fLog.clicked.connect(self.updteGraph)
        self.fLin.clicked.connect(self.updteGraph)
        self.fWidth.valueChanged.connect(self.updteGraph)
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
        self.fPeak.valueChanged.connect(self.updteGraph)
        self.fFreq.valueChanged.connect(self.updteGraph)
        self.fPhase.valueChanged.connect(self.updteGraph)
        self.fAmp.valueChanged.connect(self.updteGraph)
        self.gain.valueChanged.connect(self.updteGraph)
        self.PLOT.showGrid(x=True,y=True)
        self.PLOT.showButtons()
        self.updteGraph()
    # -----------------------
    # tab Advanced
    # -----------------------
        self.frameBase.setEnabled(False)
        self.butVapourNames.clicked.connect(lambda: self.browse_file(self.vap_names))
        self.butVapour.clicked.connect(lambda: self.browse_file(self.vap_props))
        self.vap_logical.clicked.connect(lambda: self.grayIfNotChecked(self.vap_logical,self.frameVap))
        self.resolve_base.clicked.connect(lambda: self.grayIfNotChecked(self.resolve_base,self.frameBase))
        self.use_dmps_special.clicked.connect(lambda: self.toggle_gray(self.use_dmps_special,self.gridLayout_11))
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
    # Class methods
    # -----------------------
    def checkPauseScroll(self):
        if self.pauseScroll.isChecked() == True:
            self.pauseScrolling = True
        else:
            self.pauseScrolling = False

    def resetSlider(self, slider, pos):
        slider.setProperty("value", pos)

    def updteGraph(self):
        gain = 10**(self.gain.value()/50.-1)
        rt = self.runtime.value()
        wScale = rt/2/200.0 * gain
        pScale = rt*1.1905/200.0 * gain
        aScale = 0.1 * gain
        phScale = rt/0.2/200.0 * gain
        ampScale = 1/20.0 * gain

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
        self.PLOT.plot(x,norm,pen=pg.mkPen('y', width=2), clear=True, name='Current edit')

        if self.show_extra_plots != '':
            y=self.gauss(vars.mods[self.show_extra_plots], yscale,rt)
            self.PLOT.plot(x,y,pen=pg.mkPen(color=(200, 200, 255), width=1,style=QtCore.Qt.DotLine), name=self.show_extra_plots)

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
            self.updteGraph()
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
            self.updteGraph()
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

    def load_file(self):
        print('loads file')
        pass

    def browse_folder(self, target):
        dialog = QtWidgets.QFileDialog()
        dialog.setFileMode(QtWidgets.QFileDialog.Directory)
        dialog.setOption(QtWidgets.QFileDialog.ShowDirsOnly)
        options = dialog.Options()
        options |= dialog.DontUseNativeDialog
        directory = dialog.getExistingDirectory(self, 'Choose Directory', options=options)
        if directory != '':
            target.clear()
            target.insert(directory)

    def browse_file(self, target):
        dialog = QtWidgets.QFileDialog()
        options = dialog.Options()
        options |= dialog.DontUseNativeDialog
        file = dialog.getOpenFileName(self, 'Choose File', options=options)[0]
        if file != '':
            target.clear()
            target.insert(file)

    def save_file(self):
        self.update_nml()
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


    def add_new_line(self, name, unit_ind):
        """adds items to variable table"""
        self.selected_vars.setSortingEnabled(False);
        row = self.selected_vars.rowCount()
        self.selected_vars.insertRow(row)
        item = self.namesdat.item(namesPyInds[name])
        item.setFlags(item.flags() & ~QtCore.Qt.ItemIsEnabled & ~QtCore.Qt.ItemIsSelectable)

        text = '%s'%(name)
        pmInUse = QtWidgets.QComboBox()
        pmInUse.addItems(['No','Yes'])
        unit = QtWidgets.QComboBox()
        unit.addItems(all_units[unit_ind])
        markBut = QtWidgets.QPushButton()
        # markBut.setFixedSize(column_widths[6],30)
        markBut.setCheckable(True)
        if name == 'TEMPK' or name == 'PRESSURE' :
            markBut.setEnabled(False)
            markBut.setToolTip("Temperature and pressure are always required and cannot be marked for removal")
        else:
            markBut.setToolTip("Mark variable for removal")

        markBut.setText('mark')
        cols = [text, '-1','1.0', '0.0','no']
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
        self.boxProcess.kill()
        tout = self.boxProcess.wait(timeout=10)
        self.boxProcess.poll()
        self.toggle_frame(self.frameStop)
        self.toggle_frame(self.frameStart)

    def startBox(self):
        currentWait = self.wait_for.value()
        if self.python.isChecked():
            currentPython = True
            self.python.setChecked(False)
        else:
            currentPython = False
        self.wait_for.setValue(-1)
        self.print_values('input/tmp_b65d729f784bc8fcfb4beb009ac7a31d')
        self.python.setChecked(currentPython)
        self.wait_for.setValue(currentWait)
        try:
            self.boxProcess = subprocess.Popen(["./superbox.exe", " input/tmp_b65d729f784bc8fcfb4beb009ac7a31d"], stdout=subprocess.PIPE,stderr=subprocess.STDOUT,stdin=None)
            self.MonitorWindow.clear()
            self.Timer.start(10)
            self.pollTimer.start(150000)
            self.toggle_frame(self.frameStart)
            self.toggle_frame(self.frameStop)
        except:
            self.MonitorWindow.appendPlainText('\n    Could not start model executable, is it compiled?')

    def pollMonitor(self):
        status = self.boxProcess.poll()
        if status != None:
            self.pollTimer.stop()
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
        nml.ENV.TEMPUNIT=vars.mods['TEMPK'].unit
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

        return

    def resolveHelper(self):
        text = self.fill_formation_with.currentText()
        if text == 'Fixed ratio':
            return ''
        elif 'NH3' in text:
            return 'NH3'
        elif 'DMA' in text:
            return 'DMA'


dummy = Comp()
default = Comp()

if __name__ == '__main__':
    app = QtWidgets.QApplication([])
    qt_box = QtBoxGui()
    qt_box.show()
    app.exec_()
    # sys.exit(app.exec_())
