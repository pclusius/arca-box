# try:
#     from PySide2 import QtCore, QtWidgets, QtGui
#     #from PySide2.QtCore import QProcess
# except:
#     try:
#         from PyQt5 import QtCore, QtWidgets
#     except:
#         print("You need PySide2 or PyQt5 to run this program")
#
from PyQt5 import QtCore, QtWidgets, QtGui, uic
import pyqtgraph as pg
from gui import vars, gui5
import subprocess
from numpy import linspace,log10,sqrt,exp,pi,sin

#
# import matplotlib
# # Make sure that we are using QT5
# matplotlib.use('Qt5Agg')
# from numpy import linspace, sin
# from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
# from matplotlib.figure import Figure

# For guitest.py:
# try:
#     from PySide2 import QtCore, QtGui, QtWidgets
# except:
#     try:
#         from PyQt5 import QtCore, QtGui, QtWidgets
#     except:
#         print("You need PySide2 or PyQt5 to run this program")
# #

## All variables stored here ---------------------------------

## Some constants --------------------------------------------
# X = linspace(0,3,100)
# Y = sin(X)
# ## -----------------------------------------------------------

## Some constants --------------------------------------------
column_widths = [140,70,70,70,70,50,60, 10]
all_units = [['K','C'],['Pa','hPa','bar','kPa', 'mbar'],['#','ppm','ppb','ppt','ppq'], ['as is']]
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


# class MyMplCanvas(FigureCanvas):
#     """Ultimately, this is a QWidget (as well as a FigureCanvasAgg, etc.)."""
#
#     def __init__(self, parent, width, height, dpi):
#         fig = Figure(figsize=(width, height), dpi=dpi)
#         self.axes = fig.add_subplot(111)
#         # self.axes.set_yscale('log')
#         self.axes.set_position([0.1, 0.1, 0.88, 0.88])
#         self.axes.grid()
#         self.compute_initial_figure()
#
#         FigureCanvas.__init__(self, fig)
#         self.setParent(parent)
#         FigureCanvas.setSizePolicy(self,
#                QtWidgets.QSizePolicy.Expanding,
#                QtWidgets.QSizePolicy.Expanding)
#         FigureCanvas.updateGeometry(self)
#
#     # def compute_initial_figure(self):
#     #     pass
#
#
#
# class MyStaticMplCanvas(MyMplCanvas):
#     """Simple canvas with a sine plot."""
#     def compute_initial_figure(self):
#         self.axes.plot(X, Y)
#

class Comp:
    def __init__(self):
        self.index  = 0
        self.mode  = 0
        self.col  = -1
        self.multi = 1e0    # Multiplication factor in MODE0
        self.shift = 0e0    # Constant to be addded in MODE0
        self.min = 1e1      # Minimum value for the parametrized concentration OR constant value if max <= min
        self.max = 1e5      # Peak value
        self.sig = 1e0      # Standard deviation for the Gaussian=sig of the bell curve
        self.mju = 12e0     # Time of peak value
        self.fv  = 0e0      # Angular frequency [hours] of modifying sine function
        self.ph  = 0e0      # Angular frequency [hours] of modifying sine function
        self.am  = 1e0      # Amplitude of modificaion
        self.name = 'NONAME'# Human readable name for modified variable
        self.unit = '#'     # unit name
        self.Find = 1
        self.pmInUse = 'no'

class QtBoxGui(gui5.Ui_MainWindow,QtWidgets.QMainWindow):
    """Main program window."""
    def __init__(self):
        super(QtBoxGui,self).__init__()
        # pg.mkPen('y', width=3, style=QtCore.Qt.DashLine)          ## Make a dashed yellow line 2px wide
        # pg.mkPen(0.5)                                             ## solid grey line 1px wide
        # pg.mkPen(color=(200, 200, 255), style=QtCore.Qt.DotLine)  ## Dotted pale-blue line
        # pg.setConfigOption('foreground', 'r')

        self.setupUi(self)
        # uic.loadUi('gui/LSD.ui', self)


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
        self.checkBox_aer.clicked.connect(lambda: self.grayIfNotChecked(self.checkBox_aer,self.frame_4))
        self.fsave_division.valueChanged.connect(self.toggle_printtime)

    # -----------------------
    # tab Input variables
    # -----------------------
        self.butMoveToSelVars.clicked.connect(self.select_compounds)
        self.butRemoveSelVars.clicked.connect(self.remv_item)
        self.selected_vars.setColumnHidden(7, True)

    # -----------------------
    # tab Function creator
    # -----------------------
        # plt_layout = QtWidgets.QVBoxLayout(self.PLOT)
        # scanvas = MyStaticMplCanvas(self.PLOT, width=5, height=4, dpi=100)
        # plt_layout.addWidget(scanvas)
        self.list_of_input = [self.runtime, self.run_name,self.case_dir]
        self.plotTo.clicked.connect(self.select_compounds_for_plot)
        self.plotTo_clear.clicked.connect(self.clearPlot)
        self.fMin.editingFinished.connect(self.updteGraph)
        self.fMax.editingFinished.connect(self.updteGraph)
        self.fLog.clicked.connect(self.updteGraph)
        self.fLin.clicked.connect(self.updteGraph)
        self.fWidth.valueChanged.connect(self.updteGraph)
        self.resW.clicked.connect(lambda: self.resetSlider(self.fWidth, 50))
        self.resP.clicked.connect(lambda: self.resetSlider(self.fPeak, 83))
        self.resA.clicked.connect(lambda: self.resetSlider(self.fFreq, 0))
        self.resPh.clicked.connect(lambda: self.resetSlider(self.fPhase, 0))
        self.resAm.clicked.connect(lambda: self.resetSlider(self.fAmp, 20))
        self.resG.clicked.connect(lambda: self.resetSlider(self.gain, 50))
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
        self.frameStop.setEnabled(False)
        self.startButton.clicked.connect(self.startBox)
        self.stopButton.clicked.connect(self.stopBox)
        self.boxProcess = 0 # Superbox run handle
        # self.r = 0 # Superbox output file handle
        self.Timer = QtCore.QTimer(self);
        self.pollTimer = QtCore.QTimer(self);
        self.Timer.timeout.connect(self.updateOutput)
        self.pollTimer.timeout.connect(self.pollMonitor)


    # -----------------------
    # Class methods
    # -----------------------
    def resetSlider(self, slider, pos):
        slider.setProperty("value", pos)

    def updteGraph(self):
        gain = 10**(self.gain.value()/50.-1)
        rt = self.runtime.value()
        wScale = rt/2/200.0 * gain
        pScale = rt*1.2/200.0 * gain
        aScale = 0.1 * gain
        phScale = rt/0.2/200.0 * gain
        ampScale = 1/20.0 * gain

        x = linspace(0,rt,200)
        yscale = self.radio(self.fLin, self.fLog)

        sigma = self.fWidth.value()*wScale
        if abs(sigma)<0.01: sigma = 0.01
        try:
            mini = float(self.fMin.text())
        except:
            mini = 0
        try:
            maxi = float(self.fMax.text())
        except:
            maxi=0
        peak = self.fPeak.value()*pScale
        freq = self.fFreq.value()*aScale
        phase = self.fPhase.value()*phScale
        amp = self.fAmp.value()*ampScale
        self.monW.setValue(sigma)
        self.monP.setValue(peak)
        self.monA.setValue(freq)
        self.monPh.setValue(phase)
        self.monAm.setValue(amp)
        mini  = dummy.min
        maxi  = dummy.max
        sigma = dummy.sig
        peak  = dummy.mju
        freq  = dummy.fv
        phase = dummy.ph
        amp   = dummy.am

        f = 1/sqrt(2*pi*sigma**2)
        D = peak + sin((x-peak)*freq)*amp + phase
        norm = 1/sqrt(2*pi*sigma**2)*exp(-(x-D)**2/(2*sigma**2))
        if yscale == 'lin':
            f = (maxi-mini)/f
            norm = norm*f + mini
        else:
            f = (log10(maxi-mini+1))/f
            norm = 10**(norm*f)-1 + mini

        try: # delete old legend if it exists:
            if self.plots > 0:
                self.legend.scene().removeItem(self.legend)
            self.plots = 1
        except Exception as e:
            print(e)

        self.legend = self.PLOT.addLegend()
        self.PLOT.plot(x,norm,pen=pg.mkPen('y', width=2), clear=True, name='work')
        if self.show_extra_plots != '':
            y=self.gauss(vars.mods[self.show_extra_plots], yscale,rt)
            self.PLOT.plot(x,y,pen=pg.mkPen(color=(200, 200, 255), width=1,style=QtCore.Qt.DotLine), name=self.show_extra_plots)

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
        close = QtWidgets.QPushButton()
        close.setFixedSize(60,30)
        close.setCheckable(True)

        close.setText('remove')
        cols = [text, '-1','1.0', '0.0','no']
        for i in range(4):
            self.selected_vars.setItem(row, i, QtWidgets.QTableWidgetItem(cols[i]))
        self.selected_vars.setCellWidget(row, i+1, pmInUse )
        self.selected_vars.setCellWidget(row, i+2, unit )
        self.selected_vars.setCellWidget(row, i+3, close )
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

    def clearPlot(self):
        self.show_extra_plots = ''
        self.updteGraph()

    def select_compounds_for_plot(self):
        compound = self.names_sel_2.selectedItems()[0]
        self.show_extra_plots = compound.text()
        self.updteGraph()

    def stopBox(self):
        self.Timer.stop()
        self.boxProcess.kill()
        tout = self.boxProcess.wait(timeout=10)
        self.boxProcess.poll()
        self.toggle_frame(self.frameStop)
        self.toggle_frame(self.frameStart)
        # try:
        #     self.r.close()
        # except:
        #     pass

    # def startBox(self):
    #     tmpfile = subprocess.Popen(['tee', 'process_diary.txt'], stdin=subprocess.PIPE, stdout=None).stdin
    #     self.boxProcess = subprocess.Popen(["./superbox.exe", " input/test"], stdout=tmpfile,stdin=subprocess.PIPE)
    #     self.Timer.start(10)
    #     self.r = open('process_diary.txt', 'r')
    #     self.MonitorWindow.clear()
    #
    #
    # def updateOutput(self):
    #     text2=self.r.readlines()
    #     fulltext = ''.join(text2)
    #     self.MonitorWindow.insertPlainText(fulltext)
    #     self.MonitorWindow.verticalScrollBar().setSliderPosition(self.MonitorWindow.verticalScrollBar().maximum());
    #
    #     if 'SIMULATION HAS ENDED' in fulltext[-50:]:
    #         self.MonitorWindow.setPlainText(self.MonitorWindow.toPlainText())
    #         self.MonitorWindow.verticalScrollBar().setSliderPosition(self.MonitorWindow.verticalScrollBar().maximum());
    #         self.stopBox()


    def startBox(self):
        self.print_values('input/tmp_b65d729f784bc8fcfb4beb009ac7a31d')
        try:
            self.boxProcess = subprocess.Popen(["./superbox.exe", " input/tmp_b65d729f784bc8fcfb4beb009ac7a31d"], stdout=subprocess.PIPE,stderr=subprocess.STDOUT,stdin=None)
            self.MonitorWindow.clear()
            self.Timer.start(10)
            self.pollTimer.start(1500)
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
        # nml.ENV.TEMPUNIT=0
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

if __name__ == '__main__':
    app = QtWidgets.QApplication([])
    qt_box = QtBoxGui()
    qt_box.show()
    app.exec_()
    # sys.exit(app.exec_())
