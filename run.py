try:
    from PySide2 import QtCore, QtWidgets
except:
    try:
        from PyQt5 import QtCore, QtWidgets
    except:
        print("You need PySide2 or PyQt5 to run this program")
from gui import guitest
import subprocess

import matplotlib
# Make sure that we are using QT5
matplotlib.use('Qt5Agg')
from numpy import linspace, sin
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

# For guitest.py:
# try:
#     from PySide2 import QtCore, QtGui, QtWidgets
# except:
#     try:
#         from PyQt5 import QtCore, QtGui, QtWidgets
#     except:
#         print("You need PySide2 or PyQt5 to run this program")
# #

## Some constants --------------------------------------------
X = linspace(0,3,100)
Y = sin(X)
## -----------------------------------------------------------

## Some constants --------------------------------------------
cw = [140,80,80,80,70,50,60]
all_units = [['K','C'],['Pa','hPa','bar','kPa', 'mbar'],['#','ppm','ppb','ppt','ppq']]

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

class MyMplCanvas(FigureCanvas):
    """Ultimately, this is a QWidget (as well as a FigureCanvasAgg, etc.)."""

    def __init__(self, parent, width, height, dpi):
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = fig.add_subplot(111)
        # self.axes.set_yscale('log')
        self.axes.set_position([0.1, 0.1, 0.88, 0.88])
        self.axes.grid()
        self.compute_initial_figure()

        FigureCanvas.__init__(self, fig)
        self.setParent(parent)
        FigureCanvas.setSizePolicy(self,
               QtWidgets.QSizePolicy.Expanding,
               QtWidgets.QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)

    # def compute_initial_figure(self):
    #     pass


class MyStaticMplCanvas(MyMplCanvas):
    """Simple canvas with a sine plot."""
    def compute_initial_figure(self):
        self.axes.plot(X, Y)


class MyQtApp(guitest.Ui_MainWindow, QtWidgets.QMainWindow):
    """Main program window."""
    def __init__(self):
        super(MyQtApp,self).__init__()
        self.setupUi(self)
        self.setWindowTitle('Superbox configurator')
        self.get_Names()
        # self.create_plot_area()
        # self.showMaximized()
        l = QtWidgets.QVBoxLayout(self.PLOT)
        sc = MyStaticMplCanvas(self.PLOT, width=5, height=4, dpi=100)
        l.addWidget(sc)

        self.pushButton_10.clicked.connect(self.select_compounds)
        for i in range(7):
            self.selected_vars.setColumnWidth(i, cw[i])

        self.add_new_line('TEMPK', 0)
        self.add_new_line('PRESSURE', 1)
        self.list_of_input = [self.lineEdit, self.lineEdit_2,self.lineEdit_3,self.lineEdit_4,self.lineEdit_5, self.spinBox_2]
        self.pushButton_2.clicked.connect(self.print_values)

        self.toolButton_4.clicked.connect(lambda: self.browse_folder(self.lineEdit_12))
        self.toolButton_5.clicked.connect(lambda: self.browse_folder(self.lineEdit_13))

        self.toolButton_2.clicked.connect(lambda: self.browse_file(self.lineEdit_10))
        self.toolButton_3.clicked.connect(lambda: self.browse_file(self.lineEdit_11))
        self.toolButton_6.clicked.connect(lambda: self.browse_file(self.lineEdit_14))
        self.toolButton_9.clicked.connect(lambda: self.browse_file(self.lineEdit_17))

        self.toolButton_7.clicked.connect(lambda: self.browse_file(self.lineEdit_15))
        self.toolButton_8.clicked.connect(lambda: self.browse_file(self.lineEdit_16))
        self.pushButton_3.clicked.connect(self.save_file)
        self.actionSave_2.triggered.connect(self.save_file)
        self.actionPrint.triggered.connect(self.print_values)
        self.actionOpen.triggered.connect(self.browse_file)
        self.actionQuit_Ctrl_Q.triggered.connect(self.close)

        self.checkBox_aer.clicked.connect(lambda: self.toggle_frame(self.checkBox_aer,self.frame_4))
        self.checkBox_che_4.clicked.connect(lambda: self.toggle_gray(self.checkBox_che_4,self.gridLayout_4))
        self.checkBox_acd_7.clicked.connect(lambda: self.toggle_gray(self.checkBox_acd_7,self.gridLayout_11))
        self.spinBox.valueChanged.connect(self.toggle_printtime)

        self.readOut.clicked.connect(self.startBox)
        self.readOut_2.clicked.connect(self.stopBox)
        self.boxProcess = 0
        self.r = 0
        self.Timer = QtCore.QTimer(self);
        self.Timer.timeout.connect(self.updateOutput)

    def browse_folder(self, target):
        dialog = QtWidgets.QFileDialog()
        dialog.setFileMode(QtWidgets.QFileDialog.Directory)
        dialog.setOption(QtWidgets.QFileDialog.ShowDirsOnly)
        options = dialog.Options()
        options |= dialog.DontUseNativeDialog
        directory = dialog.getExistingDirectory(self, 'Choose Directory', options=options)
        target.insert(directory)

    def browse_file(self, target):
        dialog = QtWidgets.QFileDialog()
        options = dialog.Options()
        options |= dialog.DontUseNativeDialog
        file = dialog.getOpenFileName(self, 'Choose File', options=options)[0]
        if file != '':
            target.insert(file)

    def save_file(self):
        dialog = QtWidgets.QFileDialog()
        options = dialog.Options()
        options |= dialog.DontUseNativeDialog
        file = dialog.getSaveFileName(self, 'Save INITFILE', options=options)[0]
        if file != '':
            f = open(file, 'w')
            f.write('Hello World!')
            f.close()

    def get_Names(self):
        self.namesdat.clear()
        self.namesdat.addItems(NAMES)

    def add_new_line(self, name, unit_ind):
        row = self.selected_vars.rowCount()
        self.selected_vars.insertRow(row)
        item = self.namesdat.item(namesPyInds[name])
        item.setFlags(item.flags() & ~QtCore.Qt.ItemIsEnabled & ~QtCore.Qt.ItemIsSelectable)
        text = '%s(%d)'%(name, namesFoInds[name] )
        unit = QtWidgets.QComboBox()
        unit.addItems(all_units[unit_ind])
        close = QtWidgets.QPushButton()
        close.setFixedSize(60,30)
        close.setText('remove')
        cols = [text, '-1','0', '1.0','no']
        for i in range(5):
            self.selected_vars.setItem(row, i, QtWidgets.QTableWidgetItem(cols[i]))
        self.selected_vars.setCellWidget(row, i+1, unit )
        self.selected_vars.setCellWidget(row, i+2, close )


    def toggle_gray(self, guard, group):
        index = group.count()
        for i in range(index):
            grayWidget = group.itemAt(i).widget()
            if guard.isChecked() == True:
                grayWidget.setEnabled(True)
            else:
                grayWidget.setEnabled(False)

    def toggle_frame(self, guard, frame):
        if guard.isChecked() == True:
            frame.setEnabled(True)
        else:
            frame.setEnabled(False)

    def toggle_printtime(self):
        if self.spinBox.value() != 0:
            self.spinBox_3.setEnabled(False)
        else:
            self.spinBox_3.setEnabled(True)

    def print_values(self):
        for i,v in enumerate(self.list_of_input):
            print(v.text())

    def select_compounds(self):
        compounds = self.namesdat.selectedItems()
        for c in compounds:
            self.add_new_line(c.text(), 2)

    def stopBox(self):
        self.Timer.stop()
        self.boxProcess.kill()
        tout = self.boxProcess.wait(timeout=10)
        self.boxProcess.poll()

    def startBox(self):
        tmpfile = subprocess.Popen(['tee', 'process_diary.txt'], stdin=subprocess.PIPE, stdout=None).stdin
        self.boxProcess = subprocess.Popen(["./superbox.exe", " input/test"], stdout=tmpfile,stdin=subprocess.PIPE)
        self.Timer.start(10)
        self.r = open('process_diary.txt', 'r')


    def updateOutput(self):
        text2=self.r.readlines()
        # f.close()
        # self.plainTextEdit_2.clear()
        self.plainTextEdit_2.insertPlainText(''.join(text2))
        self.plainTextEdit_2.ensureCursorVisible()


if __name__ == '__main__':
    app = QtWidgets.QApplication([])
    qt_app = MyQtApp()
    qt_app.show()
    app.exec_()
    # sys.exit(app.exec_())
