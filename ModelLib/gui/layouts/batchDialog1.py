# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'batchDialog1.ui'
#
# Created by: PyQt5 UI code generator 5.13.2
#
# WARNING! All changes made in this file will be lost!


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_batchDialog(object):
    def setupUi(self, batchDialog):
        batchDialog.setObjectName("batchDialog")
        batchDialog.resize(713, 457)
        self.gridLayout = QtWidgets.QGridLayout(batchDialog)
        self.gridLayout.setObjectName("gridLayout")
        self.tb_1 = QtWidgets.QPlainTextEdit(batchDialog)
        self.tb_1.setTextInteractionFlags(QtCore.Qt.TextSelectableByKeyboard|QtCore.Qt.TextSelectableByMouse)
        self.tb_1.setObjectName("tb_1")
        self.gridLayout.addWidget(self.tb_1, 1, 0, 1, 1)
        self.label_1 = QtWidgets.QLabel(batchDialog)
        self.label_1.setObjectName("label_1")
        self.gridLayout.addWidget(self.label_1, 0, 0, 1, 1)
        self.bDialogbuttonBox = QtWidgets.QDialogButtonBox(batchDialog)
        self.bDialogbuttonBox.setOrientation(QtCore.Qt.Horizontal)
        self.bDialogbuttonBox.setStandardButtons(QtWidgets.QDialogButtonBox.Abort|QtWidgets.QDialogButtonBox.SaveAll)
        self.bDialogbuttonBox.setCenterButtons(True)
        self.bDialogbuttonBox.setObjectName("bDialogbuttonBox")
        self.gridLayout.addWidget(self.bDialogbuttonBox, 2, 0, 1, 1)

        self.retranslateUi(batchDialog)
        QtCore.QMetaObject.connectSlotsByName(batchDialog)

    def retranslateUi(self, batchDialog):
        _translate = QtCore.QCoreApplication.translate
        batchDialog.setWindowTitle(_translate("batchDialog", "Batch Preview"))
        self.label_1.setText(_translate("batchDialog", "label_1"))
