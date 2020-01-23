# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'batchDialog.ui'
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
        self.bDialogDirTextEdit = QtWidgets.QPlainTextEdit(batchDialog)
        self.bDialogDirTextEdit.setTextInteractionFlags(QtCore.Qt.TextSelectableByKeyboard|QtCore.Qt.TextSelectableByMouse)
        self.bDialogDirTextEdit.setObjectName("bDialogDirTextEdit")
        self.gridLayout.addWidget(self.bDialogDirTextEdit, 1, 0, 1, 1)
        self.label_2 = QtWidgets.QLabel(batchDialog)
        self.label_2.setObjectName("label_2")
        self.gridLayout.addWidget(self.label_2, 2, 0, 1, 1)
        self.bDialogFileTextEdit = QtWidgets.QPlainTextEdit(batchDialog)
        self.bDialogFileTextEdit.setTextInteractionFlags(QtCore.Qt.TextSelectableByKeyboard|QtCore.Qt.TextSelectableByMouse)
        self.bDialogFileTextEdit.setObjectName("bDialogFileTextEdit")
        self.gridLayout.addWidget(self.bDialogFileTextEdit, 3, 0, 1, 1)
        self.label = QtWidgets.QLabel(batchDialog)
        self.label.setObjectName("label")
        self.gridLayout.addWidget(self.label, 0, 0, 1, 1)
        self.bDialogbuttonBox = QtWidgets.QDialogButtonBox(batchDialog)
        self.bDialogbuttonBox.setOrientation(QtCore.Qt.Horizontal)
        self.bDialogbuttonBox.setStandardButtons(QtWidgets.QDialogButtonBox.Abort|QtWidgets.QDialogButtonBox.SaveAll)
        self.bDialogbuttonBox.setCenterButtons(True)
        self.bDialogbuttonBox.setObjectName("bDialogbuttonBox")
        self.gridLayout.addWidget(self.bDialogbuttonBox, 4, 0, 1, 1)

        self.retranslateUi(batchDialog)
        QtCore.QMetaObject.connectSlotsByName(batchDialog)

    def retranslateUi(self, batchDialog):
        _translate = QtCore.QCoreApplication.translate
        batchDialog.setWindowTitle(_translate("batchDialog", "Batch Preview"))
        self.label_2.setText(_translate("batchDialog", "The following files will be created:"))
        self.label.setText(_translate("batchDialog", "The following directories will be created (existing directories are not listed):"))
