# -*- coding: utf-8 -*-

"""
    © Ihor Mirzov, UJV Rez, June 2019.
    Distributed under GNU General Public License, version 2.

    Dialog window to edit CalculiX model
"""


from PyQt5 import QtWidgets


class Dialog(QtWidgets.QDialog):  
    def __init__(self, dinput):
        super(Dialog, self).__init__()
        self.createFormGroupBox(dinput)        

        buttonBox = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.Ok | QtWidgets.QDialogButtonBox.Cancel)
        buttonBox.accepted.connect(self.accept)
        buttonBox.rejected.connect(self.reject)

        mainLayout = QtWidgets.QVBoxLayout(self)
        mainLayout.addWidget(self.formGroupBox)        
        mainLayout.addWidget(buttonBox)     

        self.setWindowTitle("Form Layout")    

    def createFormGroupBox(self, dinput):
        layout = QtWidgets.QFormLayout()
        self.linedit1 = QtWidgets.QLineEdit('TestName')
        self.combox1 = QtWidgets.QComboBox()
        self.combox1.setToolTip('Hello')
        self.combox1.addItems(['India','France','UK','USA','Germany'])
        self.spinbox1 = QtWidgets.QSpinBox()  

        for text, w in zip(dinput, (self.linedit1, self.combox1, self.spinbox1)):
            layout.addRow(text, w)     

        self.formGroupBox = QtWidgets.QGroupBox("Form layout")        
        self.formGroupBox.setLayout(layout)

    def accept(self):        
        self._output = self.linedit1.text(), self.combox1.currentText(), self.spinbox1.value()    
        super(Dialog, self).accept()

    def get_output(self):
        return self._output

