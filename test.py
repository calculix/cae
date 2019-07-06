"""
    Non modal dialod example
"""

import sys
from PyQt5 import QtWidgets

class Dialog2(QtWidgets.QDialog):
    def __init__(self, parent):
        QtWidgets.QDialog.__init__(self, parent)
        self.setModal(0)
        self.show()


def window():
    app = QtWidgets.QApplication(sys.argv)
    w = QtWidgets.QWidget()

    c = QtWidgets.QPushButton(w)
    c.setText("Nonmodal Dialog")
    c.move(15, 75)
    c.clicked.connect(lambda: Dialog2(w))

    w.resize(200, 150)
    w.show()
    sys.exit(app.exec_())


if __name__ == '__main__':
    window()