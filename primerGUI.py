import sys
from PySide6 import QtWidgets

import qtawesome as qta
from PySide6 import QtCore, QtGui
from PySide6.QtWidgets import (QApplication, QComboBox, QDialog, QFormLayout,
                               QGridLayout, QLabel, QLineEdit, QMessageBox,
                               QPushButton, QSpinBox, QStyle, QTableView,
                               QTextEdit, QWidget)
from tinydb import Query, TinyDB

db = TinyDB('db.json')

#logo
logo = 'logo.svg'

# add att to DB

def loadfromdb(db: TinyDB,):
    atttable = db.table('ATTs')
    attlist: list[str] = []
    for att in atttable.all():
        attlist.append(att['name'])
    return attlist

class StartWin(QWidget):   
    def __init__(self):
        super().__init__()
        self.initUI()

    def initUI(self):

        #Options of window
        self.setWindowTitle("Primer Designer")
        self.setWindowIcon(QtGui.QIcon(logo))
            
        #Main layout
        layout = QGridLayout(self)
        self.setLayout(layout)

        self.wellab = QLabel("Welcome to primer designer")
        layout.addWidget(self.wellab, 0, 0),

        self.specattf = QLabel(f"Selected ATT part for F:")
        layout.addWidget(self.specattf, 1, 0)
        self.choosattf = QComboBox()
        self.choosattf.addItems(loadfromdb(db))
        layout.addWidget(self.choosattf, 2, 0)
        self.editbtf = QPushButton("Add/Edit ATTs")
        self.editbtf.setIcon(self.style().standardIcon(getattr(QStyle, 'SP_FileDialogContentsView')))
        layout.addWidget(self.editbtf, 2, 1, 3, 1)     

        self.specattr = QLabel(f"Selected ATT part for R:")
        layout.addWidget(self.specattr, 3, 0)
        self.choosattr = QComboBox()
        self.choosattr.addItems(loadfromdb(db))
        layout.addWidget(self.choosattr, 4, 0)
            
        self.enterl = QLabel("Target sequence:")
        layout.addWidget(self.enterl, 5, 0)
        self.textEdit = TextEdit()
        layout.addWidget(self.textEdit, 6, 0, 1, 2)
            
        self.compbt = QPushButton("Compute")
        layout.addWidget(self.compbt, 7, 0, 1, 2)
            
        #btn1 press event
        self.editbtf.clicked.connect(self.btnpress1_clicked)

        # add att but
    def btnpress1_clicked(self):            
        self.dialog = EditATTDB()
        self.dialog.exec_()
        self.choosattf.clear()       
        self.choosattr.clear()
        self.choosattf.addItems(loadfromdb(db)) 
        self.choosattr.addItems(loadfromdb(db))
        
    def closeEvent(self, event):
        reply = QMessageBox.question(self, 'Quit', 'Are you sure you want to quit?', QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
        if reply == QMessageBox.Yes:
            event.accept()
        else:
            event.ignore()

class ATTTableModel(QtCore.QAbstractTableModel):
    header_labels = ['Name', 'Seq', 'Addon']
    def __init__(self, data):
        super(ATTTableModel, self).__init__()
        self._data = data

    def data(self, index, role):
        if role == QtCore.Qt.DisplayRole:
            # See below for the nested-list data structure.
            # .row() indexes into the outer list,
            # .column() indexes into the sub-list
            return self._data[index.row()][index.column()]
    
    def headerData(self, section, orientation, role=QtCore.Qt.DisplayRole):
        if role == QtCore.Qt.DisplayRole and orientation == QtCore.Qt.Horizontal:
            return self.header_labels[section]
        return QtCore.QAbstractTableModel.headerData(self, section, orientation, role)

    def rowCount(self, index):
        # The length of the outer list.
        return len(self._data)

    def columnCount(self, index):
        # The following takes the first sub-list, and returns
        # the length (only works if all rows are an equal length)
        return len(self._data[0])

class EditATTDB(QDialog):
    def __init__(self, parent=None):
        super(EditATTDB, self).__init__(parent)
        self.setWindowTitle("Edit ATT DB")
        self.setWindowIcon(QtGui.QIcon(logo))

        ltab = QGridLayout(self)
        self.setLayout(ltab)

        self.table = QTableView()
        self.table.setSizeAdjustPolicy(QtWidgets.QAbstractScrollArea.AdjustToContents)
        self.atttable = db.table('ATTs')
        self.data = []
        for att in self.atttable:
            self.data.append([att['name'], att['seq'], att['miss']])
        if len(self.data) == 0:
            self.data.append(["nothing", "please add", "some data"])
        self.model = ATTTableModel(self.data)
        self.table.setModel(self.model)
        self.table.clicked.connect(self.viewclicked)
        self.table.setSelectionBehavior(QTableView.SelectRows)
        ltab.addWidget(self.table, 0, 0, 1, 3)

        self.addbt = QPushButton("Add")
        self.addbt.clicked.connect(self.addbtn)
        ltab.addWidget(self.addbt, 1, 0)
        
        self.edtbt = QPushButton("Edit")
        ltab.addWidget(self.edtbt, 1, 1)
        
        self.edtbt = QPushButton("Delete")
        ltab.addWidget(self.edtbt, 1, 2)

    def viewclicked(self, clickedIndex):
        row=clickedIndex.row()
        model=clickedIndex.model()

    def addbtn(self):
        self.diag = AddATT()
        self.diag.exec_()
        self.data.clear()
        for att in self.atttable:
            self.data.append([att['name'], att['seq'], att['miss']])
        self.model = ATTTableModel(self.data)
        self.table.setModel(self.model)


class AddATT(QDialog):
    def __init__(self, parent=None):
        super(AddATT, self).__init__(parent)
        self.setWindowTitle("Add ATT to DB")
        self.setWindowIcon(QtGui.QIcon(logo))

        self.lt = QFormLayout(self)

        self.name: str
        self.naml = QLabel("ATT name:")
        self.namet = QLineEdit()
        self.lt.addRow(self.naml, self.namet)

        self.seql = QLabel("ATT seq:")
        self.seqet = QLineEdit()
        self.lt.addRow(self.seql, self.seqet)

        self.misl = QLabel("# of nuklids to be added:")
        self.miset = QSpinBox()
        self.miset.setMaximum(2)
        self.miset.setMinimum(0)
        self.lt.addRow(self.misl, self.miset)

        self.subbt = QPushButton("Add")
        self.lt.addRow(self.subbt)
        self.subbt.clicked.connect(self.addatttodb)

    def addatttodb(self):
        atttable = db.table('ATTs')
        if self.namet.text().strip() and self.seqet.text().strip():
            atttable.insert({'name': self.namet.text(), 'seq': self.seqet.text(), 'miss': self.miset.value()})
        else:
            self.infoempty()
        self.close()

    def infoempty(self):        
        self.msg = QMessageBox()
        self.msg.setWindowIcon(QtGui.QIcon(logo))
        self.msg.setIcon(QMessageBox.Information)
        self.msg.setText("Name and seq can't be empty! \nNothing was saved.")
        self.msg.setWindowTitle("Empty field")
        self.msg.setStandardButtons(QMessageBox.Ok)
        self.msg.exec_()

class TextEdit(QTextEdit):
    def __init__(self, parent=None):
        super(TextEdit, self).__init__(parent)

        self.setText("Enter target sequence in FASTA")

    def focusInEvent(self, event):
        self.clear()
        QTextEdit.focusInEvent(self, event)

def main():
    app = QApplication(sys.argv)
    win = StartWin()
    win.show()
    sys.exit(app.exec_())

if __name__ == "__main__":
    main()
