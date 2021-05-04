import sys
import re
from PySide6 import QtWidgets
from primers import primers, Primer
from PySide6 import QtCore, QtGui
from PySide6.QtWidgets import (QApplication, QComboBox, QDialog, QFormLayout,
                               QGridLayout, QLabel, QLineEdit, QMessageBox,
                               QPushButton, QSpinBox, QStyle, QTableView,
                               QPlainTextEdit, QWidget, QTabWidget, QTextEdit, QVBoxLayout, QHBoxLayout, QCheckBox, QRadioButton)
from tinydb import Query, TinyDB

db = TinyDB('db.json')

nuklpatt = re.compile(r"[ATGCatgc]+")

#logo
logo = 'assets/logo.svg'

#load database for ATTs
def loadfromdb(db: TinyDB,):
    atttable = db.table('ATTs')
    attlist: list[str] = []
    for att in atttable.all():
        attlist.append(att['name'])
    return attlist

#Pick lowest option for filling
def generatebest(target, fatt, ratt):
    #Generete all possible fillings
    mis1 = ["A", "T", "G", "C"]
    mis2 = ["AA", "AT", "AG", "AC", "TA", "TT", "TG", "TC", "GA", "GT", "GG", "GC", "CA", "CT", "CG", "CC"]
    #dummy hight temp primers
    prev_fwd = Primer(seq='', tm=500.50, tm_total=550, gc=0.5, dg=-0.1, fwd=True, offtargets=0, penalty=1.55)
    prev_rev = Primer(seq='', tm=500.50, tm_total=550, gc=0.5, dg=-0.1, fwd=True, offtargets=0, penalty=1.55)
    #F primer fillings
    if fatt['miss'] == 2:
        for add in mis2:
            fwd, rev = primers(target[:-3], fatt['seq']+add, ratt['seq'])
            if fwd.tm_total < prev_fwd.tm_total:
                prev_fwd = fwd
    elif fatt['miss'] == 1:
        for add in mis1:
            fwd, rev = primers(target[:-3], fatt['seq']+add, ratt['seq'])
            if fwd.tm_total < prev_fwd.tm_total:
                prev_fwd = fwd
    elif fatt['miss'] == 0:
        fwd, rev = primers(target[:-3], fatt['seq'], ratt['seq'])
        prev_fwd = fwd
    #R primer fillings
    if ratt['miss'] == 2:
        for add in mis2:
            fwd, rev = primers(target[:-3], fatt['seq'], ratt['seq']+add)
            if rev.tm_total < prev_rev.tm_total:
                prev_rev = rev
    elif ratt['miss'] == 1:
        for add in mis1:
            fwd, rev = primers(target[:-3], fatt['seq'], ratt['seq']+add)
            if rev.tm_total < prev_rev.tm_total:
                prev_rev = rev
    elif ratt['miss'] == 0:
        fwd, rev = primers(target[:-3], fatt['seq'], ratt['seq'])
        prev_rev = rev
    return(prev_fwd, prev_rev)

#Main window
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
            
        self.editbtf.clicked.connect(self.addeditbtn_press)
        self.compbt.clicked.connect(self.computeprimers)        

    def addeditbtn_press(self):            
        self.dialog = EditATTDB()
        self.dialog.exec_()
        self.choosattf.clear()       
        self.choosattr.clear()
        self.choosattf.addItems(loadfromdb(db)) 
        self.choosattr.addItems(loadfromdb(db))

    def computeprimers(self):
        self.inputtxt = self.textEdit.toPlainText().split("\n")
        self.tgseqs: list = self.parse_fasta(self.inputtxt)
        self.att = Query()
        self.atttable = db.table('ATTs')
        self.radd = self.atttable.search(self.att.name == self.choosattr.currentText())[0]
        self.fadd = self.atttable.search(self.att.name == self.choosattf.currentText())[0]
        if self.tgseqs:
            self.resultw = Result(targets=self.tgseqs, fadd=self.fadd, radd=self.radd)
            self.resultw.show()
        else:
            self.fastaempty()
        
    def closeEvent(self, event):
        reply = QMessageBox.question(self, 'Quit', 'Are you sure you want to quit?', QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
        if reply == QMessageBox.Yes:
            event.accept()
        else:
            event.ignore()

    #Empty FASTA error
    def fastaempty(self):
        msgfasta = QMessageBox()
        msgfasta.setWindowIcon(QtGui.QIcon(logo))
        msgfasta.setIcon(QMessageBox.Information)
        msgfasta.setText(f"FASTA text can't be empty")
        msgfasta.setWindowTitle("Empty FASTA")
        msgfasta.setStandardButtons(QMessageBox.Ok)
        msgfasta.exec_()

    #Unknown nuclide error
    def msgatt(self, index):         
        msgatt = QMessageBox()
        msgatt.setWindowIcon(QtGui.QIcon(logo))
        msgatt.setIcon(QMessageBox.Information)
        msgatt.setText(f"In sequence #{index}, is a unknown nuclide")
        msgatt.setWindowTitle("Unknown nuclide")
        msgatt.setStandardButtons(QMessageBox.Ok)
        msgatt.exec_()

    #parsing FASTA to dict of targets
    def parse_fasta (self, lines):
        self.descs: list = []
        self.seqs: list = []
        self.data = ''
        self.dic = []
        for line in lines:
            if line.startswith('>'):
                if self.data:   # have collected a sequence, push to seqs
                    self.seqs.append(self.data)
                    self.data = ''
                self.descs.append(line[1:])  # Trim '>' from beginning
            else:
                self.data += line.rstrip('\r\n')
        # there will be yet one more to push when we run out
        self.seqs.append(self.data)
        for i in self.descs:
            indx: int = self.descs.index(i)            
            if not re.fullmatch(nuklpatt, self.seqs[indx]):
                self.msgatt(indx+1)
            if re.fullmatch(nuklpatt, self.seqs[indx]):
                self.dic.append({'name': self.descs[indx], 'seq': self.seqs[indx]})
        return self.dic

#Model of ATT tab
class ATTTableModel(QtCore.QAbstractTableModel):
    header_labels = ['Name', 'Seq', 'Addon']
    def __init__(self, data):
        super(ATTTableModel, self).__init__()
        self._data = data

    def data(self, index, role):
        if role == QtCore.Qt.DisplayRole:
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

#Dialog of table of ATTs in DB
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
        if not self.data:
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
        self.edtbt.setEnabled(False)
        ltab.addWidget(self.edtbt, 1, 1)

        self.delbt = QPushButton("Delete")
        self.delbt.setEnabled(False)
        ltab.addWidget(self.delbt, 1, 2)

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

#Dialog for adding (editing) of ATT
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
        self.attregex = QtCore.QRegularExpression("[ATGCatgc]+")
        self.attvalidator = QtGui.QRegularExpressionValidator(self.attregex)
        self.seqet.setValidator(self.attvalidator)
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
        att = Query()
        atttable = db.table('ATTs')
        if self.namet.text().strip() and self.seqet.text().strip():
            if not atttable.search(att.name == self.namet.text()):
                atttable.insert({'name': self.namet.text(), 'seq': self.seqet.text().upper(), 'miss': self.miset.value()})
                self.close()
            else:
                self.infoduplicate(self.namet.text())
        else:
            self.infoempty()
            self.close()

    def infoduplicate(self, name: str):
        self.msg = QMessageBox()
        self.msg.setWindowIcon(QtGui.QIcon(logo))
        self.msg.setIcon(QMessageBox.Information)
        self.msg.setText(f"Theres is sequence with name: {name}")
        self.msg.setWindowTitle("Duplicate name")
        self.msg.setStandardButtons(QMessageBox.Ok)
        self.msg.exec_()

    def infoempty(self):        
        self.msg = QMessageBox()
        self.msg.setWindowIcon(QtGui.QIcon(logo))
        self.msg.setIcon(QMessageBox.Information)
        self.msg.setText("Name and seq can't be empty! \nNothing was saved.")
        self.msg.setWindowTitle("Empty field")
        self.msg.setStandardButtons(QMessageBox.Ok)
        self.msg.exec_()

#FASTA Text edit field
class TextEdit(QPlainTextEdit):
    def __init__(self, parent=None):
        super(TextEdit, self).__init__(parent)

        self.setPlainText("Enter target sequence in FASTA")
        self.ffocus: bool = True

    def focusInEvent(self, event):
        if self.ffocus is True:
            self.clear()
            self.ffocus = False
        QPlainTextEdit.focusInEvent(self, event)

#Result window, shows in tab primers sequences and their temp
class Result(QTabWidget):
    def __init__(self, targets, fadd, radd, parent = None):
        super(Result, self).__init__(parent)        
        self.setWindowTitle("Results")
        self.setWindowIcon(QtGui.QIcon(logo))
        self.resize(400, 150)
        self.prims = []
        for targ in targets:
            fwd, rev = generatebest(targ["seq"], fadd, radd)
            self.prims.append({'target': targ['name'], 'Fprimer': fwd.seq, 'Ftemp': fwd.tm_total, 'Rprimer': rev.seq, 'Rtemp': rev.tm_total})
		
        for prim in self.prims:
            self.tab = QWidget()
            self.addTab(self.tab, prim['target'])
            self.tabUI(prim)
		
    def tabUI(self, prim):
        self.fprim = prim['Fprimer']
        self.ftemp = str(f"{prim['Ftemp']} °C")        
        self.rprim = prim['Rprimer']
        self.rtemp = str(f"{prim['Rtemp']} °C")

        layout = QGridLayout(self)
        self.tab.setLayout(layout)

        #R primers
        self.fplb = QLabel("Fprimer:")
        layout.addWidget(self.fplb, 0, 0)
        self.fpte = QLineEdit(self.fprim)
        self.fpte.setReadOnly(True)
        layout.addWidget(self.fpte, 0, 1)
        self.fplbtm = QLabel("Fprimer TM_total:")
        layout.addWidget(self.fplbtm, 1, 0)
        self.fptemp = QLabel(self.ftemp)
        layout.addWidget(self.fptemp, 1, 1)

        #R primers
        self.rplb = QLabel("Rprimer:")
        layout.addWidget(self.rplb, 2, 0)
        self.rpte = QLineEdit(self.rprim)
        self.rpte.setReadOnly(True)
        layout.addWidget(self.rpte, 2, 1)
        self.rplbtm = QLabel("Rprimer TM_total:")
        layout.addWidget(self.rplbtm, 3, 0)
        self.rptemp = QLabel(self.rtemp)
        layout.addWidget(self.rptemp, 3, 1)
        

def main():
    app = QApplication(sys.argv)
    win = StartWin()
    win.show()
    sys.exit(app.exec_())

if __name__ == "__main__":
    main()
