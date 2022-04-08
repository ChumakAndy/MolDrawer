import sys
from PyQt5.QtWidgets import QApplication, QMainWindow, QAction, QActionGroup, QErrorMessage, QMessageBox
from engine import *
from rdkit import Chem


class Wins(QMainWindow):

    def __init__(self):
        super().__init__()

        self.wid = Grf()
        toolbar = self.addToolBar('Experiment')
        toolbar.setContextMenuPolicy(Qt.PreventContextMenu)
        toolbar.setMovable(False)
        getSm = QAction(get_icon('smile'), ':)',self)
        splin = QAction(get_icon('sanitize'), 'Sanitize',self)
        un = QAction(get_icon('undo'), 'Undo',self)
        re = QAction(get_icon('redo'), 'Redo',self)
        cl = QAction(get_icon('clear'), 'Clean',self)

        toolbar.addAction(getSm)
        toolbar.addAction(splin)
        toolbar.addAction(un)
        toolbar.addAction(re)
        toolbar.addAction(cl)

        self.gr = QActionGroup(self)
        self.actions = {}
        action_list = ['drag', '1', '3', '4', '5', '6', 'benz']
        for i in action_list:
            self.actions[i] = QAction(get_icon(i), i , self, checkable=True)
            self.gr.addAction(self.actions[i])
            toolbar.addAction(self.actions[i])
            self.actions[i].triggered.connect(self.callBot)

        self.actions['1'].setChecked(True)

        getSm.triggered.connect(self.getSmiles)
        splin.triggered.connect(self.spl)
        re.triggered.connect(self.redo)
        un.triggered.connect(self.undo)
        cl.triggered.connect(self.clean_)

        self.setCentralWidget(self.wid)
        self.setWindowIcon(get_icon('flask'))
        self.setWindowTitle("MolDrawer")
    
    def clean_(self):
        self.wid.scene().clear()
        self.wid.atomList = set()
        self.wid.boundSet = set()
        self.wid.history.data  = []
        self.wid.history.pos = -1
        self.wid.history.update()
    
    def callBot(self):
        sender = self.sender().text()
        self.wid.mode = sender
    
    def getSmiles(self):
        if self.wid.ErrorList:
            self.error_d = QErrorMessage()
            self.error_d.showMessage('Bad structure!')
        else:
            mol = self.getS()
            print(Chem.MolToSmiles(mol))
            self.msg = QMessageBox()
            self.msg.setIcon(QMessageBox.Information)
            self.msg.setText("Smiles")
            self.msg.setInformativeText(Chem.MolToSmiles(mol))
            self.msg.setWindowTitle("MolDrawer")
            self.msg.setWindowIcon(get_icon('flask'))
            self.msg.show()

    def getS(self):
        try:
            bkinds = {1: Chem.BondType.SINGLE,
                    2: Chem.BondType.DOUBLE,
                    3: Chem.BondType.TRIPLE}
            mol = Chem.RWMol()
            list_atoms = list(self.wid.atomList)

            for atom in list_atoms:
                atom.IDX = list_atoms.index(atom)
                rdatom = Chem.Atom(atom.kind)
                if atom.charge != 0:
                    rdatom.SetFormalCharge(atom.charge)
                mol.AddAtom(rdatom)

            for bond in self.wid.boundSet:
                at1 = bond.atoms[0].IDX
                at2 = bond.atoms[1].IDX
                mol.AddBond(at1, at2, bkinds[bond.multiplicity])
            
            mol = mol.GetMol()

            Chem.SanitizeMol(mol)
            Chem.Kekulize(mol, clearAromaticFlags=True)
            return mol
        except:
            self.error_d = QErrorMessage()
            self.error_d.showMessage('Bad structure!')

    def spl(self):
        if self.wid.ErrorList:
            self.error_d = QErrorMessage()
            self.error_d.showMessage('Bad structure!')
        else:
            mol = self.getS()
            self.wid.addMol(mol)

    def undo(self):
        self.wid.history.undo()

    def redo(self):
        self.wid.history.redo()





if __name__ == '__main__':
    app = QApplication(sys.argv)
    app.setStyle('Fusion')
    ex = Wins()
    ex.show()
    app.exec_()