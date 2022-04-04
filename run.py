import sys
from PyQt5.QtWidgets import QApplication, QMainWindow, QAction, QActionGroup
from engine import Grf
from rdkit import Chem


class Wins(QMainWindow):

    def __init__(self):
        super().__init__()

        self.wid = Grf()
        toolbar = self.addToolBar('Experiment')
        getSm = QAction(':)',self)
        splin = QAction('~',self)
        un = QAction('<',self)
        re = QAction('>',self)
        cl = QAction('Clean',self)

        toolbar.addAction(getSm)
        toolbar.addAction(splin)
        toolbar.addAction(un)
        toolbar.addAction(re)
        toolbar.addAction(cl)

        self.gr = QActionGroup(self)
        self.actions = {}
        action_list = ['1', '3', '4', '5', '6', 'benz']
        for i in action_list:
            self.actions[i] = QAction(i,self, checkable=True)
            self.gr.addAction(self.actions[i])
            toolbar.addAction(self.actions[i])
            self.actions[i].triggered.connect(self.callBot)

        self.actions['1'].setChecked(True)

        getSm.triggered.connect(self.getS)
        splin.triggered.connect(self.spl)
        re.triggered.connect(self.redo)
        un.triggered.connect(self.undo)
        cl.triggered.connect(self.clean_)

        self.setCentralWidget(self.wid)
    
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

    def getS(self):
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

        print(Chem.MolToSmiles(mol))
        return mol

    def spl(self):
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