import sys
from PyQt5.QtWidgets import QApplication, QMainWindow, QDockWidget, QAction, QActionGroup, QErrorMessage, QMessageBox
try:
    from .engine import *
except:
    from engine import *
from rdkit import Chem
import traceback


class MolDrawer(QDockWidget):

    def __init__(self, caller, smiles , parent = None):
        super().__init__(parent)

        self.caller = caller

        self.window = QMainWindow()

        self.wid = Grf()
        toolbar = self.window.addToolBar('Experiment')
        toolbar.setContextMenuPolicy(Qt.PreventContextMenu)
        toolbar.setMovable(False)
        getSm = QAction(get_icon('update'), 'Update',self)
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
        self.window.setCentralWidget(self.wid)
        self.setWidget(self.window)
        self.setWindowIcon(get_icon('flask'))
        self.setWindowTitle("MolDrawer")

        mol = Chem.MolFromSmiles(smiles)

        try:
            self.wid.addMol(mol)
        except:
            self.error_d = QErrorMessage()
            self.error_d.showMessage(traceback.format_exc())
    
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
            self.caller(mol)
            self.close()

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
            self.error_d.showMessage('Bad structure!')#'Все хуйня! Давай по новой!!!')

    def spl(self):
        if self.wid.ErrorList:
            self.error_d = QErrorMessage()
            self.error_d.showMessage('Bad structure!')#'Все хуйня! Давай по новой!!!')
        else:
            mol = self.getS()
            self.wid.addMol(mol)

    def undo(self):
        self.wid.history.undo()
        self.wid.ErrorList = self.wid.findErrors()

    def redo(self):
        self.wid.history.redo()
        self.wid.ErrorList = self.wid.findErrors()


if __name__ == '__main__':
    app = QApplication(sys.argv)
    app.setStyle('Fusion')
    ex = MolDrawer(print, '')
    ex.show()
    app.exec_()