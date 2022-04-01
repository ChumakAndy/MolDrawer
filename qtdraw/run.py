import sys
from PyQt5.QtWidgets import QApplication, QMainWindow, QAction
from engine import Grf
from rdkit import Chem


class Wins(QMainWindow):

    def __init__(self):
        super().__init__()

        self.wid = Grf()
        toolbar = self.addToolBar('Experiment')
        getSm = QAction('Smiles',self)
        toolbar.addAction(getSm)
        getSm.triggered.connect(self.getS)

        self.setCentralWidget(self.wid)
    
    def getS(self):
        bkinds = {1: Chem.BondType.SINGLE,
                  2: Chem.BondType.DOUBLE,
                  3: Chem.BondType.TRIPLE}
        mol = Chem.RWMol()
        list_atoms = list(self.wid.atomList)

        for atom in list_atoms:
            atom.IDX = list_atoms.index(atom)
            rdatom = Chem.Atom(atom.kind)
            # if atom.charge != '':
            #     rdatom.SetFormalCharge(ckinds[atom.charge])
            mol.AddAtom(rdatom)

        for bond in self.wid.boundSet:
            at1 = bond.atoms[0].IDX
            at2 = bond.atoms[1].IDX
            mol.AddBond(at1, at2, bkinds[bond.multiplicity])
        
        mol = mol.GetMol()
        Chem.SanitizeMol(mol)
        Chem.Kekulize(mol, clearAromaticFlags=True)

        print(Chem.MolToSmiles(mol))
        self.wid.addMol(mol)





if __name__ == '__main__':
    app = QApplication(sys.argv)
    app.setStyle('Fusion')
    ex = Wins()
    ex.show()
    app.exec_()