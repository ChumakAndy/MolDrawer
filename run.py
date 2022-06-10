from PyQt5.QtWidgets import QApplication, QMainWindow, QDockWidget, QAction, QActionGroup, QErrorMessage, QMessageBox

from rdkit import Chem
import traceback

if __name__ == '__main__':
    from engine import *
else:
    from .engine import *



class MolDrawer(QDockWidget):

    def __init__(self, caller, smiles , parent = None):
        super().__init__(parent)

        self.caller = caller

        self.window = QMainWindow()

        self.view = GrV()
        toolbar = self.window.addToolBar('Experiment')
        toolbar.setContextMenuPolicy(Qt.PreventContextMenu)
        toolbar.setMovable(False)
        upd = QAction(get_icon('update'), 'Update',self)
        splin = QAction(get_icon('sanitize'), 'Sanitize',self)
        un = QAction(get_icon('undo'), 'Undo',self)
        re = QAction(get_icon('redo'), 'Redo',self)
        cl = QAction(get_icon('clear'), 'Clean',self)

        toolbar.addAction(upd)
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

        upd.triggered.connect(self.getMol)
        splin.triggered.connect(self.sanitize)
        re.triggered.connect(self.redo)
        un.triggered.connect(self.undo)
        cl.triggered.connect(self.clean_)
        self.window.setCentralWidget(self.view)
        self.setWidget(self.window)
        self.setWindowIcon(get_icon('flask'))
        self.setWindowTitle("MolDrawer")
        try:
            self.view.addMol(Chem.MolFromSmiles(smiles))
        except:
            self.error_d = QErrorMessage()
            self.error_d.showMessage(traceback.format_exc())

    def undo(self):
        self.view.history.undo()
        self.view.ErrorList = self.wid.findErrors()

    def redo(self):
        self.view.history.redo()
        self.view.ErrorList = self.wid.findErrors()

    def clean_(self):
        self.view.scene().clear()
        self.view.history.data  = []
        self.view.history.pos = -1
        self.view.history.update()

    def getMol(self):
        if self.view.ErrorList:
            self.error_d = QErrorMessage()
            self.error_d.showMessage('Bad structure!')
        else:
            mol = self.view.makeMol()
        self.caller(mol)
        self.close()

    def sanitize(self):
        try:
            mol = self.view.makeMol()
            self.view.addMol(mol)
        except:
            self.error_d = QErrorMessage()
            self.error_d.showMessage('Bad structure!')

    def callBot(self):
        sender = self.sender().text()
        self.view.mode = sender



if __name__ == '__main__':
    app = QApplication(sys.argv)
    app.setStyle('Fusion')
    ex = MolDrawer(print, '')
    ex.show()
    app.exec_()
