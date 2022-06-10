import sys
from PyQt5.QtWidgets import QApplication, QMainWindow, QGraphicsScene, QGraphicsView, QGraphicsLineItem, QGraphicsObject, QErrorMessage
from PyQt5.QtGui import QPen, QColor, QPainter
from PyQt5.QtCore import QRectF, QPointF, QLine
from PyQt5.QtCore import Qt
import numpy as np
from math import pi
import time
from rdkit import Chem
from rdkit.Chem import AllChem
import traceback

from .drawer import Atom, Bound, History
from .drFuncs import *
from .settings import *


class MyScene(QGraphicsScene):

    def __init__(self):
        super().__init__()
        l = BOUND_LENGHT
        self.setSceneRect(-7*l, -7*l, 14*l, 14*l)


class GrV(QGraphicsView):

    def __init__(self):
        super().__init__()

        self.sc = MyScene()
        self.setScene (self.sc)

        self.setFixedSize(750, 750)

        self.setRenderHints(QPainter.Antialiasing | QPainter.HighQualityAntialiasing)
        self.scale(0.071, 0.071)
        self.setMouseTracking(True)

        self.history = History(self)

        self.mousePos = None
        self.xStart = None
        self.yStart = None
        self.xFinish = None
        self.yFinish = None
        self._clickstart = None
        self.BoundSelected = None
        self.AtomSelected = None
        self.isClicked = False
        self.onAtom_click = None
        self.atom_change = None
        self.mode = '1'
        self.ErrorList = []

    def deleteAtom(self, atom):
        bounds = atom.boundlist()
        for bound in bounds:
            self.deleteBound(bound)
        self.scene().removeItem(atom)

    def deleteBound(self, bound):
        bound.atoms = []
        self.scene().removeItem(bound)

    def findAtom(self, point):
        if not point:
            return None
        items = self.items(self.mapFromScene(*point))
        if items:
            for item in items:
                if item.itemType == 'atom':
                    return item
        return None

    def addPolygon(self, item, N, benz=False):
        itemAt = None
        points = makePointsOfPolygon(self, item, N)
        if benz:
            itemAt = item
            if not item:
                itemAt = True
        self.makeMolFromCoords(points, itemAt=itemAt)

    def mousePressEvent(self, event):
        self.xStart = self.mapToScene(event.pos()).x()
        self.yStart = self.mapToScene(event.pos()).y()
        item = self.itemAt(event.pos())

        if event.button() == Qt.RightButton:
            if item and self.findAtom((self.xStart, self.yStart)):
                item = self.findAtom((self.xStart, self.yStart))
                self.deleteAtom(item)
            elif item:
                self.deleteBound(item)

        elif event.button() == Qt.LeftButton and self.mode == '1':
            self._clickstart = time.time()
            if not item:
                atom = Atom(self, (self.xStart, self.yStart))
                bound = Bound(self, (*atom.point, self.xStart + BOUND_LENGHT, self.yStart))
                bound.atoms.append(atom)
                self.BoundSelected = bound
            elif item and self.findAtom((self.xStart, self.yStart)):
                atom = self.findAtom((self.xStart, self.yStart))
                line = lineFromAtom(atom)
                bound = Bound(self, line)
                bound.atoms.append(atom)
                self.BoundSelected = bound
            elif item.itemType == 'bound':
                self.BoundSelected = item
                self._clickstart = None

        elif event.button() == Qt.LeftButton and self.mode == 'benz':
            self.addPolygon(item, 6, benz=True)
            self._clickstart = None

        elif event.button() == Qt.LeftButton and self.mode in ['3', '4', '5', '6']:
            self.addPolygon(item, int(self.mode))
            self._clickstart = None
        
        elif event.button() == Qt.LeftButton and self.mode == 'drag':
            self._clickstart = time.time()
            item = self.findAtom((self.xStart, self.yStart))
            if item:
                self.AtomSelected = item

        self.isClicked = True

    def mouseMoveEvent(self, event):
        self.mousePos = (self.mapToScene(event.pos()).x(), self.mapToScene(event.pos()).y())
        if self.isClicked and self._clickstart is not None:
            if time.time() - self._clickstart > 0.2:
                if self.BoundSelected:
                    self.BoundSelected.line = (*self.BoundSelected.atoms[0].point, *self.mousePos)
                elif self.AtomSelected:
                        self.AtomSelected.point = self.mousePos
        self.upd()

    def mouseReleaseEvent(self, event):
        self.xFinish = self.mapToScene(event.pos()).x()
        self.yFinish = self.mapToScene(event.pos()).y()
        AtomAtFinish = self.findAtom((self.xFinish, self.yFinish))
        if self.BoundSelected and self._clickstart is None:
            self.BoundSelected.change_multiplicity()
        elif self.BoundSelected and self._clickstart is not None:
            time_ = time.time() - self._clickstart
            if not AtomAtFinish or time_ < 0.2 or self.BoundSelected.atoms[0] is AtomAtFinish:
                point = tuple(self.BoundSelected.line[2: 4])
                atomAt = self.findAtom(point)
                if atomAt:
                    for bd in atomAt.boundlist():
                        if self.BoundSelected.atoms[0] in bd.atoms:
                            self.deleteBound(self.BoundSelected)
                            self.BoundSelected = None
                            self.isClicked = False
                            return
                    if self.BoundSelected:
                        self.BoundSelected.atoms.append(atomAt)
                elif not atomAt:
                    atom = Atom(self, point)
                    self.BoundSelected.atoms.append(atom)
            elif AtomAtFinish:
                for bound in AtomAtFinish.boundlist():
                    if self.BoundSelected.atoms[0] in bound.atoms:
                        self.deleteBound(self.BoundSelected)
                        self.BoundSelected = None
                        self.isClicked = False
                        return
                if self.BoundSelected:
                    self.BoundSelected.atoms.append(AtomAtFinish)
        elif self.AtomSelected:
            self.AtomSelected = None

        self.BoundSelected = None
        self.isClicked = False
        self.xStart = None
        self.yStart = None
        self.xFinish = None
        self.yFinish = None
        self.history.update()
        self.upd()
        try:
            self.ErrorList = self.findErrors()
        except:
            self.error_d = QErrorMessage()
            self.error_d.showMessage(traceback.format_exc())

    def keyPressEvent(self, event):
        alphabet = 'abcdefghijklmnopqrstuvwxyz'
        modifiers = QApplication.keyboardModifiers()
        atom = self.findAtom(self.mousePos)
        if modifiers == Qt.ControlModifier:
            if event.key() == Qt.Key_Z:
                self.history.undo()
                self.ErrorList = self.findErrors()
            elif event.key() == Qt.Key_Y:
                self.history.redo()
                self.ErrorList = self.findErrors()
        elif atom:
            if not self.onAtom_click or time.time() - self.onAtom_click > 1:
                if event.text() and event.text() in alphabet:
                    self.atom_change = atom
                    atom_text = ''
                    self.onAtom_click = time.time()
                    atom_text += event.text()
                    atom.kind = atom_text.upper()
                elif event.text() and event.text() in '0-=+':
                    if event.text() in '=+':
                        atom.charge += 1
                    elif event.text() == '-':
                        atom.charge -= 1
                    elif event.text() == '0':
                        atom.charge = 0
            elif time.time() - self.onAtom_click < 1:
                if atom == self.atom_change:
                    if event.text() in alphabet:
                        atom_text = self.atom_change.kind + event.text().lower()
                        self.atom_change.kind = atom_text
                        self.onAtom_click = None
                        self.atom_change = None
                elif atom != self.atom_change:
                    if event.text() in alphabet:
                        self.atom_change = atom
                        atom_text = ''
                        self.onAtom_click = time.time()
                        atom_text += event.text()
                        atom.kind = atom_text.upper()
            self.history.update()
            self.upd()
            self.ErrorList = self.findErrors()

    def upd(self):
        for item in self.scene().items():
            item.update()

    def addMol(self, mol):
        if Chem.MolToSmiles(mol) == '[H][H]':
            point1 = (-BOUND_LENGHT/2, 0)
            point2 = (BOUND_LENGHT/2, 0)
            atom1 = Atom(self, point1)
            atom1.kind = 'H'
            atom2 = Atom(self, point2)
            atom2.kind = 'H'
            bo = Bound(self, (*point1, *point2))
            bo.atoms.append(atom1)
            bo.atoms.append(atom2)
            self.history.update()
            self.upd()
            return

        Chem.Kekulize(mol, clearAromaticFlags=True)
        AllChem.Compute2DCoords(mol)
        bkinds = {Chem.BondType.SINGLE: 1,
                  Chem.BondType.DOUBLE: 2,
                  Chem.BondType.TRIPLE: 3}
        c = mol.GetConformers()[0]
        coordinates = BOUND_LENGHT*c.GetPositions()/1.5 #НЕЗАБУДЬ УБРАТЬ ЭТО!!!
        atomN = list(range(c.GetNumAtoms()))
        self.scene().clear()
        atomL = []
        for i, k in zip(atomN, mol.GetAtoms()):
            point = (int(coordinates[i][0]), int(coordinates[i][1]))
            atom = Atom(self, point)
            atom.kind = k.GetSymbol()
            atom.charge = k.GetFormalCharge()
            atomL.append(atom)

        for bound in mol.GetBonds():
            a1 = atomL[bound.GetBeginAtomIdx()]
            a2 = atomL[bound.GetEndAtomIdx()]
            line = (*a1.point, *a2.point)
            bo = Bound(self, line)
            bo.atoms.append(a1)
            bo.atoms.append(a2)
            bo.multiplicity = bkinds[bound.GetBondType()]
        self.history.update()
        self.upd()
    
    def makeMol(self):
        bkinds = {1: Chem.BondType.SINGLE,
                    2: Chem.BondType.DOUBLE,
                    3: Chem.BondType.TRIPLE}
        mol = Chem.RWMol()
        list_atoms = []
        for item in self.items():
            if item.itemType == 'atom':
                list_atoms.append(item)

        for atom in list_atoms:
            atom.IDX = list_atoms.index(atom)
            rdatom = Chem.Atom(atom.kind)
            if atom.charge != 0:
                rdatom.SetFormalCharge(atom.charge)
            mol.AddAtom(rdatom)

        list_bounds = []
        for item in self.items():
            if item.itemType == 'bound':
                list_bounds.append(item)

        for bond in list_bounds:
            at1 = bond.atoms[0].IDX
            at2 = bond.atoms[1].IDX
            mol.AddBond(at1, at2, bkinds[bond.multiplicity])

        mol = mol.GetMol()

        Chem.SanitizeMol(mol)
        Chem.Kekulize(mol, clearAromaticFlags=True)
        return mol

    def makeMolFromCoords(self, points, itemAt=None):
        atoms =[]
        for point in points:
            item = self.findAtom(point)
            if item:
                atom = item
            elif not item:
                atom = Atom(self, point)
            atoms.append(atom)
        bounds = []
        for item in self.items():
            if item.itemType == 'bound':
                bounds.append(item)
        newBounds = []
        for i in range(len(points)):
            nb = True
            line = (*atoms[i - 1].point, *atoms[i].point)
            for bd in bounds:
                if line == bd.line or (*line[2:4], *line[0:2]) == bd.line:
                    newBounds.append(bd)
                    nb = False
            if nb:
                bound = Bound(self, line)
                bound.atoms.append(atoms[i - 1])
                bound.atoms.append(atoms[i])
                newBounds.append(bound)

        if itemAt == True:
            for i in range(3):
                newBounds[i*2 + 1].multiplicity = 2
        elif itemAt and itemAt.itemType == 'bound':
            if itemAt.multiplicity == 1:
                check = False
                interBounds = itemAt.atoms[0].boundlist() + itemAt.atoms[1].boundlist()
                for bd in interBounds:
                    if bd.multiplicity != 1:
                        check = True
                if not check:
                    for i in range(3):
                        newBounds[i*2 + 1].multiplicity = 2
                if check:
                    for i in range(3):
                        newBounds[i*2].multiplicity = 2
                    itemAt.multiplicity = 1
            else:
                for i in range(3):
                    newBounds[i*2].multiplicity = 2
        elif itemAt and itemAt.itemType == 'atom':
            for i in range(3):
                newBounds[i*2].multiplicity = 2
        self.upd()
    
    def findErrors(self):
        points = []

        bkinds = {1: Chem.BondType.SINGLE,
                  2: Chem.BondType.DOUBLE,
                  3: Chem.BondType.TRIPLE}
        mol = Chem.RWMol()
        list_atoms = []
        for item in self.items():
            if item.itemType == 'atom':
                list_atoms.append(item)

        try:
            for atom in list_atoms:
                atom.IDX = list_atoms.index(atom)
                rdatom = Chem.Atom(atom.kind)
                if atom.charge != 0:
                    rdatom.SetFormalCharge(atom.charge)
                mol.AddAtom(rdatom)
        except:
            points.append(atom)
            return points

        list_bounds = []
        for item in self.items():
            if item.itemType == 'bound':
                list_bounds.append(item)
        for bond in list_bounds:
            at1 = bond.atoms[0].IDX
            at2 = bond.atoms[1].IDX
            mol.AddBond(at1, at2, bkinds[bond.multiplicity])
        mol = mol.GetMol()

        problems = Chem.DetectChemistryProblems(mol)
        for problem in problems:
            point = list_atoms[problem.GetAtomIdx()]
            points.append(point)
        return points