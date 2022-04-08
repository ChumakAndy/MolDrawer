import sys
from PyQt5.QtWidgets import QApplication, QMainWindow, QGraphicsScene, QGraphicsView, QGraphicsLineItem, QGraphicsObject
from PyQt5.QtGui import QPen, QColor, QPainter
from PyQt5.QtCore import QRectF, QPointF, QLine
from PyQt5.QtCore import Qt
import numpy as np
from math import pi
import time
from rdkit import Chem
from rdkit.Chem import AllChem

from .draver import Atom, Bound, History
from .vsp import *


class MyScene(QGraphicsScene):

    def __init__(self):
        super().__init__()
        self.setSceneRect(-5000, -5000, 10000.00, 10000.00)
        

class Grf(QGraphicsView):

    def __init__(self):
        super().__init__()

        self.sc = MyScene()
        self.setScene (self.sc)
        self.setRenderHints(QPainter.Antialiasing | QPainter.HighQualityAntialiasing)
        self.scale(0.06, 0.06)
        self.setMouseTracking(True)

        self.atomList = set()
        self.boundSet = set()
        self.history = History(self)

        self.xStart = None
        self.yStart = None
        self.xFinish = None
        self.yFinish = None
        self.isClicked = False
        self.mousePos = ()

        self.itemSelect = None
        self.atomSelect = None
        self._clickstart = None
        self.onAtom_click = None
        self.atom_text = ''
        self.atom_change = None
        self.mode = '1'
        self.ErrorList = []
    
    def findErrors(self):
        points = []

        bkinds = {1: Chem.BondType.SINGLE,
                  2: Chem.BondType.DOUBLE,
                  3: Chem.BondType.TRIPLE}
        mol = Chem.RWMol()
        list_atoms = list(self.atomList)
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
        for bond in self.boundSet:
            at1 = bond.atoms[0].IDX
            at2 = bond.atoms[1].IDX
            mol.AddBond(at1, at2, bkinds[bond.multiplicity])
        
        mol = mol.GetMol()

        problems = Chem.DetectChemistryProblems(mol)
        for problem in problems:
            point = list_atoms[problem.GetAtomIdx()]
            points.append(point)
        return points

    def addPolyg(self, item, N, benz=False):
        itemAt = None
        n = N
        r = Width/(2*np.sin(pi/N))
        atoms =[]
        if not item:
            x = self.xStart
            y = self.yStart
            alf = pi/2
        elif item.itemType == 'atom':
            atom = item
            if atom.boundList:
                if len(atom.boundList) == 1:
                    line = atom.boundList[0].line
                else:
                    line = dLine(atom, poly=True)
            else:
                line = (*atom.point, atom.point[0] + Width, atom.point[1])
            if atom.point == line[0:2]:
                line = [*line[2:4], *line[0:2]]
            x , y , cosa = getDot(line, r)
            alf = pi - np.arccos(cosa)
            if line[3] > line[1]:
                alf = - alf
        elif item.itemType == 'bound':
            bound = item
            line = bound.line
            left, right = bound.findN()
            if left < right:
                line = (*line[2:4], *line[0:2])
            v = (line[2]- line[0], line[3]- line[1])
            wi = (v[0]**2 + v[1]**2)**0.5
            r = wi/(2*np.sin(pi/N))
            q, e, cosa = getDot(line, r)
            al = np.arccos(cosa)
            if line[3] > line[1]:
                al = - al
            alf = 3*pi/2 - al - pi/N

            re = r*np.sin(pi/2 - pi/N)
            v = (line[2]- (line[0]+line[2])/2, line[3]- (line[1]+line[3])/2)
            v2 = (-2*re*v[1]/wi, 2*re*v[0]/wi)
            x = (line[0]+line[2])/2 + v2[0]
            y = (line[1]+line[3])/2 + v2[1]

        for i in range(1, n+1):
            x1 = x + r * np.cos(2 * pi * i / n + alf)
            y1 = y + r * np.sin(2 * pi * i / n + alf)
            point = (int(x1), int(y1))
            atoms.append(point)

        if benz:
            itemAt = item
            if not item:
                itemAt = True
        self.makeMolFromCoords(atoms, itemAt=itemAt)

    def makeMolFromCoords(self, points, itemAt=None):
        atoms =[]
        for point in points:
            item = self.findAtom(point)
            if item:
                atom = item
            elif not item:
                atom = Atom(self, point)
            atoms.append(atom)
        bounds = list(self.boundSet)
        newBounds = []

        for i in range(len(points)-1):
            nb = True
            line = (*atoms[i].point, *atoms[i + 1].point)
            for bd in bounds:
                if line == bd.line or (*line[2:4], *line[0:2]) == bd.line:
                    newBounds.append(bd)
                    nb = False
            if nb:
                bound = Bound(self, line)
                atoms[i].boundList.append(bound)
                atoms[i + 1].boundList.append(bound)
                bound.atoms.append(atoms[i])
                bound.atoms.append(atoms[i + 1])
                newBounds.append(bound)
        nb = True
        line = (*atoms[0].point, *atoms[-1].point)
        for bd in bounds:
            if line == bd.line or (*line[2:4], *line[0:2]) == bd.line:
                newBounds.append(bd)
                nb = False
        if nb:
            bound = Bound(self, line)
            atoms[0].boundList.append(bound)
            atoms[-1].boundList.append(bound)
            bound.atoms.append(atoms[0])
            bound.atoms.append(atoms[-1])
            newBounds.append(bound)
        if itemAt == True:
            for i in range(3):
                newBounds[i*2 + 1].multiplicity = 2
        elif itemAt and itemAt.itemType == 'bound':
            if itemAt.multiplicity == 1:
                check = False
                interBounds = itemAt.atoms[0].boundList + itemAt.atoms[1].boundList
                for bd in interBounds:
                    if bd.multiplicity != 1:
                        check = True
                if not check:
                    for i in range(3):
                        newBounds[i*2].multiplicity = 2
                if check:
                    for i in range(3):
                        newBounds[i*2 + 1].multiplicity = 2
                    itemAt.multiplicity = 1
            else:
                for i in range(3):
                    newBounds[i*2 + 1].multiplicity = 2
        elif itemAt and itemAt.itemType == 'atom':
            for i in range(3):
                newBounds[i*2 + 1].multiplicity = 2
        self.upd()

    def addMol(self, mol):
        if mol.GetNumAtoms() == 1:
            mol = Chem.AddHs(mol)
        Chem.Kekulize(mol, clearAromaticFlags=True)
        AllChem.Compute2DCoords(mol)
        bkinds = {Chem.BondType.SINGLE: 1,
                  Chem.BondType.DOUBLE: 2,
                  Chem.BondType.TRIPLE: 3}
        c = mol.GetConformers()[0]
        coordinates = Width*c.GetPositions()/1.5
        atomN = list(range(c.GetNumAtoms()))
        self.scene().clear()
        self.atomList = set()
        self.boundSet = set()
        atomL = []
        for i, k in zip(atomN, mol.GetAtoms()):
            point = (int(coordinates[i][0]), int(coordinates[i][1]))
            atom = Atom(self, point)
            atom.kind = k.GetSymbol()
            atom.charge = k.GetFormalCharge()
            atomL.append(atom)
            self.atomList.add(atom)

        for bound in mol.GetBonds():
            a1 = atomL[bound.GetBeginAtomIdx()]
            a2 = atomL[bound.GetEndAtomIdx()]
            line = (*a1.point, *a2.point)
            bo = Bound(self, line)
            bo.atoms.append(a1)
            bo.atoms.append(a2)
            a1.boundList.append(bo)
            a2.boundList.append(bo)
            bo.multiplicity = bkinds[bound.GetBondType()]
            self.boundSet.add(bo)
        self.history.update()
        self.upd()

    def deleteAtom(self, atom):
        bdl = atom.boundList[:]
        for bound in bdl:
            self.deleteBound(bound)
        self.atomList.remove(atom)
        self.history.update()
        self.upd()

    def deleteBound(self, bound):
        for atom in bound.atoms:
            atom.boundList.remove(bound)
        bound.atoms = []
        self.boundSet.remove(bound)
        self.history.update()
        self.upd()


    def findAtom(self, point):
        items = self.items(self.mapFromScene(*point))
        if items:
            for item in items:
                if item.itemType == 'atom':
                    return item
        else:
            return None

    def lineFromAtom(self, atom):
        if len(atom.boundList) == 0:
            line = (*atom.point, atom.point[0] + 500, atom.point[1])
            return line
        elif len(atom.boundList) == 1:
            line = sLine(atom)
            return line
        elif len(atom.boundList) == 2:
            line = dLine(atom)
            return line
        else:
            line = dLine(atom)
            return line

    def mousePressEvent(self, event):
        self.xStart = self.mapToScene(event.pos()).x()
        self.yStart = self.mapToScene(event.pos()).y()
        if event.button() == Qt.RightButton:
            item = self.itemAt(event.pos())
            if item:
                item = self.findAtom((self.xStart, self.yStart))
                if item:
                    self.deleteAtom(item)
                else:
                    item = self.itemAt(event.pos())
                    self.deleteBound(item)
        elif event.button() == Qt.LeftButton and self.mode == 'drag':
            self._clickstart = time.time()
            item = self.findAtom((self.xStart, self.yStart))
            if item:
                self.atomSelect = item
            self.isClicked = True
        elif event.button() == Qt.LeftButton:
            item = self.itemAt(event.pos())
            self._clickstart = time.time()
            if not item:
                if self.mode == '1':
                    atom = Atom(self, (self.xStart, self.yStart))
                    bound = Bound(self, (*atom.point, self.xStart + 500, self.yStart))
                    bound.atoms.append(atom)
                    atom.boundList.append(bound)
                    self.itemSelect = bound
                elif self.mode == 'benz':
                    self.addPolyg(item, 6, benz=True)
                    self._clickstart = None
                else:
                    self.addPolyg(item, int(self.mode))
                    self._clickstart = None
            elif item:
                newitem = self.findAtom((self.xStart, self.yStart))
                if newitem:
                    if self.mode == '1':
                        line = self.lineFromAtom(newitem)
                        bound = Bound(self, line)
                        bound.atoms.append(newitem)
                        newitem.boundList.append(bound)
                        self.itemSelect = bound
                    elif self.mode == 'benz':
                        self.addPolyg(item, 6, benz=True)
                        self._clickstart = None
                    else:
                        self.addPolyg(newitem, int(self.mode))
                        self._clickstart = None
                elif item.itemType == 'bound':
                    if self.mode == '1':
                        self.itemSelect = item
                        self._clickstart = None
                    elif self.mode == 'benz':
                        self.addPolyg(item, 6, benz=True)
                        self._clickstart = None
                    else:
                        self.addPolyg(item, int(self.mode))
                        self._clickstart = None
            self.isClicked = True


    def mouseMoveEvent(self, event):
        self.mousePos = (self.mapToScene(event.pos()).x(), self.mapToScene(event.pos()).y())
        if self.isClicked:
            if self._clickstart is not None:
                if time.time() - self._clickstart > 0.2:
                    if self.atomSelect:
                        self.atomSelect.point = self.mousePos
                    elif self.itemSelect:
                        if self.itemSelect.itemType == 'bound':
                            self.itemSelect.line = (*self.itemSelect.atoms[0].point, *self.mousePos)
                    self.upd()

    def mouseReleaseEvent(self, event):
        self.xFinish = self.mapToScene(event.pos()).x()
        self.yFinish = self.mapToScene(event.pos()).y()
        item = self.findAtom((self.xFinish, self.yFinish))
        if self.itemSelect:
            if self._clickstart is None:
                self.itemSelect.change_multiplicity()
                self.itemSelect.update()
            elif self._clickstart is not None:
                time_ = time.time() - self._clickstart
                if not item or time_ < 0.2 or self.itemSelect.atoms[0] is item:
                    point = tuple(self.itemSelect.line[2: 4])
                    atomAt = self.findAtom(point)
                    if atomAt:
                        for bd in atomAt.boundList:
                            if self.itemSelect.atoms[0] in bd.atoms:
                                self.deleteBound(self.itemSelect)
                                self.itemSelect = None
                                self.isClicked = False
                                return
                        if self.itemSelect:
                            self.itemSelect.atoms.append(atomAt)
                            atomAt.boundList.append(self.itemSelect)
                            self.itemSelect.update()
                            self.itemSelect = None
                    elif not atomAt:
                        atom = Atom(self, point)
                        self.itemSelect.atoms.append(atom)
                        atom.boundList.append(self.itemSelect)
                        self.itemSelect.update()
                        self.itemSelect = None
                elif item:
                    for bd in item.boundList:
                        if self.itemSelect.atoms[0] in bd.atoms:
                            self.deleteBound(self.itemSelect)
                            self.itemSelect = None
                            self.isClicked = False
                            return
                    if self.itemSelect:
                        self.itemSelect.atoms.append(item)
                        item.boundList.append(self.itemSelect)
                        self.itemSelect.update()
                        self.itemSelect = None
        elif self.atomSelect:
            self.atomSelect.point = self.mousePos
            self.atomSelect = None
        self.history.update()
        self.upd()
        self.ErrorList = self.findErrors()


        self.itemSelect = None
        self.isClicked = False
        self.xStart = None
        self.yStart = None
        self.xFinish = None
        self.yFinish = None
    
    def keyPressEvent(self, event):
        alphabet = 'abcdefghijklmnopqrstuvwxyz'
        modifiers = QApplication.keyboardModifiers()
        atom = self.findAtom(self.mousePos)
        if modifiers == Qt.ControlModifier:
            if event.key() == Qt.Key_Z:
                self.history.undo()
            elif event.key() == Qt.Key_Y:
                self.history.redo()
        elif atom:
            if not self.onAtom_click or time.time() - self.onAtom_click > 1:
                if event.text() and event.text() in alphabet:
                    self.atom_change = atom
                    self.atom_text = ''
                    self.onAtom_click = time.time()
                    self.atom_text += event.text()
                    atom.kind = self.atom_text.upper()
                    atom.update()
                elif event.text() and event.text() in '0-=+':
                    if event.text() in '=+':
                        atom.charge = +1
                    elif event.text() == '-':
                        atom.charge = -1
                    elif event.text() == '0':
                        atom.charge = 0
                    atom.update()
            elif time.time() - self.onAtom_click < 1:
                if atom == self.atom_change:
                    if event.text() in alphabet:
                        self.atom_text = self.atom_change.kind + event.text().lower()
                        self.atom_change.kind = self.atom_text
                        self.atom_change.update()
                        self.onAtom_click = None
                        self.atom_change = None
                elif atom != self.atom_change:
                    if event.text() in alphabet:
                        self.atom_change = atom
                        self.atom_text = ''
                        self.onAtom_click = time.time()
                        self.atom_text += event.text()
                        atom.kind = self.atom_text.upper()
                        atom.update()
            self.history.update()
            self.upd()
            self.ErrorList = self.findErrors()

    def upd(self):
        for item in self.scene().items():
            self.scene().removeItem(item)

        for bound in self.boundSet:
            self.scene().addItem(bound)
        for ato in self.atomList:
            self.scene().addItem(ato)
        for item in self.scene().items():
            item.update()


