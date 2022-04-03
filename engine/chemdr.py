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
        self.setSceneRect(-2500, -2500, 5000.00, 5000.00)
        

class Grf(QGraphicsView):

    def __init__(self):
        super().__init__()

        self.sc = MyScene()
        self.setScene (self.sc)
        self.setRenderHints(QPainter.Antialiasing | QPainter.HighQualityAntialiasing)
        self.scale(0.1, 0.1)
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
        self._clickstart = None
        self.onAtom_click = None
        self.atom_text = ''
        self.atom_change = None
    
    def addPolyg(self, item, N):
        n = N
        r = Width/(2*np.sin(pi/N))
        atoms =[]
        if item.itemType == 'atom':
            atom = item
            line = atom.boundList[0].line
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
            q, e, cosa = getDot(line, r)
            alf = pi + np.arccos(cosa)
            y = line[1] + r*np.sin(pi - (pi - 2*pi/N) + np.arccos(cosa))
            x = line[0] + r*np.cos(pi - (pi - 2*pi/N) + np.arccos(cosa))
            if line[3] < line[1]:
                y = line[1] + r*np.sin(pi - (pi - 2*pi/N) - np.arccos(cosa))
                x = line[0] + r*np.cos(pi - (pi - 2*pi/N) - np.arccos(cosa))
                alf = pi - np.arccos(cosa)
        for i in range(1, n+1):
            x1 = x + r * np.cos(2 * pi * i / n + alf)
            y1 = y + r * np.sin(2 * pi * i / n + alf)
            point = (int(x1), int(y1))
            atoms.append(point)
        self.makeMolFromCoords(atoms)

    def makeMolFromCoords(self, points):
        atoms =[]
        for point in points:
            item = self.findAtom(point)
            if item:
                atom = item
            elif not item:
                atom = Atom(self, point)
            atoms.append(atom)
        bounds = list(self.boundSet)
        
        for i in range(len(points)-1):
            nb = True
            line = (*atoms[i].point, *atoms[i + 1].point)
            for bd in bounds:
                if line == bd.line or (*line[2:4], *line[0:2]) == bd.line:
                    nb = False
            if nb:
                bound = Bound(self, line)
                atoms[i].boundList.append(bound)
                atoms[i + 1].boundList.append(bound)
                bound.atoms.append(atoms[i])
                bound.atoms.append(atoms[i + 1])
        nb = True
        line = (*atoms[0].point, *atoms[-1].point)
        for bd in bounds:
            if line == bd.line or (*line[2:4], *line[0:2]) == bd.line:
                nb = False
        if nb:
            bound = Bound(self, line)
            atoms[0].boundList.append(bound)
            atoms[-1].boundList.append(bound)
            bound.atoms.append(atoms[0])
            bound.atoms.append(atoms[-1])

        self.upd()

    def addMol(self, mol):
        if mol.GetNumAtoms() == 1:
            mol = Chem.AddHs(mol)
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
            # item = self.itemAt(event.pos())
            # if item:
            #     item = self.findAtom((self.xStart, self.yStart))
            #     if item:
            #         self.deleteAtom(item)
            #     else:
            #         item = self.itemAt(event.pos())
            #         self.deleteBound(item)
            item = self.itemAt(event.pos())
            if item:
                item = self.findAtom((self.xStart, self.yStart))
                if item:
                    self.addPolyg(item, 6)
                else:
                    item = self.itemAt(event.pos())
                    self.addPolyg(item, 6)


        elif event.button() == Qt.LeftButton:
            item = self.itemAt(event.pos())
            self._clickstart = time.time()
            if not item:
                atom = Atom(self, (self.xStart, self.yStart))
                bound = Bound(self, (*atom.point, self.xStart + 500, self.yStart))
                bound.atoms.append(atom)
                atom.boundList.append(bound)
                self.itemSelect = bound
            elif item:
                newitem = self.findAtom((self.xStart, self.yStart))
                if newitem:
                    line = self.lineFromAtom(newitem)
                    bound = Bound(self, line)
                    bound.atoms.append(newitem)
                    newitem.boundList.append(bound)
                    self.itemSelect = bound
                elif item.itemType == 'bound':
                    self.itemSelect = item
                    self._clickstart = None
            self.isClicked = True
            self.upd()

    def mouseMoveEvent(self, event):
        self.mousePos = (self.mapToScene(event.pos()).x(), self.mapToScene(event.pos()).y())
        if self.isClicked:
            if self._clickstart is not None:
                if time.time() - self._clickstart > 0.2:
                    point = (self.mapToScene(event.pos()).x(), self.mapToScene(event.pos()).y())
                    if self.itemSelect.itemType == 'bound':
                        self.itemSelect.line = (*self.itemSelect.atoms[0].point,
                         self.mapToScene(event.pos()).x(),
                          self.mapToScene(event.pos()).y())
                        self.upd()

    def mouseReleaseEvent(self, event):
        self.xFinish = self.mapToScene(event.pos()).x()
        self.yFinish = self.mapToScene(event.pos()).y()
        item = self.findAtom((self.xFinish, self.yFinish))
        if self.itemSelect:
            if self.itemSelect.itemType == 'bound':
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
                        elif not atomAt:
                            atom = Atom(self, point)
                            self.itemSelect.atoms.append(atom)
                            atom.boundList.append(self.itemSelect)
                            self.itemSelect.update()
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
        self.history.update()
        self.upd()

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

    def upd(self):
        for item in self.scene().items():
            self.scene().removeItem(item)
        for bound in self.boundSet:
            self.scene().addItem(bound)
        for ato in self.atomList:
            self.scene().addItem(ato)
        for item in self.scene().items():
            item.update()


