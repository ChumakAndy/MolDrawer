from PyQt5.QtGui import QPen, QColor, QPainter, QPolygonF, QPainterPath, QFont, QBrush
from PyQt5.QtCore import QRectF, QPoint, QLine
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QGraphicsPathItem, QGraphicsItem
from rdkit import Chem
from rdkit.Chem import AllChem

from .vsp import *


class Atom(QGraphicsPathItem):

    def __init__(self, view, point, parent = None):
        super().__init__(parent)

        self.itemType = 'atom'
        self.kind = 'C'
        self.view = view
        self.point = point

        self.boundList = []

        self.view.scene().addItem(self)
        self.view.atomList.add(self)

    def paint(self, painter, options, widget):
        painter.setPen(Redpen)
        painter.setOpacity(1.0)

        if not self.boundList or self.kind != 'C':

            br = QBrush(Qt.white)
            br.setStyle(Qt.SolidPattern)
            path = QPainterPath()
            point = QPoint(*self.point)
            path.addEllipse(point, 70, 70)

            painter.fillPath(path, br)
            font = QFont()
            font.setPixelSize(125)
            font.setWeight(100)
            painter.setFont(font)
            painter.drawText(self.boundingRect(), Qt.AlignCenter, self.kind)
        self.shape()
    
    def shape(self):
        path = QPainterPath()
        point = QPoint(*self.point)
        path.addEllipse(point, 100, 100)
        return path

    def boundingRect(self):
        x = self.point[0] - 100
        y = self.point[1] - 100
        return QRectF(x, y, 200, 200)


class Bound(QGraphicsPathItem):
    
    def __init__(self, view, line, parent = None):
        super().__init__(parent)
        self.setFlag(QGraphicsPathItem.ItemIsMovable)

        self.itemType = 'bound'

        self.pen = QPen(QColor(0,0,0))
        self.pen.setWidth(30)
        self.pen.setCapStyle(Qt.RoundCap)

        self.view = view

        self.line = line
        self.atoms = []
        self.view.scene().addItem(self)
        self.view.boundSet.add(self)
        self.multiplicity = 1

    def change_multiplicity(self):
        self.multiplicity += 1
        if self.multiplicity > 3:
            self.multiplicity = 1

    def paint(self, painter, options, widget):
        if len(self.atoms) == 2:
            self.line = (*self.atoms[0].point, *self.atoms[1].point)
        else:
            try:
                line = makeLine(self.line)
                self.line = line
            except:
                pass
        poly = makePline(self.line)
        
        lines = {
            1 : [QLine(*self.line)],
            2 : [QLine(*poly[0:2]), QLine(*poly[2:4])],
            3 : [QLine(*poly[0:2]), QLine(*poly[2:4]), QLine(*self.line)]
        }

        painter.setPen(self.pen)
        painter.drawLines(*lines[self.multiplicity])
        self.shape()

    def shape(self):
        poly = makePline(self.line)
        poly = QPolygonF(poly)
        path = QPainterPath()
        path.addPolygon(poly)
        return path

    def boundingRect(self):
        w = int(Width*2)
        return QRectF(-2500, -2500, 5000, 5000)


class History:

    def __init__(self, view):
        self.view = view
        self.data  = []
        self.pos = -1
        self.update()

    def update(self):
        mol = Molecul(self.view.atomList, self.view.boundSet)

        self.data = self.data[:self.pos+1]
        self.data.append(mol)
        self.pos = len(self.data) - 1

    def restore(self, mol):
        self.view.scene().clear()
        self.view.atomList = set()
        self.view.boundSet = set()
        for atom in mol.atoms:
            ato = Atom(self.view, atom['point'])
            ato.kind = atom['kind']
            self.view.atomList.add(ato)
        for bound in mol.bounds:
            bd = Bound(self.view, bound['line'])
            bd.multiplicity = bound['multiplicity']
            self.view.boundSet.add(bd)
            point1 = bd.line[0:2]
            point2 = bd.line[2:4]
            for a in self.view.atomList:
                if point1 == ato.point:
                    bd.atoms.append(a)
                    a.boundList.append(bd)
                if point2 == ato.point:
                    bd.atoms.append(a)
                    a.boundList.append(bd)


    def undo(self):
        if self.pos == 0:
            return None
        self.pos -= 1
        self.restore(self.data[self.pos])
        self.view.upd()

    def redo(self):
        if self.pos == len(self.data)-1:
            return None
        self.pos += 1
        self.restore(self.data[self.pos])
        self.view.upd()

class Molecul:

    def __init__(self, atoms, bounds):

        self.atoms = []
        self.bounds = []
        for atom in atoms:
            a = {}
            a['point'] = atom.point[:]
            a['kind'] = atom.kind
            a['bounds'] = []
            self.atoms.append(a)

        for bound in bounds:
            b = {}
            b['multiplicity'] = bound.multiplicity
            b['line'] = bound.line[:]
            b['atoms'] = []
            self.bounds.append(b)


