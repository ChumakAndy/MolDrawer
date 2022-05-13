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
        self.setFlag(QGraphicsPathItem.ItemIsMovable)

        self.itemType = 'atom'
        self.kind = 'C'
        self.charge = 0
        self.view = view
        self.point = point

        self.boundList = []

        self.view.scene().addItem(self)
        self.view.atomList.add(self)
        self.setZValue(1)


    def paint(self, painter, options, widget):
        list_atoms = ['O', 'S', 'N', 'Cl', 'Br', 'F', 'I', 'P']
        if self.kind in list_atoms:
            painter.setPen(COLORS[self.kind])
        painter.setOpacity(1.0)

        if not self.boundList or self.kind != 'C' or self.charge:

            br = QBrush(Qt.white)
            br.setStyle(Qt.SolidPattern)
            path = QPainterPath()
            point = QPoint(*self.point)
            path.addEllipse(point, 150, 150)
 
            painter.fillPath(path, br)
            font = QFont()
            font.setPixelSize(150)
            font.setWeight(100)
            painter.setFont(font)
            painter.drawText(self.boundingRect(), Qt.AlignCenter, self.kind)
        if self.charge:
            x = self.point[0] + 75
            y = self.point[1] - 200
            cR =  QRectF(x, y, 200, 200)
            br = QBrush(Qt.white)
            br.setStyle(Qt.SolidPattern)
            painter.fillRect(cR, br)
            if self.charge >0:
                text = '+'
            if self.charge <0:
                text = '-'
            if abs(self.charge) != 1:
                text = str(abs(self.charge)) + text
            painter.drawText(cR, Qt.AlignCenter, text)

        if self in self.view.ErrorList:
            pen = QPen(Qt.red)
            pen.setWidth(15)
            painter.setPen(pen)
            x = self.point[0] - 100
            y = self.point[1] - 100
            painter.drawEllipse(x, y, 200, 200)

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
        self.pen.setWidth(Width/15)
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
            3 : [makeShortLine(poly[0:2]), makeShortLine(poly[2:4]), QLine(*self.line)]
        }
        painter.setPen(self.pen)
        if self.multiplicity == 1 or self.multiplicity == 3:
            painter.drawLines(*lines[self.multiplicity])
        elif self.multiplicity == 2:
            left, right = self.findN()
            if left > right:
                painter.drawLines(QLine(*self.line), makeShortLine(poly[0:2]))
            else:
                painter.drawLines(QLine(*self.line), makeShortLine(poly[2:4]))

        self.shape()

    def shape(self):
        poly = makePline(self.line)
        poly = QPolygonF(poly)
        path = QPainterPath()
        path.addPolygon(poly)
        return path

    def boundingRect(self):
        w = int(Width*2)
        return QRectF(-5000, -5000, 10000, 10000)
    
    def NPath(self, beta):
        path = QPainterPath()
        q,e, cosa = getDot(self.line, Width)
        alf = 180*acos(cosa)/pi
        if self.line[3] > self.line[1]:
            alf = -alf

        path.moveTo(*self.line[0:2])
        qr = QRectF(self.line[0:2][0] - Width/2, self.line[0:2][1] - Width/2, Width, Width)
        path.arcTo(qr, alf, beta)
        path.moveTo(*self.line[2:4])
        qr = QRectF(self.line[2:4][0] - Width/2, self.line[2:4][1] - Width/2, Width, Width)
        path.arcTo(qr, alf, beta)

        pa = QPainterPath()
        point = QPoint(*self.line[0:2])
        pa.addEllipse(point, 150, 150)
        path -= pa

        pa = QPainterPath()
        point = QPoint(*self.line[2:4])
        pa.addEllipse(point, 150, 150)
        path -= pa

        return path

    def findN(self):
        path = self.NPath(180)
        l = len(self.view.items(self.view.mapFromScene(path)))
        path = self.NPath(-180)
        r = len(self.view.items(self.view.mapFromScene(path)))
        return l, r


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
            ato.charge = atom['charge']
            self.view.atomList.add(ato)
        for bound in mol.bounds:
            bd = Bound(self.view, bound['line'])
            bd.multiplicity = bound['multiplicity']
            self.view.boundSet.add(bd)
            point1 = bd.line[0:2]
            point2 = bd.line[2:4]
            for a in self.view.atomList:
                if point1 == a.point:
                    bd.atoms.append(a)
                    a.boundList.append(bd)
                if point2 == a.point:
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
            self.addAtom(atom)

        for bound in bounds:
            self.addBound(bound)

    def addAtom(self, atom):
        a = {}
        a['point'] = atom.point[:]
        a['kind'] = atom.kind
        a['bounds'] = []
        a['charge'] = atom.charge
        self.atoms.append(a)
    
    def addBound(self, bound):
        b = {}
        b['multiplicity'] = bound.multiplicity
        b['line'] = bound.line[:]
        b['atoms'] = []
        self.bounds.append(b)
