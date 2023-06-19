from PyQt5.QtGui import QPen, QColor, QPainter, QPolygonF, QPainterPath, QFont, QBrush
from PyQt5.QtCore import QRectF, QPoint, QLine
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QGraphicsPathItem, QGraphicsItem
from rdkit import Chem
from rdkit.Chem import AllChem
from math import pi, acos

from .drFuncs import *
from .settings import *


class Region(QGraphicsPathItem):

    def __init__(self, view, point, parent = None):
        super().__init__(parent)

        self.view = view
        self.itemType = 'region'

        self.point1 = point
        self.point2 = point
        self.start_point = None
        self.complite = False
        self.atoms = []
        self.view.scene().addItem(self)
        self.pen = QPen(Qt.blue)
        self.pen.setWidth(15)
        self.pen.setStyle(Qt.DashLine)
        self.pen.setJoinStyle(Qt.RoundJoin)


    def move(self, point):
        x = point[0] - self.start_point[0]
        y = point[1] - self.start_point[1]
        self.point1 = (self.point1[0] + x, self.point1[1] + y)
        self.point2 = (self.point2[0] + x, self.point2[1] + y)
        for atom in self.atoms:
            atom.point = (atom.point[0] + x, atom.point[1] + y)
        self.start_point = point

    def getRect(self):
        path = QPainterPath()
        w = self.point2[0] - self.point1[0]
        h = self.point2[1] - self.point1[1]
        path.addRect(*self.point1, w, h)
        return path

    def paint(self, painter, options, widget):
        w = self.point2[0] - self.point1[0]
        h = self.point2[1] - self.point1[1]
        painter.setPen(self.pen)
        painter.drawRect(*self.point1, w, h)

    def shape(self):
        path = self.getRect()
        return path

    def boundingRect(self):
        return self.view.scene().sceneRect()


class Atom(QGraphicsPathItem):

    def __init__(self, view, point, parent = None):
        super().__init__(parent)

        self.itemType = 'atom'
        self.kind = 'C'
        self.charge = 0
        self.view = view
        self.point = point

        self.view.scene().addItem(self)
        self.setZValue(1)
        self.setAcceptHoverEvents(True)

    def boundlist(self):
        Bounds = []
        path = QPainterPath()
        point = QPoint(*self.point)
        path.addEllipse(point, ATOM_RADIUS - 50, ATOM_RADIUS - 50)
        items = self.view.items(self.view.mapFromScene(path))
        # items = self.view.items(self.view.mapFromScene(*self.point))
        for item in items:
            if item.itemType == 'bound':
                if self in item.atoms:
                    Bounds.append(item)
        return Bounds

    def paint(self, painter, options, widget):
        list_atoms = ['O', 'S', 'N', 'Cl', 'Br', 'F', 'I', 'P']
        if self.kind in list_atoms:
            painter.setPen(COLORS[self.kind])
        if not self.boundlist() or self.kind != 'C' or self.charge:

            br = QBrush(Qt.white)
            br.setStyle(Qt.SolidPattern)
            path = QPainterPath()
            point = QPoint(*self.point)
            path.addEllipse(point, 150, 150)
            painter.fillPath(path, br)
            font = QFont()
            font.setPixelSize(ATOM_SYMBOL_SIZE)
            font.setWeight(3)
            painter.setFont(font)
            painter.drawText(self.boundingRect(), Qt.AlignCenter, self.kind)

        if self.charge:
            x = self.point[0] + ATOM_SYMBOL_SIZE/2
            y = self.point[1] - ATOM_SYMBOL_SIZE*1.3
            cR =  QRectF(x, y, ATOM_SYMBOL_SIZE*1.3, ATOM_SYMBOL_SIZE*1.3)
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
            painter.setPen(Redpen)
            x = self.point[0] - ATOM_RADIUS
            y = self.point[1] - ATOM_RADIUS
            painter.drawEllipse(x, y, ATOM_RADIUS*2, ATOM_RADIUS*2)

        if self.view.region and self in self.view.region.atoms:
            painter.setPen(Bluepen)
            x = self.point[0] - ATOM_RADIUS
            y = self.point[1] - ATOM_RADIUS
            painter.drawEllipse(x, y, ATOM_RADIUS*2, ATOM_RADIUS*2)

        if self in self.view.itemsIn:
            painter.setPen(Bluepen)
            x = self.point[0] - ATOM_RADIUS
            y = self.point[1] - ATOM_RADIUS
            painter.drawEllipse(x, y, ATOM_RADIUS*2, ATOM_RADIUS*2)

    def shape(self):
        path = QPainterPath()
        point = QPoint(*self.point)
        path.addEllipse(point, ATOM_RADIUS, ATOM_RADIUS)
        return path

    def boundingRect(self):
        x = self.point[0] - 100
        y = self.point[1] - 100
        return QRectF(x, y, ATOM_RADIUS + 100, ATOM_RADIUS + 100)


class Bound(QGraphicsPathItem):

    def __init__(self, view, line, parent = None):
        super().__init__(parent)
        self.setFlag(QGraphicsPathItem.ItemIsMovable)

        self.itemType = 'bound'
        self.multiplicity = 1
        self.atoms = []

        self.pen = QPen(QColor(0,0,0))
        self.pen.setWidth(BOUND_WIDTH)
        self.pen.setCapStyle(Qt.RoundCap)

        self.view = view

        self.line = line
        self.view.scene().addItem(self)

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
        frame = makeFrame(self.line)
        
        lines = {
            1 : [QLine(*self.line)],
            2 : [QLine(*frame[0:2]), QLine(*frame[2:4])],
            3 : [makeShortLine(frame[0:2]), makeShortLine(frame[2:4]), QLine(*self.line)]
        }
        painter.setPen(self.pen)
        if self.multiplicity == 1 or self.multiplicity == 3:
            painter.drawLines(*lines[self.multiplicity])
        elif self.multiplicity == 2:
            left, right = self.findN()
            if left > right:
                painter.drawLines(QLine(*self.line), makeShortLine(frame[0:2]))
            else:
                painter.drawLines(QLine(*self.line), makeShortLine(frame[2:4]))
        if self in self.view.itemsIn:
            poly = makeFrame(self.line)
            poly = QPolygonF(poly)
            path = QPainterPath()
            path.addPolygon(poly)
            br = QBrush(Qt.blue)
            br.setStyle(Qt.SolidPattern)
            painter.setOpacity(0.5)
            painter.fillPath(path, br)

        self.shape()

    def shape(self):
        poly = makeFrame(self.line)
        poly = QPolygonF(poly)
        path = QPainterPath()
        path.addPolygon(poly)
        return path

    def boundingRect(self):
        w = int(BOUND_LENGHT*2)
        return self.view.scene().sceneRect()

    def NPath(self, beta):
        lenght = BOUND_LENGHT
        path = QPainterPath()
        q,e, cosa = getCenterPolygon(self.line, lenght)
        alf = 180*acos(cosa)/pi
        if self.line[3] > self.line[1]:
            alf = -alf

        path.moveTo(*self.line[0:2])
        qr = QRectF(self.line[0:2][0] - lenght/2, self.line[0:2][1] - lenght/2, lenght, lenght)
        path.arcTo(qr, alf, beta)
        path.moveTo(*self.line[2:4])
        qr = QRectF(self.line[2:4][0] - lenght/2, self.line[2:4][1] - lenght/2, lenght, lenght)
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
        molecul = Molecul(self.view)

        self.data = self.data[:self.pos+1]
        self.data.append(molecul)
        self.pos = len(self.data) - 1

    def restore(self,molecul):
        self.view.scene().clear()

        for atom in molecul.atoms:
            ato = Atom(self.view, atom['point'])
            ato.kind = atom['kind']
            ato.charge = atom['charge']

        for bound in molecul.bounds:
            bd = Bound(self.view, bound['line'])
            bd.multiplicity = bound['multiplicity']
            point1 = bd.line[0:2]
            point2 = bd.line[2:4]
            atomlist = []
            for item in self.view.items():
                if item.itemType == 'atom':
                    atomlist.append(item)
            for a in atomlist:
                if point1 == a.point:
                    bd.atoms.append(a)
                if point2 == a.point:
                    bd.atoms.append(a)

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

    def __init__(self, view):

        self.view = view
        self.atoms = []
        self.bounds = []
        for item in view.items():
            if item.itemType == 'atom':
                self.addAtom(item)
            elif item.itemType == 'bound':
                self.addBound(item)

    def addAtom(self, atom):
        a = {}
        a['point'] = atom.point[:]
        a['kind'] = atom.kind[:]
        a['charge'] = atom.charge
        self.atoms.append(a)

    def addBound(self, bound):
        b = {}
        b['multiplicity'] = bound.multiplicity
        b['line'] = bound.line[:]
        self.bounds.append(b)