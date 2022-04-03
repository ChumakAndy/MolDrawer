import numpy as np
from math import pi, acos

from PyQt5.QtCore import QRectF, QPoint
from PyQt5.QtGui import QPen, QColor
from PyQt5.QtCore import Qt

Redpen = QPen(QColor(255,0,0))
Redpen.setWidth(10)

Whitepen = QPen(QColor(255,255,255))


Width = 500
COLORS = {
    'O' : QPen(Qt.red),
    'S' : QPen(Qt.yellow),
    'N' : QPen(Qt.blue),
    'S' : QPen(Qt.yellow),
    'Cl' : QPen(Qt.green),
    'Br' : QPen(Qt.darkBlue),
    'F' : QPen(Qt.darkGreen),
    'I' : QPen(Qt.cyan),
    'P' : QPen(Qt.darkCyan)
}
def sLine(atom, alfa = 2*pi/3):
    line = atom.boundList[0].line
    point1 = atom.point
    if atom.point == line[2:4]:
        line = [*line[2:4], *line[0:2]]

    x1 = line[0]
    y1 = line[1]
    x2 = line[2]
    y2 = line[3]
    v = (x2 - x1, y2 - y1)
    if v[0] and v[1]:
        a = (v[1]**2)/(v[0]**2) + 1
        b = - ((2*v[1])/(v[0]**2)) *(Width**2)*np.cos(alfa)
        c = (((Width**2)*np.cos(alfa))**2)/(v[0]**2) - (Width**2)

        d = (b**2) - (4*a*c)
        l = 2 * a
        y = (-b - (d**(1/2))) / l
        x = ((Width**2)*np.cos(alfa) - v[1]*y)/v[0]
    elif v[0] == 0:
        y = (Width**2)*np.cos(alfa)/v[1]
        x = (Width**2 - y**2)**0.5
    elif v[1] == 0:
        x = (Width**2)*np.cos(alfa)/v[0]
        y = (Width**2 - x**2)**0.5

    vi = (x, y)
    point2 = (int(point1[0] + vi[0]), int(point1[1] + vi[1]))
    line = (*point1, *point2)
    return line


def dLine(atom):
    point1 = atom.point
    v = (0, 0)
    for bound in atom.boundList:
        line = bound.line
        if atom.point == line[2:4]:
            line = [*line[2:4], *line[0:2]]
        vec = (line[2] - line[0], line[3] - line[1])
        v = (v[0] + vec[0], v[1] + vec[1])

    cosa = v[0]/(v[0]**2 + v[1]**2)**(0.5)
    cosb = v[1]/(v[0]**2 + v[1]**2)**(0.5)
    vi = (-cosa*Width, -cosb*Width)

    if len(atom.boundList) % 2 != 0:
        vi = (cosa*Width, cosb*Width)
    point2 = (int(point1[0] + vi[0]), int(point1[1] + vi[1]))
    line = (*point1, *point2)
    return line


def makeLine(line, w=Width):
    x1 = line[0]
    y1 = line[1]
    x2 = line[2]
    y2 = line[3]
    v = (x2 - x1, y2 - y1)
    cosa = v[0]/(v[0]**2 + v[1]**2)**(0.5)
    cosb = v[1]/(v[0]**2 + v[1]**2)**(0.5)
    vi = (cosa*w, cosb*w)
    point1 = (line[0], line[1])
    point2 = (line[0] + vi[0], line[1] + vi[1])
    line = (*point1, *point2)
    return line


def makePline(line):
    x1 = line[0]
    y1 = line[1]
    x2 = line[2]
    y2 = line[3]
    v = (x2 - x1, y2 - y1)

    v2 = (-v[1]/12, v[0]/12)
    point1 = (x1 + v2[0], y1 + v2[1])
    p1 = QPoint(*point1)
    point2 = (x1 - v2[0], y1 - v2[1])
    p2 = QPoint(*point2)

    point3 = (point2[0] + v[0], point2[1] + v[1])
    p3 = QPoint(*point3)
    point4 = (point1[0] + v[0], point1[1] + v[1])
    p4 = QPoint(*point4)
    pol = [p2, p3, p4, p1]
    return pol

def getDot(line, w):
    x1 = line[0]
    y1 = line[1]
    x2 = line[2]
    y2 = line[3]
    v = (x2 - x1, y2 - y1)
    cosa = v[0]/(v[0]**2 + v[1]**2)**(0.5)
    cosb = v[1]/(v[0]**2 + v[1]**2)**(0.5)
    vi = (cosa*w, cosb*w)
    point2 = (line[2] + vi[0], line[3] + vi[1])
    return (*point2, cosa)
