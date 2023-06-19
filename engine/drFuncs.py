import numpy as np
from math import pi
from PyQt5.QtCore import QRectF, QPoint, QLine
from .settings import *


def makeFrame(line):
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

def makeLine(line, lenght=BOUND_LENGHT):
    x1 = line[0]
    y1 = line[1]
    x2 = line[2]
    y2 = line[3]
    v = (x2 - x1, y2 - y1)
    cosa = v[0]/(v[0]**2 + v[1]**2)**(0.5)
    cosb = v[1]/(v[0]**2 + v[1]**2)**(0.5)
    vi = (cosa*lenght, cosb*lenght)
    point1 = (line[0], line[1])
    point2 = (line[0] + vi[0], line[1] + vi[1])
    line = (*point1, *point2)
    return line

def makeArcLine(atom, alfa = 2*pi/3):
    lenght=BOUND_LENGHT
    line = atom.boundlist()[0].line
    point1 = atom.point
    if atom.point == line[2:4]:
        line = [*line[2:4], *line[0:2]]

    v = (line[2] - line[0], line[3] - line[1])
    vector = np.array(v)

    # define a rotation matrix to rotate the vector by 2/3 pi (120 degrees)
    rotation_matrix = np.array([[np.cos(alfa), -np.sin(alfa)],
                                [np.sin(alfa), np.cos(alfa)]])

    # multiply the rotation matrix by the input vector to get the new vector
    new_vector = np.dot(rotation_matrix, vector)

    vi = tuple(new_vector)
    point2 = (int(point1[0] + vi[0]), int(point1[1] + vi[1]))
    line = (*point1, *point2)
    return line

def makeCounterLine(atom, poly=False):
    lenght=BOUND_LENGHT
    
    point1 = atom.point
    v = (0, 0)
    for bound in atom.boundlist():
        line = bound.line
        if atom.point == line[2:4]:
            line = [*line[2:4], *line[0:2]]
        vec = (line[2] - line[0], line[3] - line[1])
        v = (v[0] + vec[0], v[1] + vec[1])

    cosa = v[0]/(v[0]**2 + v[1]**2)**(0.5)
    cosb = v[1]/(v[0]**2 + v[1]**2)**(0.5)
    vi = (-cosa*lenght, -cosb*lenght)

    if len(atom.boundlist()) % 2 != 0:
        vi = (cosa*lenght, cosb*lenght)
    if poly:
        vi = (-vi[0], -vi[1])
    point2 = (int(point1[0] + vi[0]), int(point1[1] + vi[1]))
    line = (*point1, *point2)
    return line

def lineFromAtom(atom):
    NomberOfBounds = len(atom.boundlist())
    if NomberOfBounds == 0:
        line = (*atom.point, atom.point[0] + 500, atom.point[1])
    elif NomberOfBounds == 1:
        line = makeArcLine(atom)
    else:
        line = makeCounterLine(atom)
    return line

def getCenterPolygon(line, r):
    x1 = line[0]
    y1 = line[1]
    x2 = line[2]
    y2 = line[3]
    v = (x2 - x1, y2 - y1)
    cosa = v[0]/(v[0]**2 + v[1]**2)**(0.5)
    cosb = v[1]/(v[0]**2 + v[1]**2)**(0.5)
    vi = (cosa*r, cosb*r)
    point2 = (line[2] + vi[0], line[3] + vi[1])
    return (*point2, cosa)

def makePointsOfPolygon(view, item, N):
    atoms =[]
    r = BOUND_LENGHT/(2*np.sin(pi/N))
    if not item:
        x = view.xStart
        y = view.yStart
        alfa = pi/2
    elif item.itemType == 'atom':
        atom = item
        if atom.boundlist and len(atom.boundlist()) == 1:
            line = atom.boundlist()[0].line
        elif not atom.boundlist():
            line = (*atom.point, atom.point[0] + BOUND_LENGHT, atom.point[1])
        else:
            line = makeCounterLine(atom, poly=True)
        if atom.point == line[0:2]:
            line = [*line[2:4], *line[0:2]]
        x , y , cosa = getCenterPolygon(line, r)
        alfa = pi - np.arccos(cosa)
        if line[3] > line[1]:
            alfa = - alfa
    elif item.itemType == 'bound':
        bound = item
        line = bound.line
        left, right = bound.findN()
        if left < right:
            line = (*line[2:4], *line[0:2])
        v = (line[2]- line[0], line[3]- line[1])
        wi = (v[0]**2 + v[1]**2)**0.5
        r = wi/(2*np.sin(pi/N))
        q, e, cosa = getCenterPolygon(line, r)
        al = np.arccos(cosa)
        if line[3] > line[1]:
            al = - al
        alfa = 3*pi/2 - al - pi/N

        re = r*np.sin(pi/2 - pi/N)
        v = (line[2]- (line[0]+line[2])/2, line[3]- (line[1]+line[3])/2)
        v2 = (-2*re*v[1]/wi, 2*re*v[0]/wi)
        x = (line[0]+line[2])/2 + v2[0]
        y = (line[1]+line[3])/2 + v2[1]
    
    for i in range(1, N+1):
        x1 = x + r * np.cos(2 * pi * i / N + alfa)
        y1 = y + r * np.sin(2 * pi * i / N + alfa)
        point = (int(x1), int(y1))
        atoms.append(point)
    return atoms

def makeShortLine(line):
    x1 = line[0].x()
    y1 = line[0].y()
    x2 = line[1].x()
    y2 = line[1].y()

    r = BOUND_LENGHT/10
    v = (x2 - x1, y2 - y1)
    if y2 < y1:
        v = (x1 - x2, y1 - y2)
    cosa = v[0]/(v[0]**2 + v[1]**2)**(0.5)
    sina = (1 - cosa**2)**0.5
    p1 = (x1 - r*cosa, y1 - r*sina)
    p2 = (x2 + r*cosa, y2 + r*sina)
    if y2>=y1:
        p1 = (x1 + r*cosa, y1 + r*sina)
        p2 = (x2 - r*cosa, y2 - r*sina)

    line = (*p1, *p2)
    line = QLine(*line)
    return line

