from PyQt5.QtGui import QPen, QColor
from PyQt5.QtCore import Qt


BOUND_LENGHT = 750
BOUND_WIDTH = 35
ATOM_RADIUS = 100
ATOM_SYMBOL_SIZE =200

Redpen = QPen(Qt.red)
Redpen.setWidth(15)

Bluepen = QPen(Qt.blue)
Bluepen.setWidth(15)

Whitepen = QPen(QColor(255,255,255))

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
