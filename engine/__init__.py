from .chemdr import *
import base64

from PyQt5.QtGui import QIcon, QPixmap

from .icons import icons


def get_icon(name):
    pixmap = QPixmap()
    pixmap.loadFromData(base64.b64decode(icons[name]))
    return QIcon(pixmap)
