
from PyQt5.QtWidgets import QWidget, QVBoxLayout
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

import numpy as np

from scipy.interpolate import griddata

class MplPlotWidget(QWidget):
    def __init__(self, parent=None):
        QWidget.__init__(self, parent)
        self.fig = Figure()
        #self.fig.subplots_adjust(left=0.2)
        self.ax1 = self.fig.add_subplot(111)
        self.canvas = FigureCanvas(self.fig)
        self.vbl = QVBoxLayout()
        self.vbl.addWidget(self.canvas)
        self.setLayout(self.vbl)

        self.line, = self.ax1.plot([], [])

        self.ax1.set_xlim(0, 1)
        self.ax1.set_ylim(-1, 1)

    def plot(self, ts, values):
        #self.ax1.cla()

        # Normalize
        values /= np.max(np.abs(values))

        # Flip if weird
        if values[0] < 0:
            values *= -1

        self.line.set_xdata(ts)
        self.line.set_ydata(values)

        self.canvas.draw()
