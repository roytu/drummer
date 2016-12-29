
from PyQt5.QtWidgets import QWidget, QVBoxLayout
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

import numpy as np

from scipy.interpolate import griddata

class MplDrumWidget(QWidget):
    def __init__(self, parent=None):
        QWidget.__init__(self, parent)
        self.fig = Figure()
        #self.fig.subplots_adjust(left=0.2)
        self.ax1 = self.fig.add_subplot(111, projection="polar")
        self.fig.axes[0].get_yaxis().set_ticklabels([])

        self.canvas = FigureCanvas(self.fig)
        self.vbl = QVBoxLayout()
        self.vbl.addWidget(self.canvas)
        self.setLayout(self.vbl)

        self.heatmap = None

    def plot(self, drum, values):
        #self.ax1.cla()

        points = drum.points
        grid_r = drum.grid_r
        grid_theta = drum.grid_theta

        # Polar heatmap
        data = griddata(points, values, (grid_r, grid_theta),
                        method='cubic', fill_value=0)
        if not self.heatmap:
            self.heatmap = self.ax1.pcolormesh(drum.thetas, drum.rs, data.T)
                                               #shading="gouraud")
        else:
            C = data.T[:-1, :-1]  # TODO WTF
            # Possible bug in matplotlib:
            # http://stackoverflow.com/questions/18797175/animation-with-pcolormesh-routine-in-matplotlib-how-do-i-initialize-the-data/31490420#31490420
            # Whatever, this fixes it.  Turn this off for gouraud shading.

            self.heatmap.set_array(C.ravel())

        self.canvas.draw()
