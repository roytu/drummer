
from plot_widget import CustomPlotWidget

class DrumPlotWidget(CustomPlotWidget):
    def __init__(self, app):
        CustomPlotWidget.__init__(self)

    def update(self, points, values):


    def on_plot(self):
        self.plot_func()
