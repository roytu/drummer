#!/usr/bin/python3

import sys
import pyaudio as pya

import numpy as np

from PyQt5.Qt import QApplication
from PyQt5.QtGui import QIntValidator
from PyQt5.QtWidgets import QMainWindow
from PyQt5.QtCore import QTimer

from scipy.io import wavfile

from gui_ import Ui_Drummer
from drum import Drum

class Gui(QApplication):
    TIMESTEP = float(1) / 44100  # seconds
    TIMELENGTH = 1.0  # seconds

    def __init__(self, args):
        QApplication.__init__(self, args)

        # Initialize window
        self.mainWindow = QMainWindow()
        self.ui = Ui_Drummer()
        self.ui.setupUi(self.mainWindow)
        self.mainWindow.show()

        # Initialize timers
        self.tickTimer = QTimer()
        self.tickTimer.setInterval(10)  # 10 ms per tick

        # Connect signals
        # Drum sliders
        self.ui.radiusSlider.valueChanged.connect(self._on_drum_change)
        self.ui.decaySlider.valueChanged.connect(self._on_drum_change)
        self.ui.tensionSlider.valueChanged.connect(self._on_drum_change)
        self.ui.densitySlider.valueChanged.connect(self._on_drum_change)

        # Simulation sliders
        self.ui.rStepsSlider.valueChanged.connect(self._on_drum_change)
        self.ui.thetaStepsSlider.valueChanged.connect(self._on_drum_change)
        self.ui.maxMSlider.valueChanged.connect(self._on_drum_change)
        self.ui.maxNSlider.valueChanged.connect(self._on_drum_change)

        self.ui.timeSlider.valueChanged.connect(self._on_time_change)
        self.ui.runningCheckBox.stateChanged.connect(self._on_running_change)
        self.tickTimer.timeout.connect(self._on_tick)
        self.hit_cid = self.ui.drumWidget.canvas.mpl_connect("button_press_event", self._on_hit)
        self.aboutToQuit.connect(self._on_quit)

        # Slider -> LineEdit signals
        def link(slider, lineedit):
            """ Corresponds a Slider with a LineEdit.  An update to either
                updates both.  Also updates lineedit with the current
                slider value.
            """
            def update_lineedit():
                lineedit.setText(str(slider.value()))
            def update_slider():
                slider.setValue(int(lineedit.text()))
            min_ = slider.minimum()
            max_ = slider.maximum()
            validator = QIntValidator(min_, max_)
            lineedit.setValidator(validator)

            slider.sliderMoved.connect(update_lineedit)
            lineedit.textEdited.connect(update_slider)

            # Update line edit
            lineedit.setText(str(slider.value()))

        link(self.ui.radiusSlider, self.ui.radiusLineEdit)
        link(self.ui.decaySlider, self.ui.decayLineEdit)
        link(self.ui.tensionSlider, self.ui.tensionLineEdit)
        link(self.ui.densitySlider, self.ui.densityLineEdit)
        link(self.ui.rStepsSlider, self.ui.rStepsLineEdit)
        link(self.ui.thetaStepsSlider, self.ui.thetaStepsLineEdit)
        link(self.ui.maxMSlider, self.ui.maxMLineEdit)
        link(self.ui.maxNSlider, self.ui.maxNLineEdit)
        link(self.ui.timeSlider, self.ui.timeLineEdit)
        link(self.ui.speedSlider, self.ui.speedLineEdit)
        link(self.ui.distanceSlider, self.ui.distanceLineEdit)

        # Initialize drum
        self.drum = None
        self._on_drum_change()
        self.drum.hit((0.8, 0, 0), 1)


        # Initialize audio stufff
        self.pyaudio = pya.PyAudio()
        self.stream = self.pyaudio.open(format=pya.paFloat32,
                                        channels=1,
                                        rate=44100,
                                        output=True)


        # Empty arrays
        #self.values = np.zeros((int(Gui.TIMELENGTH / Gui.TIMESTEP),
        #                        self.ui.rStepsSlider.value() * self.ui.thetaStepsSlider.value()))
        self.values = np.array([])

        # Start timer
        self.tickTimer.start()

    def _on_drum_change(self):
        a = float(self.ui.radiusSlider.value()) / 100
        K = 1 / (float(self.ui.decaySlider.value()) / 1000)
        tension = float(self.ui.tensionSlider.value())
        density = float(self.ui.densitySlider.value()) / 1000
        c = (tension / density) ** 0.5

        rcount = float(self.ui.rStepsSlider.value())
        thetacount = float(self.ui.thetaStepsSlider.value())
        m_max = int(self.ui.maxMSlider.value())
        n_max = int(self.ui.maxNSlider.value())

        print("Settings")
        print("=" * 20)
        print("a: {0}".format(a))
        print("K: {0}".format(K))
        print("c: {0}".format(c))
        print("rcount: {0}".format(rcount))
        print("thetacount: {0}".format(thetacount))
        print("m_max: {0}".format(m_max))
        print("n_max: {0}".format(n_max))

        self.drum = Drum(a=a, K=K, c=c, rcount=rcount, \
                         thetacount=thetacount, m_max=m_max, n_max=n_max)

    def _on_time_change(self):
        # Update plot
        if self.values.size != 0:
            t = float(self.ui.timeSlider.value()) / 1000
            i = int(t / Gui.TIMESTEP)
            self.ui.drumWidget.plot(self.drum, self.values[i])

    def _on_tick(self):
        if self.ui.runningCheckBox.isChecked():
            max_ = self.ui.timeSlider.maximum()
            t = float(self.ui.timeSlider.value())
            t = (t + 10 * (float(self.ui.speedSlider.value()) / 100)) % max_
            self.ui.timeSlider.setValue(t)

    def _on_running_change(self):
        if self.ui.runningCheckBox.isChecked():
            self.tickTimer.start()
        else:
            self.tickTimer.stop()

    def _on_hit(self, event):
        r = event.ydata
        theta = event.xdata

        # TODO variable force
        self.drum.hit((r, theta, 0), 1)

        # Reset time
        self.ui.timeSlider.setValue(0)

        # Calculate drum values
        times = np.arange(0, Gui.TIMELENGTH, Gui.TIMESTEP)
        self.values = self.drum.value(times)

        # Wave values and plot
        d = self.ui.distanceSlider.value()
        wave_values = self.drum.wave_value_from_values(self.values, d)
        self.ui.plotWidget.plot(times, wave_values)

        # Play sound if on
        if self.ui.soundGroupBox.isChecked():
            # Calculate the hit for next second
            SAMPLE_RATE = 44100  # Hz
            SAMPLE_TIME = 1  # seconds
            CHUNK = 1024

            dts = np.linspace(0, SAMPLE_TIME, SAMPLE_RATE * SAMPLE_TIME)
            samples = wave_values

            # Normalize
            samples = np.array(samples / np.max(np.abs(samples)), dtype=np.float32)
            fname = self.ui.wavfnameLineEdit.text()
            if fname != "":
                wavfile.write(fname, SAMPLE_RATE, samples)

            data = samples.tostring()
            self.stream.write(data)

    def _on_quit(self):
        self.ui.drumWidget.canvas.mpl_disconnect(self.hit_cid)

        self.stream.close()    
        self.pyaudio.terminate()

if __name__ == "__main__":
    gui = Gui(sys.argv)
    gui.exec_()
