# drummer
![drum_demo_1.png](/images/drum_demo_1.png)

Drum simulation.  Call:

    $ ./gui.py

or alternatively:

    $ python3 gui.py

# Dependencies:
* Python 3
* PyQt5
* Pyaudio (>=0.2.7)
* Numpy (>=1.10.4)
* Scipy (>= 0.17.0)
* Matplotlib (>= 1.5.1)

You may find that older versions also work; these are just the version numbers I tested with.

# How it works:
The drum is modeled as a single, flat circular membrane with a fixed radius.  The half-decay time, surface tension, and mass density fully determine the time evolution of the vibrating waves (which follow a damped polar wave equation).  See [this wiki page](https://en.wikipedia.org/wiki/Vibrations_of_a_circular_membrane) for more info on the math (plus some nice visualizations).

An impact of the drum is modeled as a Dirac spike at the point hit, which is then decomposed into a linear superposition of Bessel functions (of the first kind).  The coefficients are determined by transforming the Dirac spike via a [Hankel transform](https://en.wikipedia.org/wiki/Hankel_transform).

The vibrating membrane is then converted into sound by modeling an observer some specified distance east of the drum.  At every point on the membrane, the height of the manifold is weighted multiplicatively by its distance squared to the observer (to model the inverse-squared decay of sound).
