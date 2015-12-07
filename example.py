# -*- coding: utf-8 -*-
from __future__ import print_function
import numpy as np
from matplotlib import pyplot as plt
from skyplot_snippet import plot_sky

fig = plt.figure(figsize=(14, 8))

# a single subplot filling the whole figure w/ mollweide projection
sky_ax = fig.add_subplot(1, 1, 1, projection='mollweide')

# make a fake catalog (coords in degrees)
npoints = 10000
el, eb = 360 * np.random.rand(npoints), 180 * (np.random.rand(npoints) - 0.5)

plot_sky(sky_ax, el, eb)
plt.show()
