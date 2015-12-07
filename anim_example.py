# -*- coding: utf-8 -*-
from __future__ import print_function
import datetime
import subprocess
import numpy as np
from matplotlib import pyplot as plt
from skyplot_snippet import plot_sky

fig = plt.figure(figsize=(14, 8))

# a single subplot filling the whole figure w/ mollweide projection
sky_ax = fig.add_subplot(1, 1, 1, projection='mollweide')

# make a fake catalog (coords in degrees)
npoints = 100
el, eb = 360 * np.random.rand(npoints), 180 * (np.random.rand(npoints) - 0.5)

start_date = datetime.date(2018, 10, 1)
# number of days (frames) in animation
ndays = 100

for i in range(ndays):
    sky_ax.clear()
    plot_sky(sky_ax, el, eb, for_date=start_date + datetime.timedelta(days=i))
    sky_ax.set_title('Ecliptic Sky at {} + {}'.format(start_date, i))
    fn = './frame_{:03}.png'.format(i)
    plt.savefig(fn)
    print("Saved {}".format(fn))

# Use a tool like avconv or ffmpeg to convert the resulting PNG frames:
#
# avconv -framerate 25 -f image2 -i frame_%03d.png -c:v h264 -crf 1 out.mov