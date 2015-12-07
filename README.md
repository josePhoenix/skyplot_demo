# skyplot_snippet

Example of how to use matplotlib and some spherical geometry to plot a field of regard on the sky in ecliptic coordinates. Useful when you have a spaceborne observatory that has a certain sun-avoidance angle and a catalog whose availability you want to visualize.

## Installation

This doesn't really need installation, since it's just a snippet. If you want to use it as-is, drop `skyplot_snippet.py` into the same folder as your script and use `from skyplot_snippet import ...` as in the examples provided.

The example depends on NumPy, Matplotlib, and PyEphem. It's been tested with NumPy 1.9.3, Matplotlib 1.5.0, and PyEphem 3.7.6.0. Ensure you have at least those versions using this `pip install` command:

    pip install "matplotlib>=1.5.0" "numpy>=1.9.3" "pyephem>=3.7.6.0"
