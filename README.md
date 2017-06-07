### crustpy

This code follows the thermal evolution of a neutron star crust. To compile the eos and crust modules with cython (`pip install cython`), use

`python setup.py build_ext --inplace`

An [IPython notebook](https://github.com/andrewcumming/crustpy/blob/master/crustpy.ipynb) is included that shows how to run a model of the MXB 1659 lightcurve (see [Brown & Cumming 2009](http://iopscience.iop.org/0004-637X/698/2/1020)). Alternatively, use `python crustcool.py` which plots the results in `crustcool.pdf`.

![MXB 1659 lightcurve](https://github.com/andrewcumming/crustpy/raw/master/lc1659.png)

During accretion, the temperature at the top of the grid can be held constant, or allowed to cool and a shallow heat source included at low density.

Interactive version using [bokeh](http://bokeh.pydata.org/en/latest/):
`bokeh serve --show crustcool_bokeh.py`

To-do list:
* add entrainment of neutrons by the nuclei in the inner crust (no entrainment assumed at the moment)
