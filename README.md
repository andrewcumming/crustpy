### crustpy

This code follows the thermal evolution of a neutron star crust.
An [IPython notebook](https://github.com/andrewcumming/crustpy/blob/master/crustpy.ipynb) is included that shows how to run a model of the MXB 1659 lightcurve (see [Brown & Cumming 2009](http://iopscience.iop.org/0004-637X/698/2/1020)). Alternatively, use `python crustcool.py` which plots the results in `crustcool.pdf`.

![MXB 1659 lightcurve](https://github.com/andrewcumming/crustpy/raw/master/lc1659.png)

During accretion, the temperature at the top of the grid can be held constant, or allowed to cool and a shallow heat source included at low density.

Note the limited microphysics included in this code (to-do list):
* it is assumed the neutrons do not contribute to heat capacity (density-dependent gap not included)
* Coulomb log for electron-phonon scattering is set to unity
* no entrainment of neutrons by the nuclei in the inner crust
* neutron pressure taken from Mackie & Baym

