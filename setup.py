from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

# for notes on compiler flags e.g. using
# export CFLAGS=-O2
# http://docs.python.org/install/index.html

setup(
    cmdclass={'build_ext': build_ext},
    ext_modules=[Extension("eos", ["eos.pyx"]),Extension("crust", ["crust.pyx"])]
)
