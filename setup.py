from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
from os import system

# for notes on compiler flags e.g. using
# export CFLAGS=-O2
# http://docs.python.org/install/index.html

system('gfortran -m64 -O3 -c condegin13.f -o condegin13.o')

setup(
    cmdclass={'build_ext': build_ext},
    ext_modules=[Extension("eos", ["eos.pyx"], extra_objects=['condegin13.o'], 
    	libraries=['gfortran'],library_dirs=['/Applications/mesasdk/lib']),
    	Extension("crust", ["crust.pyx"])]
)
