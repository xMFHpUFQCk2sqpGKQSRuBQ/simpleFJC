from distutils.core import setup
from Cython.Build import cythonize
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy

ext_modules=[
    Extension("gaussxw",
              ["gaussxw/gaussxw.pyx"],
              include_dirs = [numpy.get_include()],
              extra_compile_args = [ ]
              ) 
]

setup(
    name='gaussxw',
    author='Mark Newman',
    description='''
    Functions to calculate integration points and weights for Gaussian quadrature
    ''',
    version='1.0-b1',
    license='Creative Commons Attribution-Noncommercial-Share Alike license',
    long_description=open('README.txt').read(),

    cmdclass = {"build_ext": build_ext},
    ext_modules = ext_modules
)
