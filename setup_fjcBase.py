from distutils.core import setup
from Cython.Build import cythonize
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy

ext_modules=[
    Extension("fjcBase",
              ["fjcBase.pyx"],
              include_dirs = [numpy.get_include()],
              extra_compile_args = [ ]
              ) 
]

setup(
    name='fjcBase',
    author='Gnintendo and Zzz',
    description='''
    The FJC functions for the integrals
    ''',
    version='2.0-b97',
    license='Creative Commons Attribution-Noncommercial-Share Alike license',
    long_description=open('README.txt').read(),
    
    cmdclass = {"build_ext": build_ext},
    ext_modules = ext_modules
)
