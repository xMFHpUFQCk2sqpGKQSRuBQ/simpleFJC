from distutils.core import setup
from Cython.Build import cythonize
from distutils.extension import Extension
from Cython.Distutils import build_ext

import numpy

ext_modules=[
    Extension("simpleFJC.integral"
              ,["simpleFJC/integral.pyx"]
              ,libraries=[""]
              #,extra_compile_args = ["-O3", "-ffast-math", "-march=native", "-fopenmp"]
              ,extra_link_args=['-fopenmp']
              ,include_dirs=[numpy.get_include()]
              ),
    Extension("simpleFJC.fjcBase"
              ,["simpleFJC/fjcBase.pyx"]
              ,libraries=[""]
              #,extra_compile_args = ["-O3", "-ffast-math", "-march=native", "-fopenmp"]
              ,extra_link_args=['-fopenmp']
              ,include_dirs=[numpy.get_include()]
              ),

]

setup(
    name='fjcBase',
    author='Gnintendo and Zzz',
    description='''
    The FJC functions for the integrals
    ''',
    version='2.0-b97',
    py_modules=['gaussxw'],
    license='Creative Commons Attribution-Noncommercial-Share Alike license',
    long_description=open('README.txt').read(),
    
    cmdclass = {"build_ext": build_ext},
    ext_modules = ext_modules
)
