from distutils.core import setup
from Cython.Build import cythonize
from distutils.extension import Extension
from Cython.Distutils import build_ext

setup(
    name='gaussxw',
    author='Mark Newman',
    description='''
    A simple tool to perform Freely-Joint Chain calculations
    using numerical integration using Monte Carlo techniques.
    Based on the original mcint modul of Tristan Snowsill.
    Multi-process CYthon version.
    ''',
    version='1.0-b1',
    py_modules=['gaussxw'],
    license='Creative Commons Attribution-Noncommercial-Share Alike license',
    long_description=open('README.txt').read(),
    
    cmdclass = {"build_ext": build_ext},
    ext_modules = ext_modules
)
