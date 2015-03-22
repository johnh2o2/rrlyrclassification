from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy
ext_modules = [Extension('fastperiod',['fastperiod.pyx'], include_dirs=[numpy.get_include()]), 
				Extension('fastlombscargle',['lsp_faster.pyx'], include_dirs=[numpy.get_include()])]


setup(
    name = 'Cython versions of base period finding algorithms',
    cmdclass = {'build_ext':build_ext},
    ext_modules = ext_modules
    )

