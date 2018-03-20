# python setup.py build_ext --inplace
from distutils.core import setup, Extension

interpolateModule = Extension('interpolateGF2', sources = ['gf2/gf2.c', 'interpolateGF2.c'], extra_compile_args=["-fPIC"])
ElGamalModule = Extension('ElGamalGF2', sources = ['gf2/gf2.c', 'ElGamal.c'], extra_compile_args=["-fPIC"])

setup (name = 'coinFlipping',
       version = '1.0',
       description = 'Multiparty coin flipping',
       py_modules = ['coinFlipping'])

setup (name = 'interpolateGF2',
       version = '1.0',
       description = 'GF2 interpolation package',
       ext_modules = [interpolateModule])

setup (name = 'ElGamalGF2',
       version = '1.0',
       description = 'ElGamal in GF2 package',
       ext_modules = [ElGamalModule])
