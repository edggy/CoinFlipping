# python setup.py build_ext --inplace
from distutils.core import setup, Extension

module1 = Extension('gf2c', sources = ['gf2/gf2.c'])
module2 = Extension('interpolateGF2', sources = ['gf2/gf2.c', 'interpolateGF2.c'], extra_compile_args=["-fPIC"])
module3 = Extension('ElGamalGF2', sources = ['gf2/gf2.c', 'ElGamal.c'], extra_compile_args=["-fPIC"])

setup (name = 'gf2c',
       version = '1.0',
        description = 'Galois field package',
        ext_modules = [module1])

setup (name = 'interpolateGF2',
       version = '1.0',
        description = 'GF2 interpolation package',
        ext_modules = [module2])

setup (name = 'ElGamalGF2',
       version = '1.0',
        description = 'ElGamal in GF2 package',
        ext_modules = [module3])
