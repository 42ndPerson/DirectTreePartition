from setuptools import setup, Extension
from Cython.Build import cythonize

extensions = [
    Extension("CythonCore.Euclid", ["CythonCore/Euclid.pyx"]),
    Extension("CythonCore.TwinGraph", ["CythonCore/TwinGraph.pyx"]),
    Extension("CythonCore.RegionTree", ["CythonCore/RegionTree.pyx"]),
    Extension("CythonCore.GraphNav", ["CythonCore/GraphNav.pyx"]),
    Extension("CythonCore.TwinGraphWilson", ["CythonCore/TwinGraphWilson.pyx"]),
]

setup(
    name="CythonCore",
    ext_modules=cythonize(extensions, compiler_directives={'language_level': "3"}),
)
