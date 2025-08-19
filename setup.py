from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy

extensions = [
    Extension(
        name="Procedure",
        sources=["Procedure.pyx"],
        include_dirs=[numpy.get_include()],
        # If you hit MSVC warnings, you can add extra_compile_args here
        # extra_compile_args=["/O2"],  # Windows (MSVC)
        # extra_compile_args=["-O3"],  # GCC/Clang
    )
]

setup(
    name="Procedure",
    version="0.2.0",
    description="Cython module for k-mer pattern matching and sliding sums",
    ext_modules=cythonize(extensions, language_level="3"),
)
