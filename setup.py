# setup.py

from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy

# Define the extension module with NumPy include directory.
extensions = [
    Extension(
        name="Procedure",
        sources=["Procedure.pyx"],
        include_dirs=[numpy.get_include()],
        # You can define NPY_NO_DEPRECATED_API to silence warnings (optional)
        # define_macros=[('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION')],
    )
]

setup(
    name="Procedure",
    version="0.1",
    description="Cython module for k-mer pattern matching and sliding sums",
    ext_modules=cythonize(extensions, language_level="3"),
    python_requires=">=3.13",
)
