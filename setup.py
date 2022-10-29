#!/usr/bin/env python

# Procedure for local instalation:
# conda create --name diatom  # create new env
# conda activate diatom (or source activate diatom) # actvate env
# conda install numpy, scipy, matplotlib
# pip install pybind11
# python setup.py install

import os

try:
    from setuptools import setup, Extension
except ImportError:
    from distutils.core import setup
    from distutils.extension import Extension

here = os.path.abspath(os.path.dirname(__file__))
readme = os.path.join(here, 'README.md')


class get_pybind11_include_dirs(object):
    """Returns a list of include directories for pybind11;
    Ensure that pybind11 is installed before referencing it."""

    def __str__(self):
        import pybind11
        return pybind11.get_include()


os.environ['CC'] = 'g++'

diatomic_ext1 = Extension(
    "diatomic.identify",
    sources=["diatomic/identify.cpp"],
    extra_compile_args=['-O3', '-shared', '-std=c++11'],
    include_dirs=[get_pybind11_include_dirs()]
)

diatomic_ext2 = Extension(
    "diatomic.wavenumbers",
    sources=["diatomic/wavenumbers.cpp"],
    extra_compile_args=['-O3', '-shared', '-std=c++11'],
    include_dirs=[get_pybind11_include_dirs()]
)

setup(
    name='diatomic-computations',
    version='v0.0.1',
    description='Diatomic Molecules Computations',
    long_description=open(readme).read(),
    long_description_content_type='text/markdown',
    author='Ilvie Havalyova',
    author_email='havaliova@gmail.com',
    python_requires='>=3.7.0',
    url='https://github.com/ihavalyova/DiAtomic',
    download_url='',
    license='BSD3',
    keywords=[
        'physics',
        'diatomic-molecule',
        'coupled-channels',
        'energy-levels',
        'deperturbation',
    ],
    packages=['diatomic'],
    install_requires=[
        'numpy>=1.16.0',
        'scipy>=1.5.0',
        'matplotlib>=3.2.0',
    ],
    extras_require={
        'iminuit': ['iminuit'],
        'py3nj': ['py3nj']
    },
    package_data={
        'diatomic': ['identify.cpp', 'wavenumbers.cpp']
    },
    include_package_data=True,
    classifiers=[
        "License :: OSI Approved :: BSD License",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
    ],
    zip_safe=False,
    ext_modules=[diatomic_ext1, diatomic_ext2],
)
