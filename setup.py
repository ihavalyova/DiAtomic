#!/usr/bin/env python

import io
import os
from setuptools import setup, find_packages


name = 'diatomic-computations'
description = 'Computations of Diatomic Molecules.'
url = 'https://github.com/ihavalyova/DiAtomic'
email = 'havaliova@gmail.com'
author = 'Ilvie Havalyova'
requires_python = '>=3.6.0'
version = 'v0.0.1'
required_install = [
    'scipy>=1.5.0',
    'numpy>=1.16.0',
    'matplotlib>=3.2.0',
    'more_itertools',
    'ruamel.yaml',
    'numba'
]
required_extras = {
    'iminuit': ['iminuit'],
    'py3nj': ['py3nj']
}

here = os.path.abspath(os.path.dirname(__file__))
readme = os.path.join(here, 'README.md')

try:
    with io.open(readme, encoding='utf-8') as f:
        long_description = '\n' + f.read()
except FileNotFoundError:
    long_description = description


def run_setup():

    setup(
        name=name,
        version=version,
        description=description,
        long_description=long_description,
        long_description_content_type='text/markdown',
        author=author,
        author_email=email,
        python_requires=requires_python,
        url=url,
        download_url='',
        keywords=[
            'physics', 'diatomic-molecule', 'coupled-channels',
            'energy-levels', 'deperturbation',
        ],
        packages=find_packages(
            exclude=['test', '*.tests', 'test*']
        ),
        install_requires=required_install,
        extras_require=required_extras,
        package_data={
            'diatomic': []
        },
        include_package_data=True,
        license='BSD License',
        classifiers=[
            "License :: OSI Approved :: BSD License",
            "Programming Language :: Python :: 3.6",
            "Programming Language :: Python :: 3.7",
            "Programming Language :: Python :: 3.8",
            "Programming Language :: Python :: 3.9",
        ],
        zip_safe=False,
    )


if __name__ == '__main__':
    run_setup()
