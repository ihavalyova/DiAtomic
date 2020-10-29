#!/usr/bin/env python
# Author:  Ilvie Havalyova
# Contact: havaliova@gmail.com

import os
import pathlib
from setuptools import setup

setup_path = pathlib.Path(__file__).parent
readme_path = (os.path.join(setup_path, "README.md")).read_text()

def run_setup():

  setup(
    name="DiAtomic-Molecule-Calculations",
    version="1.0.0",
    description="Various Calculations of Diatomic Molecules",
    long_description=readme_path,
    long_description_content_type="text/markdown",
    url="https://github.com/ihavalyova/DiAtomic",
    author="Ilvie Havalyova",
    author_email="havaliova@gmail.com",
    license="BSD License",
    classifiers=[
        "License :: OSI Approved :: BSD License",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
    ],
    packages=["diatomic"],
    include_package_data=True,
    install_requires=["", ""],
#     entry_points={
#         "": [
#             "=.__main__:main",
#         ]
#     },
)

def main():
    run_setup()

if __name__ == '__main__':
    main()
