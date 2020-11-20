![](./doc/logo.png)

---

[![build](https://github.com/ihavalyova/DiAtomic/workflows/build/badge.svg?branch=master)](https://github.com/ihavalyova/DiAtomic/actions)
[![Language grade: Python](https://img.shields.io/lgtm/grade/python/g/ihavalyova/DiAtomic.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/ihavalyova/DiAtomic/context:python)
[![Total alerts](https://img.shields.io/lgtm/alerts/g/ihavalyova/DiAtomic.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/ihavalyova/DiAtomic/alerts/)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/ihavalyova/DiAtomic/master)

<!-- https://mybinder.org/v2/gh/ihavalyova/DiAtomic/master -->

DiAtomic is an open-source library written in Python with useful routines for performing various computations of diatomic molecules. It can be used to:

- compute the energy eigenvalues and eigenvectors for
  - single non-interacting electronic state
  - system of arbitrary number of coupled electronic states
  
  by solving the radial Schrodinger equation
- compute the transition frequencies and line intensities when PECs or term values are available
- fit of the computed energy levels to the experimental data

It supports object-oriented approach with easy to use functionality and efficiently vectorized code.

**To install**: call Python pip from the command line:

```pip install diatomic-computations```

**Documentation**: <a href="https://github.com/ihavalyova/DiAtomic/blob/master/doc/Diatomic.md" target="_blank">Detailed version</a>

**Example notebooks:** <a href="https://github.com/ihavalyova/DiAtomic/blob/master/doc/" target="_blank">Jupyter notebooks</a>

**Source:** <a href="https://github.com/ihavalyova/DiAtomic/tree/master/diatomic" target="_blank">Python source</a>

**Cite as:**

I. Havalyova and A. Pashov: Аn open-source library for calculating energy levels and spectra of diatomic molecules (2020) https://doi.org/10.5281/zenodo.4028548
