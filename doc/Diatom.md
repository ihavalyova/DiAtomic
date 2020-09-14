
- [**DiAtom** module: what is it used for?](#diatom-module-what-is-it-used-for)
- [**DiAtom**  module: how to install and setup?](#diatom-module-how-to-install-and-setup)
- [Diatomic molecule: basic theoretical concepts](#diatomic-molecule-basic-theoretical-concepts)
  - [The total Hamiltonian and basis functions](#the-total-hamiltonian-and-basis-functions)
  - [The Scrodinger equation for a single state and coupled system of states](#the-scrodinger-equation-for-a-single-state-and-coupled-system-of-states)
  - [The interaction terms and their matrix elements](#the-interaction-terms-and-their-matrix-elements)
  - [Methods for solution of the Schrodinger equation](#methods-for-solution-of-the-schrodinger-equation)
  - [Potential Energy function models (PECs models)](#potential-energy-function-models-pecs-models)
- [Computing Energy Eigenvalues](#computing-energy-eigenvalues)
  - [Molecule Data Object Definition](#molecule-data-object-definition)
  - [Grid Object Definition](#grid-object-definition)
  - [Channel Object Definition](#channel-object-definition)
  - [Coupling Object Definition](#coupling-object-definition)
  - [Molecule Levels Computation](#molecule-levels-computation)
  - [Examples](#examples)
- [Fitting of the Calculated Energy Levels](#fitting-of-the-calculated-energy-levels)
  - [SVD Fit](#svd-fit)
  - [Minuit Fit](#minuit-fit)
  - [Levenberg-Marquard Fit](#levenberg-marquard-fit)
- [Computing the Transition Frequencies and Intensities](#computing-the-transition-frequencies-and-intensities)
  - [States represented by channels](#states-represented-by-channels)
  - [States represented by term values](#states-represented-by-term-values)
  - [Line strength](#line-strength)
    - [Honl-London Factors](#honl-london-factors)
    - [Frank-Condon Factors](#frank-condon-factors)
  - [Einstein A coefficient and Radiative Lifetime](#einstein-a-coefficient-and-radiative-lifetime)
- [Plotting](#plotting)


# **DiAtom** module: what is it used for?
The Python package **```DiAtom```** allows various calculations for diatomic molecules to be performed. It supports single and coupled channels computations of bound rovibrational levels, intensity calculations, fitting to the experimental data.
The current functionality covered by the program includes:
* ..
* ..

Just as an example of what you can do with **```DiAtom```** module, if you have a set of potentials for a couple of molecular electronic states represented by points (abinitio, RKR and etc.) and want to see how they look you can do something like:
<!-- 
```python
p = diatom.Plotting()
p.plot_potentials_points(['p1.pot', 'p2.pot', 'p3.pot', 'p4.pot'], show=True, ipoints=50, ylim=(1e3, 1e5))
```
or even simpler:
-->

```python
import glob

p = diatom.Plotting()
p.plot_potentials_points(glob.glob('./*.pot'), show=True, ipoints=120, xlim=(2.5, 16), ylim=(9e3, 2.5e4))
```
assuming your potential files are in the current directory.


![KCs_potentials](./plotting/kcs_potential_points.svg)

# **DiAtom**  module: how to install and setup?

**```DiAtom```**  module can be installed from the Python software repository PyPI (Python Package Index) via pip. From Linux command line execute

```console
pip install diatom
```

and from Jupyter or IPython execute

````python
In [1]: ! pip install diatom
````
To quickly check whether the installation has been successful type

```console
python
>>> import diatom
>>> diatom
```

and the path to the \_\_init\_\_.py file in the install location should be outputed.

After installing create a new python file for example called main.py and import the diatom module

```python
#!/usr/bin/env python

from diatom import *
```
To execute the file from the Linux command line write
```console
python main.py
```

or type

```console
chmod u+x main.py
./main.py
```

to make the file executable and run it. To execute the file from the interactive shell of IPython (Interactive Python) type ipython then

```python
In [1]: run main.py
```

The **```DiAtom```** module is extensivly tested on Linux platform but works under Windows and MacOS as well.

# Diatomic molecule: basic theoretical concepts

## The total Hamiltonian and basis functions

The total Hamiltonian of a diatomic molecule in the rotating molecule-fixed coordinate system with origin at the center of mass of the molecule can be written as a sum of several terms:

<p align="center"><img src="/doc/tex/0befcf4d6f64295f51648ed6605604bf.svg?invert_in_darkmode&sanitize=true" align=middle width=381.66411525pt height=16.438356pt/></p>

where <img src="/doc/tex/88eff123eb81302eced70d075cfd831d.svg?invert_in_darkmode&sanitize=true" align=middle width=49.06969154999999pt height=24.65753399999998pt/> and <img src="/doc/tex/1565b8417026bdbae4f63e0a7df02e21.svg?invert_in_darkmode&sanitize=true" align=middle width=90.41667525pt height=24.65753399999998pt/> are the vibrational and rotational part of the total nuclear kinetic energy operator in spherical polar coordinates, <img src="/doc/tex/c3a1d7f8e57bd6a18adbaa8585b1b61b.svg?invert_in_darkmode&sanitize=true" align=middle width=40.47569294999999pt height=24.65753399999998pt/> is the kinetic energy of the electrons, 

<!-- <p align="center"><img src="/doc/tex/33361e7842bc9784188a404d33d61dfb.svg?invert_in_darkmode&sanitize=true" align=middle width=107.76820724999999pt height=16.438356pt/></p> -->

## The Scrodinger equation for a single state and coupled system of states

<!-- omit in toc -->
### Single channel approximation for an isolated state

The energy eigenvalues of single isolated state of a diatomic molecule can be obtained by solving the radial Schrodinger equation

<p align="center"><img src="/doc/tex/04ec33ae2ec81c77f6d478424e026fe4.svg?invert_in_darkmode&sanitize=true" align=middle width=428.47916925pt height=40.11819404999999pt/></p>

with internuclear distance labeled with <img src="/doc/tex/1e438235ef9ec72fc51ac5025516017c.svg?invert_in_darkmode&sanitize=true" align=middle width=12.60847334999999pt height=22.465723500000017pt/>, the sum of the second and the third term <img src="/doc/tex/f11fee23e73e65ca5d56841195727865.svg?invert_in_darkmode&sanitize=true" align=middle width=196.35702239999998pt height=26.76175259999998pt/> from the left hand side is called effective potential energy curve, the reduced mass is <img src="/doc/tex/9c774da583d642ab71cf6fc72e1cb939.svg?invert_in_darkmode&sanitize=true" align=middle width=166.2062853pt height=24.65753399999998pt/> where <img src="/doc/tex/6f549764f2f97bec950c14de5352994a.svg?invert_in_darkmode&sanitize=true" align=middle width=22.500061649999992pt height=22.465723500000017pt/> and <img src="/doc/tex/dced8cd0d35e2af2d3499c10d7ee6289.svg?invert_in_darkmode&sanitize=true" align=middle width=22.500061649999992pt height=22.465723500000017pt/> are the masses of the two atoms, J is the rotational quantum number; <img src="/doc/tex/311f29e88a451816c3e621e910343bf4.svg?invert_in_darkmode&sanitize=true" align=middle width=27.481438049999987pt height=22.465723500000017pt/> are the energy eigenvalues of the rovibrational levels and <img src="/doc/tex/6130302efada133fa7ea1cd1ccb12a8f.svg?invert_in_darkmode&sanitize=true" align=middle width=25.14126119999999pt height=22.831056599999986pt/> are their corresponding eigenfunctions.

<!-- omit in toc -->
### The coupled channels problem

## The interaction terms and their matrix elements

The most important operators and their matrix elements are:

- Spin-Orbit

- L-Uncoupling

- Spin-Uncoupling

- Spin-Electornic

- Spin-Rotation

- Spin-Spin

- <img src="/doc/tex/b23332f99af850a48831f80dbf681ed6.svg?invert_in_darkmode&sanitize=true" align=middle width=11.41554479999999pt height=22.465723500000017pt/> and <img src="/doc/tex/9432d83304c1eb0dcb05f092d30a767f.svg?invert_in_darkmode&sanitize=true" align=middle width=11.87217899999999pt height=22.465723500000017pt/> Doubling

## Methods for solution of the Schrodinger equation

Finite-Difference and Fourier Grid Hamiltonain (DVR type method) are the most frequently applied methods for numerical solution of the 1D Schordinger equation for single and coupled channels problems in molecular spectroscopy. In both methods the wavefunction is approximated over an equidistant or non-equidistant grid of points.

<!-- omit in toc -->
#### Uniform grid

In this case the grid points <img src="/doc/tex/6660896e4379722ff79bba94961b201c.svg?invert_in_darkmode&sanitize=true" align=middle width=17.132374049999992pt height=22.465723500000017pt/> in the interval from <img src="/doc/tex/3076bed5ea0f6f7be4dcc194b5375a3b.svg?invert_in_darkmode&sanitize=true" align=middle width=36.92324789999999pt height=22.465723500000017pt/> to <img src="/doc/tex/ea449f9e9a48e2959872aac8fa65e1ca.svg?invert_in_darkmode&sanitize=true" align=middle width=38.73108029999999pt height=22.465723500000017pt/> are determined by:

<p align="center"><img src="/doc/tex/6f277b0df057cb1fe10512647772c30a.svg?invert_in_darkmode&sanitize=true" align=middle width=292.75927724999997pt height=16.438356pt/></p>

where <img src="/doc/tex/83b80e12323c9ced450508d7bd3822a5.svg?invert_in_darkmode&sanitize=true" align=middle width=23.169583799999987pt height=22.465723500000017pt/> is the number of grid points and <img src="/doc/tex/5385180385f31e03adfa5256148b06bb.svg?invert_in_darkmode&sanitize=true" align=middle width=23.66048849999999pt height=22.465723500000017pt/> is the grid step

<p align="center"><img src="/doc/tex/f17deb669038ee0be520dfb23a4a35c0.svg?invert_in_darkmode&sanitize=true" align=middle width=145.76189265pt height=36.09514755pt/></p>

<!-- omit in toc -->
#### Nonuniform grid



<!-- omit in toc -->
  ### Finite-Difference Method (FD)

The second derivative of the wavefunction with respect to the internuclear distance is approximated by five-point central difference schema:

<p align="center"><img src="/doc/tex/ad06c975f5915bd636d67109e3ca03c9.svg?invert_in_darkmode&sanitize=true" align=middle width=465.42534059999997pt height=39.452455349999994pt/></p>

The kinetic energy matrix elements are then computed:

<p align="center"><img src="/doc/tex/f27a85c449f4aa2a88eca7afdbe3d420.svg?invert_in_darkmode&sanitize=true" align=middle width=368.15905499999997pt height=38.973783749999996pt/></p>

This is a banded symmetric matrix. The potential energy matrix is diagonal:

<p align="center"><img src="/doc/tex/8f9f8303492ff77504d46be1964f9256.svg?invert_in_darkmode&sanitize=true" align=middle width=105.12704895pt height=17.031940199999998pt/></p>


<!-- omit in toc -->
  ### Fourier Grid Hamiltonian (FGH)

FGH is a type of a collocation (pseudospectral) method for which the solution is approximated over a special grid points called collocation points.

<!-- omit in toc -->
#### **Sinc basis**

In sinc basis the kinetic energy matrix elements are computed as:

<!-- omit in toc -->
#### **Fourier basis**

## Potential Energy function models (PECs models)
- Pointwise potential

- Analytical potential
   - Morse potential

  The Morse potential function is defined as:

  <p align="center"><img src="/doc/tex/32aad865ba81d6b7c52a0cb268f1a101.svg?invert_in_darkmode&sanitize=true" align=middle width=217.99017239999998pt height=19.526994300000002pt/></p>
  where <img src="/doc/tex/d374a50b5b4324e84b1ade9e26ae25e3.svg?invert_in_darkmode&sanitize=true" align=middle width=15.84310529999999pt height=22.465723500000017pt/> is the term value, <img src="/doc/tex/1876b29402de38408e13f33010df5c24.svg?invert_in_darkmode&sanitize=true" align=middle width=19.84651514999999pt height=22.465723500000017pt/> measures the energy from the bottom of the potential to the dissociation limit (it is not the dissociation energy)
  <img src="/doc/tex/8217ed3c32a785f0b5aad4055f432ad8.svg?invert_in_darkmode&sanitize=true" align=middle width=10.16555099999999pt height=22.831056599999986pt/> is a constant and <img src="/doc/tex/7b6a7362f072c7598f6fdb54363c97e8.svg?invert_in_darkmode&sanitize=true" align=middle width=13.65323849999999pt height=14.15524440000002pt/> is the equilibrium internuclear distance.

   - EMO (Expanded Morse Oscillator) potential

  It is defined as

  <p align="center"><img src="/doc/tex/e88b8315d59ce4cb6c5fa4897772e628.svg?invert_in_darkmode&sanitize=true" align=middle width=383.9485419pt height=18.312383099999998pt/></p>

  which is a Morse potential function with a radial dependence of the exponential coefficient

  <p align="center"><img src="/doc/tex/f1e0f767c174f6661b8264ae68817fec.svg?invert_in_darkmode&sanitize=true" align=middle width=176.35144889999998pt height=47.806078649999996pt/></p>

  represented as a polynomial expansion in the powers of the dimensionless function

  <p align="center"><img src="/doc/tex/ba1da685931183a2ce899c4a8fbe37c3.svg?invert_in_darkmode&sanitize=true" align=middle width=118.4382078pt height=35.551675499999995pt/></p>

  with p being a positive integer number.

   - MLR (Morse/Long-Range) potential

<p align="center"><img src="/doc/tex/826ef0830ed7bf386d12de409f2bb662.svg?invert_in_darkmode&sanitize=true" align=middle width=106.3151496pt height=16.438356pt/></p>

# Computing Energy Eigenvalues

## Molecule Data Object Definition

Initially an object of type **```MoleculeData```** should be instantiated like:

```python
# creating MoleculeData object called mdata
mdata = MoleculeData()
```

The basic input information about the molecule should be defined by the following properties via the created **```MoleculeData```** object:

<!-- omit in toc -->
### Reduced masses and isotopes

There are two ways for specifing the reduced mass.

- **```molecule```** - defines one or more isotopic forms of the same molecule by specifing their chemical symbols i.e. each defined item corresponds to a diffrent isotope of the same molecule.
  - the molecule symbols should be in the format: 'M1A1M2A2' where M1/M2 are the mass numbers and A1/A2 are the chemical symbols of the first/second atom (spaces are allowed).
  - The reduced masses for each isotope will be computed automatically using an atomic database available from https://www.nist.gov/pml/atomic-weights-and-isotopic-compositions-relative-atomic-masses
  - this property does not determine which and how many isotopes will be included in the calculations (refer to the property **```nisotopes```** below) but only defines the isotopes by the symbols.
  - it should be an iterable of type list or tuple of strings and is not mandatory

In the following example the symbols for three of the isotopes of NiH molecule - <img src="/doc/tex/a01778636bb5ed1e24e91614821128df.svg?invert_in_darkmode&sanitize=true" align=middle width=43.15082309999999pt height=26.76175259999998pt/>, <img src="/doc/tex/f04c8e60a7fed2786b83de051f533185.svg?invert_in_darkmode&sanitize=true" align=middle width=43.15082309999999pt height=26.76175259999998pt/> and <img src="/doc/tex/bfb4aef268944889b9a31b9123d21404.svg?invert_in_darkmode&sanitize=true" align=middle width=43.15082309999999pt height=26.76175259999998pt/> are defined as: 

```python
# define the symbols for three isotopes
mdata.molecule = ['58Ni1H', '60Ni1H', '62Ni1H']

# and this also works
mdata.molecule = ['58 Ni 1 H', '60Ni 1H', '62 Ni 1H']
```

- **```masses```** - defines one or more isotopic forms of the same molecule by specifing the values of their reduced masses i.e. each defined item represents the reduced mass for each isotope computed by <img src="/doc/tex/a5f91a48b545c00d9b57f59247ed65e8.svg?invert_in_darkmode&sanitize=true" align=middle width=190.2203985pt height=24.65753399999998pt/> in amu units.
  - this property does not determine which and how many isotopes will be included in the calculations (refer to the property **```nisotopes```** below) but only defines their masses.
  - it should be an iterable of type list or tuple of float numbers and is not mandatory

```python
# define the reduced masses for three of NiH isotopes
mdata.masses = [0.990592988928, 0.991157254368, 0.991686280089]

# or from the masses of the separate atoms
mA, mB = 1.00782503223, 57.94453305
mdata.masses = [mA * mB / (mA + mB)]
```

> **_NOTE:_** At least one of the above two properties (**```molecule```** or **```masses```**) should be specified (i.e. one of them is mandatory).

- **```nisotopes```** - the isotopes which to be included in the compuations
  - each number corresponds to the index of an isotope in the list of defined isotopes; the counting starts from 1 up to the number of items in **```molecule```** or **```masses```**
  - it should be a list or tuple of integer numbers and is mandatory.

In this example the computations will be performed only for the second and the third isotope:
```python
# define the masses of 3 isotopes
mdata.molecule = ['58Ni1H', '60Ni1H', '62Ni1H']

# but compute only for the 2nd and the 3rd one
mdata.nisotopes = [2, 3]
```

<!-- omit in toc -->
### Values of the rotational quantum number J

- **```jrange```** - defines the interval of rotational qunatum numbers which to be used in the computations.
  - specify the initial and the final value of the required range of rotational quantum numbers
  - it is an iterable of type list or tuple containing only 2 elements (integer or float numbers)

```python
# the J values will be: 0, 1, 2, 3, 4 and 5
mdata.jrange = (0, 5)
```

- **```jvalues```** - defines a specific set of values for the rotational qunatum number
  - it could be a single number or list of numbers defined as type list/tuple of integer or float numbers

```python
# the J values will be: 7.5, 8.5, 11.5, 12.5
mdata.jvalues = (7.5, 8.5, 11.5, 12.5)
```

In the below example the values of the rotational quantum number that will be used in the calculations are: 0.5, 1.5, 2.5, 3.5, 4.5, 10.5, 12.5, 15.5, 16.5

```python
mdata.jrange = (0.5, 4.5)
mdata.jvalues = (10.5, 12.5, 15.5, 16.5)
```
> **_NOTE:_** At least one of the properties **```jvalues```** or **```jrange```** is mandatory i.e. either one or both should be set.

<!-- omit in toc -->
### Parity labels

- **```parities```** - parity levels which to be computed. 
  - we use a convention for designating the parity labels with numbers 1 and 0 for e and f respectivly.
  - it accepts a sinle value or tuple/list of one or two values of integer type (0 or 1).
  - not mandatory; if it is not set both e- and f-levels will be computed by default.

```python
# both e- and f-levels will be computed
mdata.parities = 0, 1

# only e-levels will be computed
mdata.parities = 1
```

<!-- omit in toc -->
### Reference level

- **```referenceJ```** - a reference J level whose energy will be used for shifting all remaining levels during the computations
  - the e/f levels and also the levels corresponding to different isotopes in general might be shifted with different energy
  - integer or float number

- **```referenceE```** - a reference energy value which will be used for shifting all levels i.e. all levels will be shifted with the same value
  - integer or float number

Examples:
```python
mdata.referenceJ = 2.5
```

```python
mdata.referenceE = 1000.
```
None of them is mandatory. **```referenceJ```** has higher proprity if both are specified simultaneously.

<!-- omit in toc -->
### Experimental Data

Experimental data can be provided if we call the method **```get_exp_data```**:

```python
mdata.get_exp_data('exp_KCs.dat', markers=[1,2,3,11,12])
```

The only required parameter is the name of the file containing the experimental data. **```markers```** is an optional parameter which allows the program to differentiate between the data of the isotopes. The markers are provided as a separate column in the experimental data file (see below).

---
<!-- omit in toc -->
#### **Structure of the experimental data file**

```python
298
      1       0     2.5        0.000000      1        0.0050   5        5   
      2       0     2.5        0.000000      0        0.0050   5        5   
      3       0     3.5       53.827000      1        0.0050   5        5   
      4       0     3.5       53.827000      0        0.0050   5        5   
      5       0     4.5      122.962000      1        0.0050   5        5   
      6       0     4.5      122.962000      0        0.0050   5        5   
      7       0     5.5      207.352000      1        0.0050   5        5   
      8       0     5.5      207.352000      0        0.0050   5        5   
      9       0     6.5      306.935000      1        0.0050   5        5   
     10       0     6.5      306.935000      0        0.0050   5        5
```

The first row is the number of data which should be included in the computations. The columns are as follows: counter, the vibrational quantum number, the rotational qunatum number, the experimental energy value, parity label (0 or 1), experimental uncertainty, marker and state. Levels with markers between 0 and 9 belong to the first isotope, levels with merkers between 10 and 19 belong to the second isotope and so on. The markers also provide convinient and easy way to temporary exculde some levels from the computations or group them by common characteristics.

---

## Grid Object Definition
We need to instanciate an object of type **```Grid```** in order to define the parameters releated to the mesh of points for the considered problem. Uniform and non-uniform type of grids are supported. The parameters that are needed to initilize the **```Grid```** object are:
- **```npoints```** - the number of grid points
    - It should be a positive integer number
  
- **```rgrid```** - the range of internuclear distancies .
    - specify the initial and the final values of the R grid.
    - a tuple/list containing two integer or float numbers

We may provide an additional (optional) parameter:
- **```solver```** - the method used for solution of the Schrodinger equation
    - The initialization of the grid depends on which solver method is selected.
    - It should be a string and the following values are available:
      - **```'sinc'```**
      - **```'fourier'```**
      - **```'fd5'```**

      If one of the values **```'sinc'```** or **```'fourier'```** is selected then FGH (Fourier Grid Hamiltonian) method will be used for the solution of the Schrodinger equation. If 'fd5' option is selected the Schrodinger equation will be solved by Finite Diffrence method.**```'sinc'```** is set by default option.

> **_NOTE:_** All parameters in program with values of type string are case-insensitive.

<!-- omit in toc -->
### Uniform grid
In order to generate an equdistant set of points **```npoints```** and **```rgrid```** are the only required parameters. 

This example shows how to initialize the grid object for uniform grid.
```python
# create a uniform grid of points
grid = Grid(npoints=100, rgrid=(1.0, 3.0))
```

<!-- omit in toc -->
### Nonuniform grid
For the generation of the nonuniform grid an analytical mapping procedure is applied. In this case **```alpha```** and **```rbar```** parameters should be set.

- **```alpha```** - the power parameter used in the analytical mapping formulas
  - It should be a positive integer or float number
  - If this parameter has a negative value then uniform grid without mapping will be generated. This allows one to easyly switch between uniform and nonuniform grid reperesentations.
- **```rbar```** - the value of this parameter is usually close to the internuclear distance

This example shows how the grid object should be initialized for nonuniform grid.

```python
# create a nonuniform grid of points
grid = Grid(npoints=100, rgrid=(1.0, 3.0), solver='sinc', alpha=3.0, rbar=1.458)
```
----
We can see how the generated grid of points looks by calling the function
**```get_grid_points```**
```python
print(grid.get_grid_points())
```
which will produce an output like

```python
[1.41729459 1.48936042 1.56142625 1.63349207 1.7055579  1.77762373
 1.84968955 1.92175538 1.99382121 2.06588703 2.13795286 2.21001869
 2.28208451 2.35415034 2.42621617 2.498282   2.57034782 2.64241365
 2.71447948 2.7865453  2.85861113 2.93067696 3.00274278 3.07480861
 3.14687444 3.21894026 3.29100609 3.36307192 3.43513774 3.50720357
 3.5792694  3.65133522 3.72340105 3.79546688 3.8675327  3.93959853
 4.01166436 4.08373018 4.15579601 4.22786184 4.29992766 4.37199349
 4.44405932 4.51612515 4.58819097 4.6602568  4.73232263 4.80438845]
```
or we can save the generated grid points in file by writing

```python
import numpy as np

np.savetxt('grid_points.dat', grid.rgrid)
```

## Channel Object Definition

A channel is defined as an electronic state with definite values of <img src="/doc/tex/b23332f99af850a48831f80dbf681ed6.svg?invert_in_darkmode&sanitize=true" align=middle width=11.41554479999999pt height=22.465723500000017pt/>, <img src="/doc/tex/e257acd1ccbe7fcb654708f1a866bfe9.svg?invert_in_darkmode&sanitize=true" align=middle width=11.027402099999989pt height=22.465723500000017pt/>, <img src="/doc/tex/813cd865c037c89fcdc609b25c465a05.svg?invert_in_darkmode&sanitize=true" align=middle width=11.87217899999999pt height=22.465723500000017pt/> and <img src="/doc/tex/9432d83304c1eb0dcb05f092d30a767f.svg?invert_in_darkmode&sanitize=true" align=middle width=11.87217899999999pt height=22.465723500000017pt/> quantum numbers. Each channel has to be defined as an object of type **```Channel```** and the parameters that are needed to instatniate objects of this type are:

- **```filep```** - the name of the file containing the parameters for the potential
  - It is a string parameter referring to an existing file or file path
  - The structure of the file depends on which model function is selected by the propery **```model```** described below.
- **```model```** - defines the potential energy function model (PEC model); it could be a pointwise function, builtin analytical potential function or custom analytical potential function (defined by the user). The possible values are:
  - **```'pointwise'```** : pointwise potential with cubic spline interpolation using the **```scipy```** builtin class <a href="https://yaml.org/" target="_blank">scipy.interpolate.CubicSpline</a>
  - **```'cspline'```** : pointwise potential with cubic spline interpolation using our own implementation
  - **```'Morse'```** : Morse potential
  - **```'EMO'```** : EMO (Expanded Morse Oscilator) potential
  - **```'MLR'```** : MLR (Morse/Long-Range) potential
  - **```'custom'```** : custom analytical potential
- **```nlambda```** - the quantum number <img src="/doc/tex/b23332f99af850a48831f80dbf681ed6.svg?invert_in_darkmode&sanitize=true" align=middle width=11.41554479999999pt height=22.465723500000017pt/> of the state
  - positive integer number
- **```sigma```** - the quantum number <img src="/doc/tex/813cd865c037c89fcdc609b25c465a05.svg?invert_in_darkmode&sanitize=true" align=middle width=11.87217899999999pt height=22.465723500000017pt/> of the state
  - integer or float number
- **```multiplicity```** - the multiplicity defined as <img src="/doc/tex/eb1977c2439009c2f58e150e1d2769c0.svg?invert_in_darkmode&sanitize=true" align=middle width=47.556992999999984pt height=22.465723500000017pt/> for the state
  - positive integer or float number
- **```rot_correction```** - correction to the diagonal rotational Hamiltonian
  - optional parameter of type integer or float


Here is an example definition of **```Channel```** object (for readability the parameters are written on separate lines):
```python
ch1 = Channel(
    filep='poten_morse.pot',
    model='morse',
    nlambda=0,
    sigma=0,
    multiplicity=1,
    rot_correction=0.
)
```

----
<!-- omit in toc -->
#### **Structure of the potential files**

- for pointwise potentials i.e. **```model='pointwise'```** or **```model='cspline'```** the potential file looks like

```python
13
    1.1168478261      14200.41785382920    0
    1.2146739130       8607.83894098507    0
    1.2880434783       5683.93800056556    0
    1.3858695652       2976.97391482846    0
    1.5081521739       1953.39920190690    0
    1.6304347826       2525.33711155168    0
    1.7527173913       4070.10001126557    0
    1.8505434783       5627.41358380562    0
    1.9728260870       7828.95851896513    0
    2.0859375000      10124.81139863041    0
    2.2265625000      12893.01493454805    0
    2.3671875000      15520.09623454644    0
    2.5078125000      17570.48597865039    0
```

It contains 3 columns. The first column determines the grid points, the second column are the values of the potential at these points, the third column is required when running the fitting procedure and tells whether the parameter is fixed or free. On the top of the file the number of points is specified.

- for Morse potential i.e. **```model='Morse'```**

```
Te = 1.00000000e-05  0
De = 1.00000000e+01  0
a  = 1.00000000      0
re = 1.27463182      0
```

- for EMO potential i.e. **```model='EMO'```**

```
Te =    1.00000   1
De =    1000.00   0
p  =    2         1
N  =    2         1
b0 =    1.00000   1 
b1 =    1.10000   1
b2 =    1.10000   1 
re =    1.27455   1
```

- for MLR potential i.e. **```model='MLR'```**

...
- for custom potential i.e. **```model='custom'```**

```
param1 = 4.5e-11   0
param2 = 4.5e-05   0
param3 = 9.9e-01   0
param4 = 2.4e+01   0
param5 = 8.4e+04   0
param6 = 5.1e+00   0
```

----

- **```custom_function```** - the name of some analytical function defined by the user

  - it is a python function defined in the current or in any other *.py file. In this case the parameter **```model```** has to be set to 'custom'.
  - it accepts exactly 2 input parameter and should return a single parameter. The first input parameter is an array containing the values of the potential parameters as defined in the potential file through the parameter **```filep```**. In this file the parameters have to be defined with the keyword 'param' followed by the number of corresponding parameter as is shown in the example below. The second input argument is an array containing the grid points. The output parameter has to be an array containing the calculated values of the potential function on the grid points. The length of the returned array should be equal to the number of grid points.

  - all input and output arguments are 1D **```numpy```** arrays of type **```'numpy.ndarray'```**.
  - it is an optional parameter

A simple example of a user defined function implementing the Morse potential is shown here:

```python
import numpy as np

def VMorse(params, rgrid):

    # unpacking the input parameters
    Te, De, a, re = params

    # returns the values of the Morse function over the grid
    return Te + De * np.power((1.0 - np.exp(-a*(rgrid-re))), 2.0)
```
In this case the parameters in the potential file are defined in the following way:

```
param1 = 4.55633525e-11   0
param2 = 4.55633525e-05   0
param3 = 9.99989591e-01   0
param4 = 2.40869887e+00   0
```
Note the special construction 'param+number' with conscuitive numbers starting from 1. In the input array the parameters will be aranged by this number.

> **_NOTE:_** Note that all parameters releated to the custom potential should be in au units!

Afterwards we need to create a list containing all  **```Channel```** obejcts which will be included the following compuations then call the function **```set_channel_parameters```** and pass the created list with channels as an argument to it:

```python
# combine the defined channels ch1 and ch2 in one list
channels = [ch1, ch2]

# then pass it to the function
Channel.set_channel_parameters(channels)
```

## Coupling Object Definition

The interactions between the channels are represented as objects of type **```Coupling```**. The parameters that we need to provide in order to initilize a **```Coupling```** object are:

- **```interact```** - the channels connected by this interaction
  - should be of type tuple or tuple of tuples
  - the numbers should correspond to the indicies of the channels in the channel list

- **```coupling```** - the type of the interation with the following available options:
  - **```spin-orbit```** : the diagonal and off-diagonal Spin-Orbit (SO) interaction
  - **```LJ```** : L-uncoupling interaction
  - **```SJ```** : Spin-uncoupling interaction
  - **```SL```** : Spin-electronic interaction
  - **```spin-rot```** : Spin-Rotation interaction
  - **```spin-spin```** : Spin-Spin interaction
  - **```LambdaDe```** : second-order <img src="/doc/tex/9432d83304c1eb0dcb05f092d30a767f.svg?invert_in_darkmode&sanitize=true" align=middle width=11.87217899999999pt height=22.465723500000017pt/> or <img src="/doc/tex/b23332f99af850a48831f80dbf681ed6.svg?invert_in_darkmode&sanitize=true" align=middle width=11.41554479999999pt height=22.465723500000017pt/> doubling effect on e-parity levels
  - **```LambdaDf```** : second-order <img src="/doc/tex/9432d83304c1eb0dcb05f092d30a767f.svg?invert_in_darkmode&sanitize=true" align=middle width=11.87217899999999pt height=22.465723500000017pt/> or <img src="/doc/tex/b23332f99af850a48831f80dbf681ed6.svg?invert_in_darkmode&sanitize=true" align=middle width=11.41554479999999pt height=22.465723500000017pt/> doubling effect on f-parity levels
  - **```LambdaD```** : second-order <img src="/doc/tex/9432d83304c1eb0dcb05f092d30a767f.svg?invert_in_darkmode&sanitize=true" align=middle width=11.87217899999999pt height=22.465723500000017pt/> or <img src="/doc/tex/b23332f99af850a48831f80dbf681ed6.svg?invert_in_darkmode&sanitize=true" align=middle width=11.41554479999999pt height=22.465723500000017pt/> doubling effect on both e- and f-parity levels

- **```model```** - the model of the coupling function wtih  possible values:
    - **```pointwise```**
    - **```cspline```**
    - **```custom```**

- **```multiplier```** - integer or float number which will multiply the defined function
  - if not provided the default value will be set to 1.0

- **```label```** - a label used to connect a certain coupling object and the corresponding parameters in the coupling file. 
  - should be of type string

All parameters are mandatory except multiplier.

```python
cp1 = diatom.Coupling(
    interact=(1,2),
    coupling='spin-orbit',
    model='pointwise',
    multiplier=1.0,
    label='cp1'
)
```
----

<!-- omit in toc -->
#### **Structure of the file with coupling parameters**

The couplings file uses a <a href="https://yaml.org/" target="_blank">YAML</a> syntax. For pointwise functions the file should consist of three columns. Before decalring each group of parameters the label to the corresponding **```Coupling```** object should be specified like:

```yaml
cp1:
- 0.750000000000      -606.45610689008595       0
- 1.250000000000      -593.55453482886696       0
- 1.500000000000      -599.03808570186402       0
- 2.000000000000      -606.92282658440399       0
- 2.500000000000      -603.70988681001700       0
- 3.000000000000      -602.90906506719705       0
- 5.000000000000      -603.00000000000000       0
cp2:
- 0.750000000000         0.02157527230695       0
- 1.500000000000         0.00110018854657       0
- 2.000000000000        -0.01268353489511       0
- 2.500000000000        -0.04232747723667       0
- 3.000000000000         0.00105314077299       0
- 5.000000000000         0.00001000000000       0
```

----
<!-- omit in toc -->
#### Shared parameters

When two or more interacting states coupled by the same or different interactions share the same set of parameters a more complicated construction of the **```Coupling```** object is possible.
This is what frequently happens in the case of states coupled by the <img src="/doc/tex/f7f72bf6b74988049786767c04a6bdf3.svg?invert_in_darkmode&sanitize=true" align=middle width=21.278616149999987pt height=22.465723500000017pt/> operator. Here is an example for the interaction <img src="/doc/tex/95986b9d0125352e11ec2ed1ae280905.svg?invert_in_darkmode&sanitize=true" align=middle width=19.70325884999999pt height=26.76175259999998pt/> ~ <img src="/doc/tex/b6ccdd799b87869964f3c3d42b645c41.svg?invert_in_darkmode&sanitize=true" align=middle width=21.07313174999999pt height=26.76175259999998pt/> in the lowest doublet states of NiH molecule:

```python
cp2 = diatom.Coupling(
    interact=((2,4), (3,5), (3,4)),
    coupling=('LJ', 'LJ', 'SL'),
    model='pointwise',
    multiplier=(2.0, 2.0, 2.0),
    label='cp2'
)
```

If we have defined the channels 2, 3, 4 and 5 as <img src="/doc/tex/431d13954ddb470c575ad747add2579f.svg?invert_in_darkmode&sanitize=true" align=middle width=39.54354689999999pt height=26.76175259999998pt/>, <img src="/doc/tex/9fa68b1d66e128b534bfe07e9855f1bf.svg?invert_in_darkmode&sanitize=true" align=middle width=39.54354689999999pt height=26.76175259999998pt/>, <img src="/doc/tex/f3523a0e299ea9fcab02794a2368f2ca.svg?invert_in_darkmode&sanitize=true" align=middle width=40.913419799999986pt height=26.76175259999998pt/> and <img src="/doc/tex/b54d59552cddec58d4e0da7405bb2952.svg?invert_in_darkmode&sanitize=true" align=middle width=40.913419799999986pt height=26.76175259999998pt/> then the pairs (2,4), (3,5) and (3,4) are connected by the same <img src="/doc/tex/f7f72bf6b74988049786767c04a6bdf3.svg?invert_in_darkmode&sanitize=true" align=middle width=21.278616149999987pt height=22.465723500000017pt/> operator in two different rotational interactions. Defined in this way they will use the same set of parameters - those labeled by 'cp2'. This type of definition is ....... 
<!-- We could define each of these pairs of interacting states as separate **```Coupling```** objects each having labels and the results will be the same. But then the fit will treat them ..... -->

## Molecule Levels Computation
<!-- omit in toc -->
### Calculate the eigenvalues and eigenvectors

To compute the eigenvalues an object of type **```MoleculeLevels```** should beinitialized with three required and one optional parameters as in the example:

```python
mlevels = MoleculeLevels(mdata, grid, channels, couplings=couplings)
```

The first parameter is the created [Molecule Data](#Molecule-Data-Definition) object, the second one is the created [Grid](#grid-definition) object and the third one is the created list with [Channel](#channel-definition) objects.
The optional parameter is created list of [Coupling](#coupling-definition) objects.

- **```eig_decomp```** - defines the  <a href="https://docs.scipy.org/" target="_blank">SciPy</a>  procedure which to be used for eigenvalues decomposition; the two possible values are:

  - **```'lapack'```** : calls the **```scipy```** builtin 
    <a href="https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.eigh.html" target="_blank">scipy.linalg.eigh</a> function to compute the eigenvalues and eigenvectors.

  - **```'arpack'```** : calls the **```scipy```** builtin
    <a href="https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.eigsh.html#scipy.sparse.linalg.eigsh" target="_blank">scipy.sparse.linalg.eigsh</a> function to compute the eigenvalues and eigenvectors.

  These two **```scipy```** procedures provide high-level interface to standard LAPACK and ARPACK routines written in Fortran and C. Default value is **```'lapack'```**

<!-- omit in toc -->
#### **Eigenvalue Decomposition with LAPACK**

This is the recommended choice when FGH is selected as a solving method.

If we need the eigenvalues and the eigenvectors only in a certain range of values or indicies then it is recommended to use one of the following two parameters whcih will likely reduce the computations time:

- **```energy_subset_index```** - defines the first and the last index of the subinterval of indicies for the returned eigenvalues

  - iterable of type list/tuple with 2 integer numbers defining the two indicies

    ```python
    mlevels.calculate_levels(energy_subset_index=(0, 6))
    ```

- **```energy_subset_value```** - defines the smallest and the largest value of subinterval of values for the returned eigenvalues

  - iterable if type list/tuple with 2 integer or float numbers defining the smallest and the largest value

    ```python
    mlevels.calculate_levels(energy_subset_value=(1000., 3000.))
    ```

  > **_NOTE:_**  If neither is used then all possible eigenvalues and thier eigenvectors will be computed and returned

- **```lapack_driver```** - the lapack routine computing the eigenvalues and eigenvectors; the possible values are:
    - **```'ev'```**  : calls dsyev routine which computes _all_ eigenvalues and eiegnvectors of a real symmetric matrix
    - **```'evd'```** : calls dsyevd routine which computes _all_ eigenvalues and eiegnvectors for real symmetric matrix using a divide and conquer algorithm
    - **```'evr'```** : calls dsyevr routine which computes a _selected_ eigenvalues and eigenvectors of a real symmetric matrix using the RRR algorithm
    - **```'evx'```** : calls dsyevx routine which computes a _selected_ eigenvalues and eiegnvectors of a real symmetric matrix

  The default value is set to **```'evr'```** which in general is the recommended choice. **```'evx'```** is faster when only a few of the eigenvalues are desired. **```energy_subset_index```** and **```energy_subset_value```** cannot be used together with **```'ev'```** and **```'evd'```** because they compute all eigenvalues.

    ```python
    mlevels.calculate_levels(eig_decomp='lapack', lap_driver='evx')
    ```

<!-- omit in toc -->
#### **Eigenvalue Decomposition with ARPACK**
ARPACK procedure is efficient and suitable for finding the _largest_ eigenvalues of a sparse matrix especially when a few of them are needed. If the _smallest_ eigenvalues need to be computed then it is recommended to use a shift-invert mode. In this mode the original eigenvalue problem will be transformed to an eqivalent problem so that the original small eigenvalues u will correspond to the transformed large eigenvalues v: v = 1 / u

For further information: https://docs.scipy.org/doc/scipy/reference/tutorial/arpack.html

The releated parameters are:
- **```arpack_k```** - the dsired number of eigenvalues and eigenvectors
  - it is a positive integer number and should be smaller than the dimension of the matrix 
- **```arpack_which```** - which eigenvalues to compute
  - The possible values are: **```'LM'```** (largest in magnitude), **```'SM'```** (smallest in magnitude), **```'LA'```** (largest algebraic), **```'SA'```** (smallest algebraic), **```'BE'```** (half (k/2) from each end of the spectrum).
- **```arpack_sigma```** - if this parameter is specified a shift-invert mode is applied. Then the procedure will return k (specified by **```arpack_k```**) shifted eigenvalues which have values around the specified value of this parameter (integer or float number).
  
  - efficient when the smallest eigenvalues are required (**```arpack_which```**=**```'SM'```**)

```python
mlevels.calculate_levels(eig_decomp='arpack', arpack_k=25, arpack_which='SM', arpack_sigma=0.1)
```
> **_NOTE:_** With ARPACK it is not possible to compute all eigenvectors of a matrix.

> **_NOTE:_** The parameters **```lapack_driver```**, **```subset_by_index```** and **```subset_by_value```** are not relevent when **```eig_decomp```**='arpack'.
Similarly **```arpack_k```**, **```arpack_which```** and  **```arpack_sigma```** are not relevent when **```eig_decomp```**='lapack'

<!-- omit in toc -->
### Identification of computed levels
- **```identify```** - defines the method for identification of the calculated energies i.e. the procedure used for the assignment of vibrational quantum number for each level and the determination of the state which it belongs to based on the available experimental data.
  - integer number with 2 possible values:
    - **```0```** : identification by energy
    - **```1```** : identification by the largest value of the coupling coefficient
    <!-- - _2_ : a slow version of _0_; defined and used only for debugging puropse
    - _3_ : a slow version of _1_; defined and used only for debugging puropse -->
  - defult value is set to **```0```**
  - it is not relevent when experimental data are not provided

<!-- omit in toc -->
### Output Format and Storing Options

- **```store_predicted```** - whether to save _all_ computed eigenvalues in file
  - the predicted  eigenvalues will be written in file called 'evalues_predicted.dat'. The name of this file can be changed if additional parameter **```predicted_file```** is added.

The file looks like:

```python
#     No  v      Ecalc          J   parity  marker   CC1   state   lambda  omega
      1   0      56.113986     1.0     0      1     1.000      1      0     0.00
      2   1     180.940849     1.0     0      1     1.000      1      0     0.00
      3   2     387.925302     1.0     0      1     1.000      1      0     0.00
      4   3     677.259248     1.0     0      1     1.000      1      0     0.00
      5   4    1049.088043     1.0     0      1     1.000      1      0     0.00
      6   5    1503.471257     1.0     0      1     1.000      1      0     0.00
```
  
- **```store_info```** - save more detailed information regarding the computations
  
- **```store_evecs```** - save the computed eigenvectors in files

- **```sort_output```** - sort the selected eigenvalues by provided column numbers in order

- **```sort_predicted```** - sort all computed eigenvalues by provided column numbers in order

## Examples

<!-- omit in toc -->
### The minimal working example: Morse potential for <img src="/doc/tex/67ee5b6741d3b4139212bbf0588a9e3f.svg?invert_in_darkmode&sanitize=true" align=middle width=18.881345999999994pt height=22.465723500000017pt/>

This is an example with the minimal set of nessecary parameters for running a single channel computations with the Morse potential for <img src="/doc/tex/67ee5b6741d3b4139212bbf0588a9e3f.svg?invert_in_darkmode&sanitize=true" align=middle width=18.881345999999994pt height=22.465723500000017pt/> molecule.

It is done with only a few lines of code.

```python
#!/usr/bin/env python

from diatom import *

mdata = MoleculeData()
mdata.molecule = ['1H1H']
mdata.nisotopes = [1]
mdata.jrange = (0, 1)

grid = Grid(npoints=170, rgrid=(0.3, 2.5))

ch1 = Channel(
    filep='morse_H2.pot',
    model='morse',
    nlambda=0,
    sigma=0,
    multiplicity=1,
)

channels = [ch1]
Channel.set_channel_parameters(channels)

mlevels = MoleculeLevels(mdata, grid, channels)
mlevels.calculate_levels()
```
The content of the potential file 'morse_H2.pot' is the following:
```python
Te =   0.0000000          0
De =   3.8276634e+04      0
a  =   1.9419600          0
re =   0.7419100          0
```
Here all eigenvalues and eigenvectors for J=0 and J=1 with both e/f parity labels will be computed with using the 'sinc' method. The computed eigenvalues will be stored in file called 'evalues_predicted.dat'.
The first few lines of this file look like

```python
#     No  v      Ecalc          J   parity  marker   CC1    state  lambda   omega
      1   0    2165.954681     0.0     0      1     1.000      1      0     0.00
      2   1    6308.624389     0.0     0      1     1.000      1      0     0.00
      3   2   10198.976338     0.0     0      1     1.000      1      0     0.00
      4   3   13837.015745     0.0     0      1     1.000      1      0     0.00
      5   4   17222.753191     0.0     0      1     1.000      1      0     0.00
      6   5   20356.205017     0.0     0      1     1.000      1      0     0.00
      7   6   23237.390221     0.0     0      1     1.000      1      0     0.00
      8   7   25866.323436     0.0     0      1     1.000      1      0     0.00
```

Then we can plot the potential:

```python
plot = diatom.Plotting()
plot.plot_potentials_on_grid(mlevels, show=True, fformat='svg')
```

![Morse potential H2](./plotting/potential_morse_H2.svg)

<!-- omit in toc -->
### Another example: KCs ...



<!-- omit in toc -->
### A more complex example: pointwise potentials for <img src="/doc/tex/8207f29c059ad1741b922ed20deb57b7.svg?invert_in_darkmode&sanitize=true" align=middle width=29.338013099999987pt height=26.76175259999998pt/>, <img src="/doc/tex/8cbc30f07cc338dcb4d7b09bb011e028.svg?invert_in_darkmode&sanitize=true" align=middle width=19.70325884999999pt height=26.76175259999998pt/> and <img src="/doc/tex/b6ccdd799b87869964f3c3d42b645c41.svg?invert_in_darkmode&sanitize=true" align=middle width=21.07313174999999pt height=26.76175259999998pt/> states of NiH

The complex of the three doublet electronic states <img src="/doc/tex/8207f29c059ad1741b922ed20deb57b7.svg?invert_in_darkmode&sanitize=true" align=middle width=29.338013099999987pt height=26.76175259999998pt/>, <img src="/doc/tex/8cbc30f07cc338dcb4d7b09bb011e028.svg?invert_in_darkmode&sanitize=true" align=middle width=19.70325884999999pt height=26.76175259999998pt/> and <img src="/doc/tex/b6ccdd799b87869964f3c3d42b645c41.svg?invert_in_darkmode&sanitize=true" align=middle width=21.07313174999999pt height=26.76175259999998pt/> for NiH molecule is an example of strongly perturbed and coupled states by various interactions. The three electronic states are described by five channels and in this example only five of all interactions are accounted for.

```python
#!/usr/bin/env python

from diatom import *

mdata = MoleculeData()
mdata.molecule = ['58Ni1H']
mdata.nisotopes = [1]
mdata.jrange = (0.5, 12.5)
mdata.referencej = 2.5
mdata.parities = 0, 1
mdata.set_exp_data('nih_exp.dat', markers=[5])

grid = Grid(npoints=170, rgrid=(0.75, 3.0), solver='sinc')

vpot = 'v_nih.pot'
wpot = 'w_nih.pot'
xpot = 'x_nih.pot'

ch1 = Channel(
    filep=vpot,
    model='pointwise',
    nlambda=0,
    sigma=0.5,
    multiplicity=2
)

ch2 = Channel(
    filep=wpot,
    model='pointwise',
    nlambda=1,
    sigma=-0.5,
    multiplicity=2
)

ch3 = Channel(
    filep=wpot,
    model='pointwise',
    nlambda=1,
    sigma=0.5,
    multiplicity=2
)

ch4 = Channel(
    filep=xpot,
    model='pointwise',
    nlambda=2,
    sigma=-0.5,
    multiplicity=2
)

ch5 = Channel(
    filep=xpot,
    model='pointwise',
    nlambda=2,
    sigma=0.5,
    multiplicity=2
)

channels = [ch1, ch2, ch3, ch4, ch5]
Channel.set_channel_parameters(channels)

cp1 = Coupling(
    interact=((2, 2), (3, 3)),
    coupling=('spin-orbit', 'spin-orbit'),
    model='pointwise',
    multiplier=(-1.0, 1.0),
    label='SO_WW'
)

cp2 = Coupling(
    interact=((4, 4), (5, 5)),
    coupling=('spin-orbit', 'spin-orbit'),
    model='pointwise',
    multiplier=(-1.0, 1.0),
    label='SO_XX'
)

cp3 = Coupling(
    interact=(1, 2),
    coupling='spin-orbit',
    model='pointwise',
    multiplier=1.0,
    label='SO_VW'
)

cp4 = Coupling(
    interact=(3, 4),
    coupling='spin-orbit',
    model='pointwise',
    multiplier=1.0,
    label='SO_VX'
)

cp5 = Coupling(
    interact=((2, 4), (3, 5), (3, 4)),
    coupling=('LJ', 'LJ', 'SL'),
    model='pointwise',
    multiplier=(2.0, 2.0, 2.0),
    label='LJ_WX'
)

couplings = [cp1, cp2, cp3, cp4, cp5]
Coupling.set_coupling_parameters('couplings.dat', couplings)

mlevels = MoleculeLevels(mdata, grid, channels, couplings=couplings)
mlevels.calculate_levels(energy_subset_value=(0, 7000.), identify=0)
```

The couplings file looks like this
```yaml
SO_WW:
- 0.750000000000      -300.71901701285799       0
- 1.250000000000      -300.59347903145402       0
- 1.500000000000      -295.33299681003803       0
- 2.000000000000      -304.13289200589099       0
- 2.500000000000      -301.78320985488801       0
- 3.000000000000      -301.49811571690401       0
SO_XX:
- 0.750000000000      -606.45610689008595       0
- 1.250000000000      -593.55453482886696       0
- 1.500000000000      -599.03808570186402       0
- 2.000000000000      -606.92282658440399       0
- 2.500000000000      -603.70988681001700       0
- 3.000000000000      -602.90906506719705       0
SO_VW:
- 0.750000000000      -620.71974693075003       0
- 1.500000000000      -645.86464080798601       0
- 2.000000000000      -716.97127978362300       0
- 2.500000000000      -735.96722850460196       0
- 3.000000000000      -740.41087084696699       0
SO_VX:
- 0.750000000000      -590.84142865184197       0
- 1.500000000000      -600.45207379396902       0
- 2.000000000000      -606.05480585440398       0
- 2.500000000000      -605.57223541150302       0
- 3.000000000000      -603.00000000000000       0
LJ_WX:
- 0.750000000000         1.23969597272800       0
- 1.200000000000         0.94539092137400       0
- 1.500000000000         0.99276646115100       0
- 2.000000000000         1.23560350766500       0
- 5.000000000000         1.00000000000000       0
```

The output looks like:

```python
#    No       v       J   omega   sigma  lambda  parity  marker   Ecalc           Eexp          delta        unc       CC1      CC2      CC3      CC4      CC5      state
      1       0     2.5     2.5     0.5      2      0      5        0.000000       0.000000     0.000000     0.0050    0.000    0.000    0.000    0.000    1.000     5
      2       0     2.5     2.5     0.5      2      1      5        0.000000       0.000000     0.000000     0.0050    0.000    0.000    0.000    0.000    1.000     5
      3       0     3.5     2.5     0.5      2      0      5       54.246138      53.827000     0.419138     0.0050    0.000    0.000    0.001    0.000    0.999     5
      4       0     3.5     2.5     0.5      2      1      5       54.246138      53.827000     0.419138     0.0050    0.000    0.000    0.001    0.000    0.999     5
      5       0     4.5     2.5     0.5      2      0      5      123.918310     122.962000     0.956310     0.0050    0.000    0.000    0.001    0.000    0.999     5
      6       0     4.5     2.5     0.5      2      1      5      123.918310     122.962000     0.956310     0.0050    0.000    0.000    0.001    0.000    0.999     5
```

We can use the program to further analyize and plot the experimental and computed data.

```python
import numpy as np
import matplotlib.pyplot as plt

# this function returns the experimental data as numpy array
exp_data = mdata.get_exp_data()

# then write some code to plot them as function of J

fdata = exp_data[exp_data[:, 4] == 0]
edata = exp_data[exp_data[:, 4] == 1]

# plot exp data on two plots for e and f
_, ax = plt.subplots(1, 2, sharex=True, sharey=True, figsize=(11, 5))

ax[0].plot(fdata[:, 2], fdata[:, 3], 'bo', fillstyle='none')
ax[0].set_xlabel('J')
ax.set_ylabel(r'Energy (cm<img src="/doc/tex/db982724014b78c648c5d008c86b3d09.svg?invert_in_darkmode&sanitize=true" align=middle width=16.82656799999999pt height=26.76175259999998pt/>)')

ax[1].plot(edata[:, 2], edata[:, 3], 'ro', fillstyle='none')
ax[1].set_xlabel('J')
plt.show()

# or plot all of them in one plot for each vibrational level and state 
_, ax1 = plt.subplots(figsize=(9, 6))
vs = np.unique(exp_data[:, 1])
ss = np.unique(exp_data[:, -1])

for v in vs:
    for s in ss:
        f = fdata[(fdata[:, 1] == v) & (fdata[:, -1] == s)]
        e = edata[(edata[:, 1] == v) & (edata[:, -1] == s)]
        ax1.plot(e[:, 2], e[:, 3], color='blue', marker='o', fillstyle='none', markersize=8, linewidth=0.2)
        ax1.plot(f[:, 2], f[:, 3], color='green', marker='X', fillstyle='none', markersize=6, linewidth=0.2)

ax1.set_xlabel('J')
ax1.set_ylabel(r'Energy (cm<img src="/doc/tex/db982724014b78c648c5d008c86b3d09.svg?invert_in_darkmode&sanitize=true" align=middle width=16.82656799999999pt height=26.76175259999998pt/>)')
plt.show()
```

![Experimental data plot](./plotting/exp_date_ef.svg)
![Experimental data plot](./plotting/exp_date_ef_vs.svg)

# Fitting of the Calculated Energy Levels

The **```DiAtom```** module has several implemented procedures for weighted least-squares fitting of the calculated to the experimental energy levels. In all cases what we want to minimize is the difference between the experimental (observed) energies and the calculated ones. Therefore we define the <img src="/doc/tex/a67d576e7d59b991dd010277c7351ae0.svg?invert_in_darkmode&sanitize=true" align=middle width=16.837900199999993pt height=26.76175259999998pt/> function as:

<p align="center"><img src="/doc/tex/8ec442e942723e367a844cc07ec2a05b.svg?invert_in_darkmode&sanitize=true" align=middle width=200.27529884999998pt height=47.8235406pt/></p>

where the computed energies are functions of the parameters that we would like to determine by minimizing the <img src="/doc/tex/a67d576e7d59b991dd010277c7351ae0.svg?invert_in_darkmode&sanitize=true" align=middle width=16.837900199999993pt height=26.76175259999998pt/> value. The dependance <img src="/doc/tex/4ca1aee6e8270689594bcea797736c2f.svg?invert_in_darkmode&sanitize=true" align=middle width=53.89543995pt height=27.91243950000002pt/> is in general nonlinear therefore an iterative procedure will be applied - starting from some trial values of the parameters a corrections will be generated and added to the current values on each iteration whcih will improve the <img src="/doc/tex/a67d576e7d59b991dd010277c7351ae0.svg?invert_in_darkmode&sanitize=true" align=middle width=16.837900199999993pt height=26.76175259999998pt/> value. The corrections will be found by solving the system:

<p align="center"><img src="/doc/tex/3d56e9c278a94b22d5667eb320678c6a.svg?invert_in_darkmode&sanitize=true" align=middle width=328.22589855pt height=47.2889208pt/></p>

where we have approximated the dependance <img src="/doc/tex/4ca1aee6e8270689594bcea797736c2f.svg?invert_in_darkmode&sanitize=true" align=middle width=53.89543995pt height=27.91243950000002pt/> by the first two terms in its Taylor expansion around <img src="/doc/tex/14fc0e68e921111d80f847f49f14fa1b.svg?invert_in_darkmode&sanitize=true" align=middle width=26.803695599999987pt height=29.190975000000005pt/>. This is a linear system of n equations with m unknowns (usually n > m) in the form <img src="/doc/tex/5e1e7c803719b446197e2edde2f0815c.svg?invert_in_darkmode&sanitize=true" align=middle width=52.89932669999998pt height=31.141535699999984pt/> where the unknown vector x is the vector with the corrections <img src="/doc/tex/3919bbc84b8079e27194efe99a1f6a80.svg?invert_in_darkmode&sanitize=true" align=middle width=23.09366069999999pt height=22.465723500000017pt/>, the right-hand side vector b is <img src="/doc/tex/f2e84668ab502e9be17614dcee78f4aa.svg?invert_in_darkmode&sanitize=true" align=middle width=124.01297039999999pt height=29.190975000000005pt/>, and the coefficient matrix A is formed by the first derivatives of the energies with respect to the parameters. 
<!-- The overall goal of the fit could be summirzied as:

<p align="center"><img src="/doc/tex/c57fae4c7444bddea5a45372a6b5fb2b.svg?invert_in_darkmode&sanitize=true" align=middle width=164.16738195pt height=32.990165999999995pt/></p> -->

As a first step we need to initialize the **```Fitting```** object for example like
```python
fit = diatom.Fitting(mlevels, progress=False)
```

The first parameter is the created **```MoleculeLevels```** object and the second parameter **```progress```** is optional and specifies whether to print some detailed output after the end of _every_ iteration. The default is **```False```** which means that a detailed output will only be printed after the end of the _last_ iteration.

## SVD Fit

In general it is not recommended to solve the above linear system by the method of the normal equations (that uses the matrix inverse) since the matrix A is singular mainly because of the following problem. Sometimes there exist two or more linear combinations of functions with the fitted parameters that can be added to the model functions without changing the <img src="/doc/tex/a67d576e7d59b991dd010277c7351ae0.svg?invert_in_darkmode&sanitize=true" align=middle width=16.837900199999993pt height=26.76175259999998pt/> value. This is an indication of a linear dependance between the model functions (the data matrix will be singular) and also means that there exist two or more sets of parameters that fit the data equally well. In this cases it is recommended to use the Singular Value Decomposition (SVD). In SVD the matrix <img src="/doc/tex/53d147e7f3fe6e47ee05b88b166bd3f6.svg?invert_in_darkmode&sanitize=true" align=middle width=12.32879834999999pt height=22.465723500000017pt/> (n x m) is represented as a product of three matrices <img src="/doc/tex/348ca4b32b3c0af83c00431b016dfad9.svg?invert_in_darkmode&sanitize=true" align=middle width=78.404007pt height=27.91243950000002pt/>, two unitary (or orthogonal in the real case) matrices U (n x n) and V (m x m) and one diagonal matrix <img src="/doc/tex/813cd865c037c89fcdc609b25c465a05.svg?invert_in_darkmode&sanitize=true" align=middle width=11.87217899999999pt height=22.465723500000017pt/> (n x m). This is known as full SVD. When <img src="/doc/tex/180688a8c2192e65c9812628d22d321f.svg?invert_in_darkmode&sanitize=true" align=middle width=46.21760714999999pt height=20.908638300000003pt/> (more data than parameters), <img src="/doc/tex/813cd865c037c89fcdc609b25c465a05.svg?invert_in_darkmode&sanitize=true" align=middle width=11.87217899999999pt height=22.465723500000017pt/> will have at most m nonzero rows and more compact representation is possible: <img src="/doc/tex/3652fbbfd63546aca5a4a0537748d81f.svg?invert_in_darkmode&sanitize=true" align=middle width=78.40398554999999pt height=31.141535699999984pt/> where <img src="/doc/tex/c827023fcc3b2c873aae402363f41328.svg?invert_in_darkmode&sanitize=true" align=middle width=13.01596064999999pt height=31.141535699999984pt/> is (n x m) submatrix of U and <img src="/doc/tex/d4261132636819ae7c6f4039dafc4016.svg?invert_in_darkmode&sanitize=true" align=middle width=11.87217899999999pt height=31.141535699999984pt/> is the (m x m) submatrix of <img src="/doc/tex/813cd865c037c89fcdc609b25c465a05.svg?invert_in_darkmode&sanitize=true" align=middle width=11.87217899999999pt height=22.465723500000017pt/>. This is known as "economy" SVD. The U and V matrices are called right and left singular vectors and the diagonal elments of <img src="/doc/tex/813cd865c037c89fcdc609b25c465a05.svg?invert_in_darkmode&sanitize=true" align=middle width=11.87217899999999pt height=22.465723500000017pt/> are called singular values. It is important that the singular values are hierarchically aranged from the largest to the smallest i.e. <img src="/doc/tex/73f03081ed485c91b9becde40212c0a5.svg?invert_in_darkmode&sanitize=true" align=middle width=139.52397194999998pt height=20.908638300000003pt/>...

<!-- The singular values are different from zero when the model functions are linearly independant.

SVD is very special matrix factorization because it can be applied to _any_ matrix, it is unique and guaranteed to exist. -->

The **```run_svd```** method has only default parameters:

- **```niter```** - the number of iterations. Default is **```niter=1```**

- **```deriv```** - the method used for computation of derivatives. Default is **```deriv='n'```**. When **```deriv='n'```** the derivatives will be computed numerically with the finite diffrence method. Another possibility is **```deriv='a'```** in which case the derivatives of the energies with respect to the parameters will be computed analitically with the Hellman-Feynman theorem:

<p align="center"><img src="/doc/tex/c3694908f44be879f0a2f10bc915062f.svg?invert_in_darkmode&sanitize=true" align=middle width=122.68708484999998pt height=33.81208709999999pt/></p>

- **```tol```** - the tolerance value. It determines which and how many linear combinations of the fitted parameters will be discarded because the matrix is singular. The singular values which are less than **```tol```** times the largest singular value are treated as zero. The rank of the matrix is determined by the number of the nonzero singular values and **```tol```** governs the effective rank of the matrix that is the number of singular values smaller than some specific value. Default is **```tol=0.1```**

- **```lapack_driver```** - The name of the LAPACK routine used to solve the least-squares problem. The three possible values **```'gelsd'```**, **```'gelsy'```** and **```'gelss'```** correspond to the names of these routines. **```'gelsd'```** uses singular value decomposition and a divide and conquer method, **```'gelsy'```** uses the orthogonal QR factorization of the matrix and **```'gelss'```** uses the singular value decomposition of the matrix. In the most cases **```'gelsd'```** will likely be the most efficient method. Default is **```lapack_driver='gelsd'```**.

- **```step_size```** - used to determine the change in the parameters during the computation of derivatives. Default is **```step_size=1.0e-4```**.

- **```is_weighted```** - whether to apply a weighted least-squares fitting with the method proposed by J. Watson (default is **```False```**) in which case the weights <img src="/doc/tex/c5e40d7bbcc7d29edc388de6341b5f66.svg?invert_in_darkmode&sanitize=true" align=middle width=26.80945079999999pt height=28.894955100000008pt/> given by the experimental uncerteinty will be replaced by the expression:

<p align="center"><img src="/doc/tex/fecf078bd951ed0d732f61d5030effdc.svg?invert_in_darkmode&sanitize=true" align=middle width=170.20721849999998pt height=38.51761815pt/></p>

- **```restart```** - not yet implemented
- **```limit```** - not yet implemented
- **```regular```** - not yet implemented


To find the least-squares solution **```run_svd```** calls <a href="https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.lstsq.html" target="_blank">scipy.linalg.lstsq</a> function. In the most widespread numerical libraries and packages the implementations which use SVD are based on the LAPACK implementation in Fortran - the routine called DGESVD.

## Minuit Fit

**```DiAtom```** provides an option to call the Minuit library in C++ via the Python frontend **```iminuit```**. It works with the Minuit Migrad subroutine for local function minimization. In addition to that it is able to estimate the uncertainty in the fitted parameters by computing the covariance and correlation matrices using two different algorithms and procedures called Hesse and Minos.

To run the Minuit Fit we should call the method **```run_minuit```** through the created **```Fitting```** object. This method has only optional parameters:

- **```niter```** - The number of iterations. Default is **```niter=0```**.
- **```step_size```** -  used to determine the initial gradient. Default is **```step_size=1.0e-4```**
- **```uncert```** - whether to compute and print the covariance and correlation matricies and the uncertainties in the in the model parameters after the fit completed. Default is **```False```**.

Example

```python
fit.run_minuit(niter=5, step_size=1.0e-3, uncert=True)
```

## Levenberg-Marquard Fit

Not yet implemented

# Computing the Transition Frequencies and Intensities

The transition frequencies between the rovibrational levels of two electronic states can be computed in two cases:
1. both states are represented by channels (as **```Channel```** objects) i.e. by potential curves 
2. both states are represented by their term values.

In either case we need first to define an object of type **```Spectrum```**:

```python
spec = diatom.Spectrum()
```

## States represented by channels

In this case we should call the method **```calculate_frequencies_by_states```** which has two required and many optional parameters.

```python
spec.calculate_frequencies_by_states(uch=1, lch=2)
```

- **```uch```** - the number of the channel corresponding to the upper electronic state
- **```lch```** - the number of the channel corresponding to the lower electronic state

All optional parameters are listed in the section below.

## States represented by term values

In this case we should call the method **```calculate_frequencies_by_term_values```** which has two required and many optional parameters.

```python
spec.calculate_frequencies_by_term_values(uch=1, lch=2)
```

- **```uterms```** - the name of the file containing the term values for the upper electronic state
- **```lterms```** - the name of the file containing the term values for the lower electronic state

-----

<!-- omit in toc -->
### Optional parameters

Both methods **```calculate_frequencies_by_states```**  and **```calculate_frequencies_by_term_values```** have the same set of optional parameters which are:


## Line strength

### Honl-London Factors

### Frank-Condon Factors

## Einstein A coefficient and Radiative Lifetime

# Plotting

The program provides some basic and simple plotting functions.

![Hamiltonian matrix colormesh](./plotting/hcolormesh.png)

<!-- ## Cross-section -->
