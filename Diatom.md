
- [**Diatom** module: what it is useful for and how to install and setup?](#diatom-module-what-it-is-useful-for-and-how-to-install-and-setup)
- [Diatomic molecule: basic theoretical concepts](#diatomic-molecule-basic-theoretical-concepts)
  - [The total Hamiltonian and basis functions](#the-total-hamiltonian-and-basis-functions)
  - [The Scrodinger equation for a single state and coupled system of states](#the-scrodinger-equation-for-a-single-state-and-coupled-system-of-states)
  - [The interaction terms and their matrix elements](#the-interaction-terms-and-their-matrix-elements)
  - [Methods for solving the Schrodinger equation](#methods-for-solving-the-schrodinger-equation)
  - [Potential Energy function models (PECs models)](#potential-energy-function-models-pecs-models)
- [Computing Energy Eigenvalues](#computing-energy-eigenvalues)
  - [Molecule Data Definition](#molecule-data-definition)
  - [Grid Definition](#grid-definition)
  - [Channel Definition](#channel-definition)
      - [The structure of the potential files](#the-structure-of-the-potential-files)
  - [Coupling Definition](#coupling-definition)
      - [Structure of the couplings file](#structure-of-the-couplings-file)
  - [Experimental Data](#experimental-data)
  - [Molecule Levels Computation](#molecule-levels-computation)
  - [The minimal working example](#the-minimal-working-example)
- [Fitting of Calculated Energy Levels](#fitting-of-calculated-energy-levels)
- [Plotting](#plotting)
- [Computing the Transition Frequencies](#computing-the-transition-frequencies)
- [Computing the Transition Intensities](#computing-the-transition-intensities)
  - [Einstein A coefficient and Radiative Lifetime](#einstein-a-coefficient-and-radiative-lifetime)
  - [Line strength](#line-strength)


# **Diatom** module: what it is useful for and how to install and setup?
The Python package **Diatom** allows various calculations for diatomic molecules to be performed. It supports single and coupled channels computations of bound rovibrational levels, intensity calculations, fitting to the experimental data.
The current functionality covered by the program includes:
* ..
* ..

Diatom package can be installed from the Python software repository PyPI (Python Package Index) via pip. From Linux command line execute

```
<img src="/tex/5b0e3276242a118e0cba11074b7127e8.svg?invert_in_darkmode&sanitize=true" align=middle width=700.2744672pt height=315.6164385pt/> chmod u+x main.py
<img src="/tex/b27545348bf905b65b0ddfb427bda915.svg?invert_in_darkmode&sanitize=true" align=middle width=1068.53799525pt height=276.1643841pt/><img src="/tex/21bfd3a02c5cb27c5e3ada6162fedac5.svg?invert_in_darkmode&sanitize=true" align=middle width=381.66411524999995pt height=24.65753399999998pt/><img src="/tex/7b342c84e914e939ac4cbd0eaaca7459.svg?invert_in_darkmode&sanitize=true" align=middle width=700.50290805pt height=203.6529759pt/><img src="/tex/797ffc6d128d85b2c09b9c23eb8e0986.svg?invert_in_darkmode&sanitize=true" align=middle width=406.34181885pt height=37.80850590000001pt/><img src="/tex/37ca9269559fef7d7504c55ca37c7020.svg?invert_in_darkmode&sanitize=true" align=middle width=268.2199014pt height=39.45205439999997pt/>R<img src="/tex/be0f0cee207d363be76cfc4530d0c67e.svg?invert_in_darkmode&sanitize=true" align=middle width=275.19620369999996pt height=22.831056599999986pt/>U(R) + ({\hbar^2}/{2\mu R^2})J(J+1)<img src="/tex/8d46e98288f18431d642a331d196bc77.svg?invert_in_darkmode&sanitize=true" align=middle width=578.6930638499999pt height=22.831056599999986pt/>\mu = M_{1} M_{2} /(M_{1} + M_{2})<img src="/tex/1da9fcafafb5c088ab6ddc5e21359621.svg?invert_in_darkmode&sanitize=true" align=middle width=44.86319144999999pt height=22.831056599999986pt/>M_1<img src="/tex/fd92a53167b3c6ae9574071613d555dc.svg?invert_in_darkmode&sanitize=true" align=middle width=27.11199479999999pt height=22.831056599999986pt/>M_2<img src="/tex/456a7ae47fb0b9caedf5550d4ec541a2.svg?invert_in_darkmode&sanitize=true" align=middle width=467.6466399pt height=22.831056599999986pt/>E_{vJ}<img src="/tex/a6c401b133327749f2d1d6015bf6e66e.svg?invert_in_darkmode&sanitize=true" align=middle width=390.9370575pt height=22.831056599999986pt/>\phi_{vJ}<img src="/tex/1e10cdb38b6f2980f14cb619e6afa95b.svg?invert_in_darkmode&sanitize=true" align=middle width=700.2744622499999pt height=322.0091391pt/><img src="/tex/14189077a38971377cea77bfabc6e843.svg?invert_in_darkmode&sanitize=true" align=middle width=217.99017239999998pt height=29.190975000000005pt/><img src="/tex/24844e386612a660188bbd69111bd424.svg?invert_in_darkmode&sanitize=true" align=middle width=246.4918137pt height=22.831056599999986pt/>D_{e}<img src="/tex/638b242b5112d037b2148b9e4ffead2d.svg?invert_in_darkmode&sanitize=true" align=middle width=824.2124484pt height=24.65753399999998pt/>\beta<img src="/tex/9acfbef9f5bc0134997832d6d47f3218.svg?invert_in_darkmode&sanitize=true" align=middle width=112.25229014999998pt height=22.831056599999986pt/>r_{e}<img src="/tex/aa4b599d89fdfa80139ea0e6db0642bd.svg?invert_in_darkmode&sanitize=true" align=middle width=327.67215855pt height=78.90410880000002pt/><img src="/tex/be22f48b5fe0e3c05c068cf6bfb428b7.svg?invert_in_darkmode&sanitize=true" align=middle width=383.94854189999995pt height=26.76175259999998pt/><img src="/tex/395c22333c4cfc54ed09e349f7a123b4.svg?invert_in_darkmode&sanitize=true" align=middle width=640.4583355499999pt height=45.84475499999998pt/><img src="/tex/07ddc82faad966e1ca7aa09c751fcd4a.svg?invert_in_darkmode&sanitize=true" align=middle width=192.07546875pt height=32.256008400000006pt/><img src="/tex/9a3645a25abf72e32fb414d1722b2cfa.svg?invert_in_darkmode&sanitize=true" align=middle width=586.57693875pt height=45.84475499999998pt/><img src="/tex/490b8cc557b393901eddeeccbbc02ee3.svg?invert_in_darkmode&sanitize=true" align=middle width=101.93435999999998pt height=33.76610489999999pt/><img src="/tex/4238e787c12f2e3e6690ef11ed1486f1.svg?invert_in_darkmode&sanitize=true" align=middle width=700.2747410999999pt height=795.4337952pt/>^{58}\textrm{NiH}<img src="/tex/24ee684c2922b0d32c54a34089c92ec0.svg?invert_in_darkmode&sanitize=true" align=middle width=4.5662248499999905pt height=14.15524440000002pt/>^{60}\textrm{NiH}<img src="/tex/fd92a53167b3c6ae9574071613d555dc.svg?invert_in_darkmode&sanitize=true" align=middle width=27.11199479999999pt height=22.831056599999986pt/>^{62}\textrm{NiH}<img src="/tex/24405ce93a16e0f497025fd5e5f7bfab.svg?invert_in_darkmode&sanitize=true" align=middle width=700.5030021pt height=282.55708469999996pt/>\mu = m_{A}*m_{B} / (m_{A} + m_{B})<img src="/tex/baca17a163a9cb6b60d7da054e52fec2.svg?invert_in_darkmode&sanitize=true" align=middle width=947.5261718999999pt height=2255.1598245pt/>R_i<img src="/tex/01a0b2843109824344cc6649405eb01d.svg?invert_in_darkmode&sanitize=true" align=middle width=138.1516059pt height=22.831056599999986pt/>R_{min}<img src="/tex/5b7ce0b11756ea13938668594c96a288.svg?invert_in_darkmode&sanitize=true" align=middle width=13.90414904999999pt height=20.221802699999984pt/>R_{max}<img src="/tex/50b0e30a0c2655fedf357501ff3566da.svg?invert_in_darkmode&sanitize=true" align=middle width=132.89916914999998pt height=22.831056599999986pt/><img src="/tex/f2473b861a7f7efddac0503d81572ae6.svg?invert_in_darkmode&sanitize=true" align=middle width=292.75927724999997pt height=24.65753399999998pt/><img src="/tex/d58337d48610f0bf4202bfdfca7f94fe.svg?invert_in_darkmode&sanitize=true" align=middle width=30.182742149999992pt height=39.45205439999997pt/>N_{R}<img src="/tex/83b907c6487854dce8c00196a3cbc702.svg?invert_in_darkmode&sanitize=true" align=middle width=213.55245284999998pt height=22.831056599999986pt/>\Delta_{R}<img src="/tex/709791a6749a4974e2d16b520812ebdb.svg?invert_in_darkmode&sanitize=true" align=middle width=96.51884055pt height=22.831056599999986pt/><img src="/tex/216bed4bde630a321c7ab38b6e7183b6.svg?invert_in_darkmode&sanitize=true" align=middle width=124.77998774999999pt height=29.205422400000014pt/><img src="/tex/b051c4d167aa8c077ee4e155947b83f7.svg?invert_in_darkmode&sanitize=true" align=middle width=604.6591122pt height=124.74886379999998pt/>R_{min}$ and the final $R_{max}$ values of the R grid.
    - It accepts a tuple/list containing two integer or float numbers

The above two arguments are mandatory.

- **```solver```** - the method used for solution of the Schrodinger equation
    - The initialization of the grid depends on which solver method is selected.
    - It should be a string and the following values are available:
      - 'sinc'
      - 'fourier'
      - 'fd5'

      If one of the values 'sinc' or 'fourier' is selected then FGH (Fourier Grid Hamiltonian) method will be used for the solution of the Schrodinger equation.

      If 'fd5' option is selected the Schrodinger equation will be solved by Finite Diffrence method using...

      These two methods are the most frequently applied when the Schrodinger equation is solving for coupled channels problems in diatomic molecular spectroscopy.
    - If this argument is not set the 'sinc' method will be used by default.
    - It is case-insensitive (i.e. 'fd5', 'FD5', 'fD5' are acceptable values).

> **_NOTE:_** All parameters in program which values should be of type string are case-insensitive.

This example shows how the grid object should be initialized for uniform grid.
```python
# create a uniform grid of points
grid = diatom.Grid(npoints=100, rgrid=(1.0, 3.0), solver='sinc')
```
<!-- omit in toc -->
### Nonuniform grid
For the generation of the nonuniform grid an analytical mapping procedure is appled.

If nonuniform grid is requested then both **```alpha```** and **```rbar```** parameters should be set.

- **```alpha```** - 
  - It should be a positive integer or float number
  - If this parameter has a negative value then uniform grid without mapping will be generated. This allows one to easyly switch between uniform and nonuniform grid reperesentations.

- **```rbar```** - 

This example shows how the grid object should be initialized for nonuniform grid.

```python
# create a nonuniform grid of points
grid = diatom.Grid(npoints=100, rgrid=(1.0, 3.0), solver='sinc', alpha=3.0, rbar=1.458)
```
----

Python does not have so restrictive means for contolling the privite/public access of the internally defined data as the other OO (Object-Oriented) languages. In other words we are allowed to access directly some of the internal properties of the declared classes and objects. For example if we would like to see how the generated grid of points looks we can simply write

```python
print(grid.rgrid)
```
That will produce an output like

```
[1.41729459 1.48936042 1.56142625 1.63349207 1.7055579  1.77762373
 1.84968955 1.92175538 1.99382121 2.06588703 2.13795286 2.21001869
 2.28208451 2.35415034 2.42621617 2.498282   2.57034782 2.64241365
 2.71447948 2.7865453  2.85861113 2.93067696 3.00274278 3.07480861
 3.14687444 3.21894026 3.29100609 3.36307192 3.43513774 3.50720357
 3.5792694  3.65133522 3.72340105 3.79546688 3.8675327  3.93959853
 4.01166436 4.08373018 4.15579601 4.22786184 4.29992766 4.37199349
 4.44405932 4.51612515 4.58819097 4.6602568  4.73232263 4.80438845
 4.87645428 4.94852011 5.02058593 5.09265176 5.16471759 5.23678341
 5.30884924 5.38091507 5.45298089 5.52504672 5.59711255 5.66917837]
```
Or if we want to save the generated grid in file we can write

```python
import numpy as np

np.savetxt('grid_points.dat', grid.rgrid)
```
This can facilitate the debugging process.

Although the program has additional external options for storing such information which will be considered in the following sections, this demonstrates how powerful and convenient the OOP approach is.

----

## Channel Definition

A channel is defined as an electronic state with definite values of <img src="/tex/b23332f99af850a48831f80dbf681ed6.svg?invert_in_darkmode&sanitize=true" align=middle width=11.41554479999999pt height=22.465723500000017pt/>, <img src="/tex/e257acd1ccbe7fcb654708f1a866bfe9.svg?invert_in_darkmode&sanitize=true" align=middle width=11.027402099999989pt height=22.465723500000017pt/>, <img src="/tex/813cd865c037c89fcdc609b25c465a05.svg?invert_in_darkmode&sanitize=true" align=middle width=11.87217899999999pt height=22.465723500000017pt/> and <img src="/tex/9432d83304c1eb0dcb05f092d30a767f.svg?invert_in_darkmode&sanitize=true" align=middle width=11.87217899999999pt height=22.465723500000017pt/> quantum numbers.

Each channel have to be defined as an object of type Channel. The parameters that we need to pass to the constructor of the Channel object are:

- **```filep```** - the name of the file containing the potential parameters
  - It is a string parameter referring to an existing file or file path
  - The structure of the file depends on which model function is selected by the propery **```model```** described below.

- **```model```** - defines the potential energy function model (PECs model)
  - It could be a pointwise, analytical (one of the internally defined potentials) or custom analytical potential (defined by the user). The possible values are:
    - _'pointwise'_ : pointwise potential with cubic spline interpolation using internal Python implementation
    - _'cspline'_ : pointwise potential with cubic spline interpolation using our own implementation
    - _'Morse'_ : analytical Morse potential
    - _'EMO'_ : analytical EMO (Expanded Morse Oscilator) potential
    - _'MLR'_ : analytical MLR (Morse/Long-Range) potential
    - _'custom'_ : custom analytical potential

  - Its should be a string and the value is case-insensitive

----

#### The structure of the potential files

- for pointwise potentials i.e. **```model```**=_'pointwise'_ or **```model```** =_'cspline'_ the potential file look like

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

It contains 2 or 3 columns. The first column
The third column is required when running the fitting procedure.
As first row on the top of the file the number of points should be specified. (TODO: remove this).

If some point is no longer needed we can comment it out like:
```python
    1.3858695652       2976.97391482846    0
  # 1.5081521739       1953.39920190690    0
    1.6304347826       2525.33711155168    0
    1.7527173913       4070.10001126557    0
  # 1.8505434783       5627.41358380562    0
```

- for Morse potential i.e. **```model```**=_'Morse'_

```
Te = 1.00000000e-05  0
De = 1.00000000e+01  0
a  = 1.00000000      0
re = 1.27463182      0
```

- for EMO potential i.e. **```model```**=_'EMO'_

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

- for MLR potential i.e. **```model```**=_'MLR'_

...
- for custom potential i.e. **```model```**=_'custom'_

```
param1 = 4.5e-11   0
param2 = 4.5e-05   0
param3 = 9.9e-01   0
param4 = 2.4e+01   0
param5 = 8.4e+04   0
param6 = 5.1e+00   0
```

----

- **```nlambda```** - the quantum number <img src="/tex/b23332f99af850a48831f80dbf681ed6.svg?invert_in_darkmode&sanitize=true" align=middle width=11.41554479999999pt height=22.465723500000017pt/> of the state
  - It is a positive integer number
- **```sigma```** - the quantum number <img src="/tex/813cd865c037c89fcdc609b25c465a05.svg?invert_in_darkmode&sanitize=true" align=middle width=11.87217899999999pt height=22.465723500000017pt/> of the state
  - It is an integer or float number
- **```multiplicity```** - the multiplicity defined as <img src="/tex/eb1977c2439009c2f58e150e1d2769c0.svg?invert_in_darkmode&sanitize=true" align=middle width=47.556992999999984pt height=22.465723500000017pt/> for the state
  - It is a positive integer or float number
- **```rot_correction```** - correction to the diagonal rotational Hamiltonian
  - Optional parameter of type integer or float

```python
ch1 = diatom.Channel(
    filep='poten_morse.pot',
    model='morse',
    nlambda=0,
    sigma=0,
    multiplicity=1,
    rot_correction=0.
)
```

For readability the parameters are written on separate lines.

- **```custom_function```** - the name of the custom analytical function defined by the user

  - Its value should be of type function - a name of already defined function in the current or in any other *.py file.
  - The user defined function should accept exactly 2 input parameter and should return a single parameter. 
  - The first input argument is an array containing the values of the potential parameters as defined in the potential file through the parameter **```filep```**. In this file the parameters have to be defined with the keyword 'param' followed by the number of corresponding parameter. An example is shown below. The second input argument is an array containing the grid points.
  - The output parameter has to be an array containing the calculated values of the potential function on the grid points. The length of the returned array should be equal to the number of grid points.

  - All input and output arguments are 1D numpy arrays of type 'numpy.ndarray'.
  - **```custom_function```** is an optional parameter
  - In this case the parameter **```model```** has to be set to 'custom'.

A simple example of a user defined function implementing the Morse potential is shown below:

```python
import numpy as np

def VMorse(params, rgrid):

    # unpacking the input parameters
    Te, De, a, re = params

    # returns the values of the Morse function over the grid
    return Te + De * np.power((1.0 - np.exp(-a*(rgrid-re))), 2.0)
```
In this case the parameters in the potential file are defined in the following way:
TODO: change this to values of real molecule

```
param1 = 4.55633525e-11   0
param2 = 4.55633525e-05   0
param3 = 9.99989591e-01   0
param4 = 2.40869887e+00   0
```
Note the special construction 'param+number' with conscuitive numbers starting from 1. In the input array the parameters will be aranged by the that number therefore the order in which they are defined in the file does not matter. 

> **_NOTE:_** Since we do the internal calculations in au units...!

After defining all channel objects we need to create a list containing the Channel obejcts which will be included the following compuations. In order for the program to know which are the channels and what the values of their parameters we need to provide that list to it. This is done by calling the class method **```set_channel_parameters()```** and passing the created list with channels as input argument to it:

```python
# combine the defined channels ch1 and ch2 in one list
channels = [ch1, ch2]

# then pass it to this function
diatom.Channel.set_channel_parameters(channels)
```

## Coupling Definition

The parameters that we need to pass to the constructor of the Coupling object:

- **```interact```** - the numbers of the channels which are connected by this interaction
  - should be of type tuple or tuple of tuples
  - the numbers should correspond to a defined channel

- **```coupling```** - the type of the interation
  - the available values are:
    - _'spin-orbit'_ : the diagonal and off-diagonal Spin-Orbit (SO) interaction
    - _'LJ'_ : the rotational L-uncoupling interaction
    - '_SJ'_ : the rotational Spin-uncoupling interaction
    - _'SL'_ : the rotational Spin-electronic interaction
    - _'spin-rot'_ : Spin-Rotational interaction
    - _'spin-spin'_ : Spin-Spin interaction
    - _'LambdaDe'_ : second-order $\Omega$- or $\Lambda-$doubling effect on the e-parity levels
    - _'LambdaDf'_ : second-order $\Omega$- or $\Lambda-$doubling effect on the f-parity levels
  - should be of type string or tuple of strings

- **```model```** - the type of the function which will be used for modelling the interacion
  - should be of type string

- **```multiplier```** - a number which will multiply the defined function
  - integer or float number
  - if not provided the default value will be set to 1.0

- **```label```** - 
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

#### Structure of the couplings file

```
cp1:
- 0.75000000  -300.731917268926   0
- 1.25000000  -300.188348980659   0
- 1.50000000  -295.722686487558   0
- 2.00000000  -305.071069853626   0
- 3.00000000  -301.500000000000   0
cp2:
- 0.75000000  -606.439130285529   0
- 1.25000000  -594.031805306883   0
- 1.50000000  -599.307655859103   0
- 2.00000000  -607.563537893363   0
- 3.00000000  -603.000000000000   0
3:
- 0.75000000     1.024663995832   0
- 1.20000000     1.223710979140   0
- 1.50000000     1.003797203439   0
- 2.00000000     0.995600549099   0
- 3.00000000     1.000000000000   0
```

----

More complicated example for coupling definition:

```python
cp2 = diatom.Coupling(
    interact=((2,4), (3,5), (3,4)),
    coupling=('LJ', 'LJ', 'SL'),
    model='pointwise',
    multiplier=(2.0, 2.0, 2.0),
    label='3'
)
```

## Experimental Data

## Molecule Levels Computation
<!-- omit in toc -->
### Define The Molecule Levels Object

To compute the eigenvalues first we need to initilize an object of type **```MoleculeLevels```**.
Three required parameters have to be provided as in the example:

```python
mlevels = diatom.MoleculeLevels(mdata, grid, channels)
```

The first parameter is the created [Molecule Data](#Molecule-Data-Definition) object, the second one is the created [Grid](#grid-definition) object and the third one is the created list with [Channel](#channel-definition) objects.


If couplings between the channels and experimental data are available then two more optional parameters have to be provided as in the example:

```python
mlevels = diatom.MoleculeLevels(mdata, grid, channels, couplings=couplings, data=exp_data)
```
The first optional parameter is created list of [Coupling](#coupling-definition) objects. The second optional parameter is the returned list of [Experimental Data](#experimental-data).

<!-- omit in toc -->
### Calculate the eigenvalues and eigenvectors

- **```eig_decomp```** - defines the Python <a href="https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.eigh.html" target="_blank">SciPy</a> package to be used for eigenvalues decomposition
  - The two possible values are:

    - _'lapack'_ : will compute the eigenvalues and eigenvectors by calling 
    <a href="https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.eigh.html" target="_blank">scipy.linalg.eigh()</a> procedure

    - _'arpack'_ : will compue the eigenvalues and eigenvectors by calling 
    <a href="https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.eigsh.html#scipy.sparse.linalg.eigsh" target="_blank">scipy.sparse.linalg.eigsh() </a> procedure

  - Default is set to 'lapack'
  - The two SciPy procedures provide high-level interface to standard LAPACK and ARPACK routines written in Fortran and C.

<!-- omit in toc -->
#### Eigenvalue Decomposition with LAPACK

This is the recommended choice when FGH method is selected as solving method since the generated Hamiltonain is a dense matrix except for the case of high number of channels.

The range of eigenvalues and eigenvectors which to be computed can be set by one of the two optional parameters 

- **```energy_subset_index```** - defines a subset of indicies for the returned eigenvalues

  -  the first and the last index of the desired eigenvalues and the corresponding eigenvectors in the whole array of computed eigenvalues
  
  - iterable with 2 integer numbers defining the two indicies

    ```python
    mlevels.calculateLevels(energy_subset_index=(0, 6))
    ```

- **```energy_subset_value```** - defines a subset of values for the returned eigenvalues

  - the smallest and the largest value of the selected range of eigenvalues - only the eigenvalues with values between these two numbers and the corresponding eigenvectors will be returned 

  - iterable with 2 integer or float numbers - the smallest and the largest value

    ```python
    mlevels.calculateLevels(energy_subset_value=(1000., 3000.))
    ```

  > **_NOTE:_**  If neither of them is set the whole range of comupted eigenvalues will be returned

- **```lapack_driver```** - the lapack procedure which will be called for the eigenvalues and eigenvectors computation
  - The possible values are:

    - _'ev'_  : call syev() procedure
    - _'evd'_ : call syevd() procedure
    - _'evr'_ : call syevr() procedure
    - _'evx'_ : call syevx() procedure

  - The default value is set to 'evr'
  > **_NOTE:_**  **```subset_by_index```** and **```subset_by_value```** cannot be used together with 'ev' and 'evd'.

    ```python
    mlevels.calculateLevels(eig_decomp='lapack', lap_driver='evx')
    ```

<!-- omit in toc -->
#### Eigenvalue Decomposition with ARPACK
ARPACK procedure is efficient and suitable for finding the _largest_ eigenvalues of a sparse matrix especially when a few of them are needed. If the _smallest_ eigenvalues need to be computed then it is recommended to use a shift-invert mode. In this mode the original eigenvalue problem will be transformed to an eqivalent problem so that the original small eigenvalues u will correspond to the transformed large eigenvalues v: v = 1 / u

For further information: https://docs.scipy.org/doc/scipy/reference/tutorial/arpack.html

The releated parameters are:
- **```arpack_k```** - the dsired number of eigenvalues and eigenvectors
  - positive integer number
  - should be smaller than the dimension of the matrix 
- **```arpack_which```** - which eigenvalues and eigenvectors to find
  - The possible values are:
    - _'LM'_ : Largest (in magnitude) eigenvalues.
    - _'SM'_ : Smallest (in magnitude) eigenvalues.
    - _'LA'_ : Largest (algebraic) eigenvalues.
    - _'SA'_ : Smallest (algebraic) eigenvalues.
    - _'BE'_ : Half (k/2) from each end of the spectrum.
- **```arpack_sigma```** - if this parameter is specified then shift-invert mode is applied. Then the procedure will return k (specified by **```arpack_k```**) shifted eigenvalues which have values around the specified value of this parameter.
  
  - efficient when the smallest eigenvalues are required (**```arpack_which```**='SM')
  - integer or float number

```python
mlevels.calculateLevels(eig_decomp='arpack', arpack_k=25, arpack_which='SM', arpack_sigma=0.1)
```
> **_NOTE:_** With ARPACK it is not possible to compute all eigenvectors of a matrix.

> **_NOTE:_** The parameters **```lapack_driver```**, **```subset_by_index```** and **```subset_by_value```** are not relevent when **```eig_decomp```**='arpack'.
Similarly **```arpack_k```**, **```arpack_which```** and  **```arpack_sigma```** are not relevent when **```eig_decomp```**='lapack'

<!-- omit in toc -->
### Identification of computed levels
- **```identify```** - defines the method for identification of the calculated energies i.e. the procedure used for the assignment of vibrational quantum number for each level and the determination of the state which it belongs to based on the available experimental data.
  - integer number with 2 possible values:
    - _0_ : identification by energy
    - _1_ : identification by the largest value of the coupling coefficient
    <!-- - _2_ : a slow version of _0_; defined and used only for debugging puropse
    - _3_ : a slow version of _1_; defined and used only for debugging puropse -->
  - defult value is set to _0_ 
  - it is not relevent when experimental data are not provided

<!-- omit in toc -->
### Storing and Formatting the Ouputs
- **```store_predicted```** - 
  
- **```store_info```** - 
  
- **```store_evecs```** - 

- **```sort_output```** - 

- **```sort_predicted```** - 

## The minimal working example
An example with the minimal set of nessecary parameters for running a single channel computations with Morse function is presented.

It is done with only a few lines of code.

```python
import diatom

mdata = diatom.MoleculeData()
mdata.masses = [0.972222222]
mdata.nisotopes = [1]
mdata.jrange = (1, 5)

grid = diatom.Grid(npoints=100, rgrid=(1.0, 3.0))

ch1 = diatom.Channel(
    filep='morse_HCl.pot',
    model='morse',
    nlambda=0,
    sigma=0,
    multiplicity=1,
)

channels = [ch1]
diatom.Channel.set_channel_parameters(channels)

mlevels = diatom.MoleculeLevels(mdata, grid, channels)
mlevels.calculateLevels()
```
The content of the potential file:
```
Te =    1.00000000e-05   0
De =    1.00000000e+01   0
a  =    1.00000000       0
re =    1.27463182       0
```
Here all eigenvalues and eigenvectors for J=1,2,3,4 and 5 with both e/f parity labels will be computed with the 'sinc' option of FGH method. The computed eigenvalues will be stored in file called 'eigenvalues_all.dat'. The first few lines of this file look like

```
#     No   v      Ecalc          J  parity  marker      CC1   state  lambda   omega
      1   0      56.113986     1.0     0      1     1.000      1      0     0.00
      2   1     180.940849     1.0     0      1     1.000      1      0     0.00
      3   2     387.925302     1.0     0      1     1.000      1      0     0.00
      4   3     677.259248     1.0     0      1     1.000      1      0     0.00
      5   4    1049.088043     1.0     0      1     1.000      1      0     0.00
      6   5    1503.471257     1.0     0      1     1.000      1      0     0.00
      7   6    2040.431414     1.0     0      1     1.000      1      0     0.00
      8   7    2659.982664     1.0     0      1     1.000      1      0     0.00
      9   8    3362.127325     1.0     0      1     1.000      1      0     0.00
     10   9    4146.873422     1.0     0      1     1.000      1      0     0.00
```
# Fitting of Calculated Energy Levels

# Plotting

![image info](./hcolormesh.png)

# Computing the Transition Frequencies

# Computing the Transition Intensities

## Einstein A coefficient and Radiative Lifetime

## Line strength

<!-- ## Cross-section -->
