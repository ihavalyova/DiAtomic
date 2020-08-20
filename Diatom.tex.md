
- [**Diatomic** module: what is it used for?](#diatomic-module-what-is-it-used-for)
- [**Diatom** module: how to install and setup?](#diatom-module-how-to-install-and-setup)
- [Diatomic molecule: basic theoretical concepts](#diatomic-molecule-basic-theoretical-concepts)
  - [The total Hamiltonian and basis functions](#the-total-hamiltonian-and-basis-functions)
  - [The Scrodinger equation for a single state and coupled system of states](#the-scrodinger-equation-for-a-single-state-and-coupled-system-of-states)
  - [The interaction terms and their matrix elements](#the-interaction-terms-and-their-matrix-elements)
  - [Methods for solving the Schrodinger equation](#methods-for-solving-the-schrodinger-equation)
  - [Potential Energy function models (PECs models)](#potential-energy-function-models-pecs-models)
- [Computing Energy Eigenvalues](#computing-energy-eigenvalues)
  - [Molecule Data Object Definition](#molecule-data-object-definition)
  - [Grid Object Definition](#grid-object-definition)
  - [Channel Object Definition](#channel-object-definition)
  - [Coupling Object Definition](#coupling-object-definition)
      - [Shared parameters](#shared-parameters)
  - [Experimental Data](#experimental-data)
  - [Molecule Levels Computation](#molecule-levels-computation)
  - [Examples](#examples)
- [Fitting of the Calculated Energy Levels](#fitting-of-the-calculated-energy-levels)
  - [SVD Fit](#svd-fit)
  - [Minuit Fit](#minuit-fit)
  - [Levenberg-Marquard Fit](#levenberg-marquard-fit)
- [Computing the Transition Frequencies](#computing-the-transition-frequencies)
  - [States represented by potentials](#states-represented-by-potentials)
  - [States represented by term values](#states-represented-by-term-values)
- [Computing the Transition Intensities](#computing-the-transition-intensities)
  - [Line strength](#line-strength)
    - [Honl-London Factors](#honl-london-factors)
    - [Frank-Condon Factors](#frank-condon-factors)
  - [Einstein A coefficient and Radiative Lifetime](#einstein-a-coefficient-and-radiative-lifetime)
- [Plotting](#plotting)


# **Diatomic** module: what is it used for?
The Python package **```Diatomic```** allows various calculations for diatomic molecules to be performed. It supports single and coupled channels computations of bound rovibrational levels, intensity calculations, fitting to the experimental data.
The current functionality covered by the program includes:
* ..
* ..

Just as an example of what you can do with **```Diatomic```** module, if you have a set of potentials for a couple of molecular electronic states represented by points (abinitio, RKR and etc.) and want to see how they look you can do something like:
```python
p = diatom.Plotting()
p.plot_potentials_points(['p1.pot', 'p2.pot', 'p3.pot', 'p4.pot'], show=True, ipoints=50, ylim=(1e3, 1e5))
```
or even simpler:
```python
import glob
p = diatom.Plotting()
p.plot_potentials_points(glob.glob('./*.pot'), show=True, ipoints=50, ylim=(1e3, 1e5))
```
assuming your potential files are in the current directory.

# **Diatom** module: how to install and setup?

**```Diatom```** module can be installed from the Python software repository PyPI (Python Package Index) via pip. From Linux command line execute

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

import diatom
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

The **```Diatom```** package is extensivly tested on Linux platform but works under Windows and MacOS as well.

# Diatomic molecule: basic theoretical concepts

## The total Hamiltonian and basis functions

The total Hamiltonian of a diatomic molecule in the rotating molecular frame with origin at the center of mass of the molecule can be written as a sum of several terms:

$$
\mathbf{H} = \mathbf{T}_{\mathrm{N}}(R) + \mathbf{H}_{\mathrm{rot}}(R, \theta, \phi) + \mathbf{T}_{\mathrm{e}}(r) + \mathbf{V}(R, r) + \mathbf{H}_{\mathrm{rel}}
$$

where the each of the terms is defined as follows:

## The Scrodinger equation for a single state and coupled system of states
...
<!-- omit in toc -->
### Single channel approximation for an isolated state

The energy eigenvalues of single isolated state of diatomic molecule can be obtained by solving the radial Schrodinger equation

$$
\left( \frac{-\hbar^{2}}{2\mu} \frac{d^{2}}{dR^{2}} + U(R) + \frac{\hbar^{2}}{2\mu R^{2}}J(J+1) \right) \phi_{vJ}(R) = E_{vJ} \phi_{vJ}(R)
$$

with internuclear distance labeled with $R$, the sum of the second and the third term $U(R) + ({\hbar^2}/{2\mu R^2})J(J+1)$ from the left hand side is called effective potential energy curve, the reduced mass is $\mu = M_{1} M_{2} /(M_{1} + M_{2})$ where $M_1$ and $M_2$ are the masses of the two atoms, J is the rotational quantum number; $E_{vJ}$ are the energy eigenvalues of the rovibrational levels and $\phi_{vJ}$ are their corresponding eigenfunctions.

<!-- omit in toc -->
### The coupled channels problem

## The interaction terms and their matrix elements

## Methods for solving the Schrodinger equation
<!-- omit in toc -->
  ### Fourier Grid Hamiltonian (FGH)

<!-- omit in toc -->
  ### Finite Difference Method (FD)

## Potential Energy function models (PECs models)
- Pointwise potential

- Analytical potential
   - Morse potential

  The Morse potential function is defined as:

  $$
  V(R) = T_{e} +D_{e}[1 - e^{\beta(r-r_{e})}]^2
  $$
  where Te is the term value energy??, $D_{e}$ measures the energy from the bottom of the potential to the dissociation limit (it different from the dissociation energy)
  $\beta$ is a constant and $r_{e}$ is the equilibrium internuclear distance.

   - EMO (Expanded Morse Oscillator) potential

  It is defined as

  $$ 
  V_{EMO}(R) = T_{e} + D_{e}[1 - exp(-\beta_{EMO}(R).(R-R_{e}))]^{2} 
  $$

  which is a Morse potential function with a radial dependence of the exponential coefficient

  $$
  \beta_{EMO}(R) = \sum_{i=0}^{N} \beta_{i} . y(R)^{i}
  $$

  represented as a polynomial expansion in the powers of the dimensionless function

  $$
  y(R) = \frac{R^{p} - R_{e}^{p}}{R^{p} + R_{e}^{p}}
  $$

  with p being a positive integer number.

   - MLR (Morse/Long-Range) potential

$$
V_{MLR}(R) = T_{e}
$$

# Computing Energy Eigenvalues

## Molecule Data Object Definition

Initially an object of type **```MoleculeData```** should be instantiated in the already created *.py file which we called main.py:

```python
# creating molecule data object called mdata
mdata = diatom.MoleculeData()
```

Part of the input information about the molecule should be defined by the following properties and attached to the created **```MoleculeData```** object:

<!-- omit in toc -->
### Reduced masses and isotopes

- **```molecule```** - defines one or more isotopic forms of the same molecule by specifing their symbols
  - it should be an iterable of type list or tuple of strings.
  - each defined item inside this iterable represents a symbol corresponding to diffrent isotope of the same molecule. 
  - the molecule symbols should be specified by thier atomic symbols and each mass number put in front of them; only spaces between the symbols are allowed.
  - The reduced masses for each isotope will be computed by looking for their atomic masses in an existing database (ascii file) available from
  https://www.nist.gov/pml/atomic-weights-and-isotopic-compositions-relative-atomic-masses
  - this property does not determine which and how many isotopes will be included in the calculations (refer to the property **```nisotopes```** below) but only defines the masses by the symbols.
  - it is not mandatory

In the following example the symbols for three of the isotopes of NiH molecule - $^{58}\textrm{NiH}$, $^{60}\textrm{NiH}$ and $^{62}\textrm{NiH}$ are defined as: 

```python
# define the symbols for three isotopes
mdata.molecule = ['58Ni1H', '60Ni1H', '62Ni1H']

# this is also a valid syntax
mdata.molecule = ['58Ni 1H']
mdata.molecule = ['58 Ni 1 H']
```

- **```masses```** - defines the reduced masses for one or more isotopes of the molecule by specifing their numerical values comupted beforehand. 

  - use when the automatic computation of the reduced masses is not preffered
  - it should be an iterable of type list or tuple of float numbers.
  - each defined item inside this iterable represents the reduced mass for a different isotope computed by $\mu = m_{A}*m_{B} / (m_{A} + m_{B})$ in amu units.
  - this property does not determine which and how many isotopes will be included in the calculations (refer to the property **```nisotopes```** below) but only defines their masses.
  - it is not mandatory

```python
# define the reduced masses for three of NiH isotopes
mdata.masses = [0.990592988928, 0.991157254368, 0.991686280089]

# or from the masses of the separate atoms
mA, mB = 1.00782503223, 34.968852682
mdata.masses = [mA * mB / (mA + mB)]
```

> **_NOTE:_** At least one of the above two properties (**```molecule```** or **```masses```**) should be specified (i.e. one of them is mandatory) otherwise the execution fails.

> **_NOTE:_** The properties **```molecule```**/**```masses```** and **```nisotopes```** are closely connected.

- **```nisotopes```** - tells which isotopes to be included in the compuations by providing their numbers
  - it should be a list or tuple of integer numbers
  - a list of numbers each corresponding to an isotope with mass defined by the property **```masses```**. In other words the values of this property are indicies of the property **```masses```** counting from 1.
  - the values starts from 1 up to the number of masses in **```masses```**
  - it is mandatory

The energies only for the second and the third isotope will be computed in this example:
```python
# define the masses of 3 isotopes
mdata.molecule = ['58Ni1H', '60Ni1H', '62Ni1H']

# but compute only for the second and the third one
mdata.nisotopes = [2,3]
```

<!-- omit in toc -->
### Values of the rotational quantum number J

- **```jrange```** - defines a continious sequence of J values by specifing the inital and the final values.

  - it is an iterable of type list or tuple containing only 2 elements (integer or float numbers)
  - the two elements are the initial and the final J value of the required range of rotational quantum numbers which to be used in the computations

```python
# the J values will be: 0, 1, 2, 3, 4 and 5
mdata.jrange = (0, 5)
```

- **```jvalues```** - defines a specific set of values for J

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

- **```parities```** - tells which parity levels to be computed. 
  - levels with one or both e/f parity labels can be computed. 
  - we use a convention for designating the parity labels with numbers 0 and 1 for f and e respectivly.
  - it accepts a sinle value or tuple/list of one or two values of integer type (0 or 1).
  - not mandatory - if it is not set both e- and f-levels will be computed by default.

```python
# both e- and f-levels will be computed
mdata.parities = 0, 1

# only e-levels will be computed
mdata.parities = 1
```

<!-- omit in toc -->
### Reference level

- **```referenceJ```** - a reference J level which energy will be used for shifting all other levels
  - the energies of all levels will be shifted by the energy of the level having this J value
  - in the general case the e/f levels and levels corresponding to diffrent isotopes might be shifted with diffrent energy
  - integer or float number

- **```referenceE```** - a reference energy value which will be used for shifting all levels
  - all levels will be shifted with the same value
  - integer or float number

Examples:
```python
mdata.referenceJ = 2.5
```

```python
mdata.referenceE = 1000
```
None of them is mandatory.

If both of these parameters are specified simultaneously, **```referenceJ```** will be used.

## Grid Object Definition

Will define a grid of points on the internuclear distance R.

To set all parameters releated to the defininiton of the grid we need to initilize the grid object. This can be done if we call the **```Grid```** constructor and pass the required arguments. Two types of grids are supported:
<!-- omit in toc -->
### Uniform grid

An equdistant set of points will be generated on which all subsequent calculations will be performed. The grid points $R_i$ in the interval from $R_{min}$ to $R_{max}$ are determined by:

$$
R_i = R_{min} + (i-1)\Delta_{R}, \qquad i = 1\dots N_{R}
$$

where $N_{R}$ is the number of grid points and $\Delta_{R}$ is the grid step

$$
\Delta_{R} = \frac{R_{max} - R_{min}}{N_{R} - 1}
$$

The required arguments in this case are:

- **```npoints```** - the number of grid points
    - It should be a positive integer number

- **```rgrid```** - the range of R values.
    - We need to specify the initial $R_{min}$ and the final $R_{max}$ values of the R grid.
    - It accepts a tuple/list containing two integer or float numbers

The above two arguments are mandatory.

- **```solver```** - the method used for solution of the Schrodinger equation
    - The initialization of the grid depends on which solver method is selected.
    - It should be a string and the following values are available:
      - **```'sinc'```**
      - **```'fourier'```**
      - **```'fd5'```**

      If one of the values **```'sinc'```** or **```'fourier'```** is selected then FGH (Fourier Grid Hamiltonian) method will be used for the solution of the Schrodinger equation.

      If 'fd5' option is selected the Schrodinger equation will be solved by Finite Diffrence method using...

      These two methods are the most frequently applied when the Schrodinger equation is solving for coupled channels problems in diatomic molecular spectroscopy.
    - **```solver='sinc'```** is set by default.

> **_NOTE:_** All parameters in program with values of type string are case-insensitive.

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
<!-- omit in toc -->
### Comment
Python does not have so restrictive means for contolling the privite/public access of the internally defined data as the other OO (Object-Oriented) languages. In other words we are allowed to access directly some of the internal properties of the declared classes and objects. For example if we would like to see how the generated grid of points looks we can simply write

```python
print(grid.rgrid)
```
That will produce an output like

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
Or if we want to save the generated grid in file we can write

```python
import numpy as np

np.savetxt('grid_points.dat', grid.rgrid)
```
This can facilitate the debugging process.

Although the program has additional external options for storing such information which will be considered in the following sections, this demonstrates how powerful and convenient the OOP approach is.

----

## Channel Object Definition

A channel is defined as an electronic state with definite values of $\Lambda$, $S$, $\Sigma$ and $\Omega$ quantum numbers.

Each channel have to be defined as an object of type **```Channel```**. The parameters that we need to pass to its constructor are:

- **```filep```** - the name of the file containing the potential parameters
  - It is a string parameter referring to an existing file or file path
  - The structure of the file depends on which model function is selected by the propery **```model```** described below.

- **```model```** - defines the potential energy function model (PECs model)
  - It could be a pointwise function, builtin analytical potential function or custom analytical potential function (defined by the user). The possible values are:
    - **```'pointwise'```** : pointwise potential with cubic spline interpolation using internal Python implementation
    - **```'cspline'```** : pointwise potential with cubic spline interpolation using our own implementation
    - **```'Morse'```** : analytical Morse potential
    - **```'EMO'```** : analytical EMO (Expanded Morse Oscilator) potential
    - **```'MLR'```** : analytical MLR (Morse/Long-Range) potential
    - **```'custom'```** : custom analytical potential

----

<!-- omit in toc -->
#### Structure of the potential files

- for pointwise potentials i.e. **```model='pointwise'```** or **```model='cspline'```** the potential file look like

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

- **```nlambda```** - the quantum number $\Lambda$ of the state
  - It is a positive integer number
- **```sigma```** - the quantum number $\Sigma$ of the state
  - It is an integer or float number
- **```multiplicity```** - the multiplicity defined as $2S+1$ for the state
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

  - All input and output arguments are 1D **```numpy```** arrays of type **```'numpy.ndarray'```**.
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

## Coupling Object Definition

The parameters that we need to provide in order to initilize a **```Coupling```** object are:

- **```interact```** - the numbers of the channels which are connected by this interaction
  - should be of type tuple or tuple of tuples
  - the numbers should correspond to a defined channel

- **```coupling```** - the type of the interation
  - the available values are:
    - **```'spin-orbit'```** : the diagonal and off-diagonal Spin-Orbit (SO) interaction
    - **```'LJ'```** : L-uncoupling interaction
    - **```'SJ'```** : Spin-uncoupling interaction
    - **```'SL'```** : Spin-electronic interaction
    - **```'spin-rot'```** : Spin-Rotation interaction
    - **```'spin-spin'```** : Spin-Spin interaction
    - **```'LambdaDe'```** : second-order $\Omega$- or $\Lambda-$doubling effect on e-parity levels
    - **```'LambdaDf'```** : second-order $\Omega$- or $\Lambda-$doubling effect on f-parity levels
    - **```'LambdaD'```** : second-order $\Omega$- or $\Lambda-$doubling effect on e- and f-parity levels
  - should be of type string or tuple of strings

- **```model```** - the type of the function which will be used for modelling the interacion
  - should be of type string
  - The possible options are:
    - **```pointwise```**
    - **```cspline```**
    - **```custom```**

- **```multiplier```** - integer or float number which will multiply the defined function
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

<!-- omit in toc -->
#### Structure of the couplings file
The couplings file uses a <a href="https://yaml.org/" target="_blank">YAML</a> syntax.

```yaml
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
cp3:
- 0.75000000     1.024663995832   0
- 1.20000000     1.223710979140   0
- 1.50000000     1.003797203439   0
- 2.00000000     0.995600549099   0
- 3.00000000     1.000000000000   0
```

----

#### Shared parameters

When two or more interacting states coupled by the same or different interactions share the same set of parameters a more complicated construction of the **```Coupling```** object is possible.
This is what frequently happens in the case of states coupled by the $L_+$ operator. Here is an example for the interaction ${^2\Pi}$ ~ $^2\Delta$ in the lowest doublet states of NiH molecule:

```python
cp2 = diatom.Coupling(
    interact=((2,4), (3,5), (3,4)),
    coupling=('LJ', 'LJ', 'SL'),
    model='pointwise',
    multiplier=(2.0, 2.0, 2.0),
    label='cp2'
)
```

If we have defined the channels 2, 3, 4 and 5 as ${^2\Pi_{1/2}}$, $^2\Pi_{3/2}$, $^2\Delta_{3/2}$ and $^2\Delta_{5/2}$ then the pairs (2,4), (3,5) and (3,4) are connected by the same $L_+$ operator in two different rotational interactions. Defined in this way they will use the same set of parameters - those labeled by 'cp2'. This type of definition is ....... Of course we could define each of these pairs of interacting states as separate **```Coupling```** objects each having labels and the results will be the same. But then the fit will treat them .....

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

- **```eig_decomp```** - defines whcih <a href="https://docs.scipy.org/" target="_blank">SciPy</a> procedure to be used for eigenvalues decomposition
  - The two possible values are:

    - **```'lapack'```** : calls the **```scipy```** builtin 
    <a href="https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.eigh.html" target="_blank">scipy.linalg.eigh</a> function to compute the eigenvalues and eigenvectors.

    - **```'arpack'```** : calls the **```scipy```** builtin
    <a href="https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.eigsh.html#scipy.sparse.linalg.eigsh" target="_blank">scipy.sparse.linalg.eigsh</a> function to compute the eigenvalues and eigenvectors.

    These two **```scipy```** procedures provide high-level interface to standard LAPACK and ARPACK routines written in Fortran and C.
  - Default is set to **```'lapack'```**


<!-- omit in toc -->
#### Eigenvalue Decomposition with LAPACK

This is the recommended choice when FGH method is selected as a solving method since the generated Hamiltonain matrix except for the case of high number of channels is not a sparse matrix.

If we need the eigenvalues and the eigenvectors only in a certain interval of values or indicies then it is possible to get only a subset of them which will likely reduce the computation time. The range of eigenvalues and eigenvectors which to be computed can be set by one of the two optional parameters:

- **```energy_subset_index```** - defines a subinterval of indicies for the returned eigenvalues

  -  the first and the last index of the desired subinterval of eigenvalues and their eigenvectors in the whole array of computed eigenvalues
  - iterable of type list/tuple with 2 integer numbers defining the two indicies

    ```python
    mlevels.calculate_levels(energy_subset_index=(0, 6))
    ```

- **```energy_subset_value```** - defines a subinterval of values for the returned eigenvalues

  - the smallest and the largest value of the selected range of eigenvalues - only the eigenvalues with values between these two numbers and the corresponding eigenvectors will be returned 

  - iterable if type list/tuple with 2 integer or float numbers - the smallest and the largest value

    ```python
    mlevels.calculate_levels(energy_subset_value=(1000., 3000.))
    ```

  > **_NOTE:_**  If neither of them is set all possible eigenvalues and thier eigenvectors will be computed and returned

- **```lapack_driver```** - the lapack routine which computes the eigenvalues and eigenvectors
  - The possible values are:
    - **```'ev'```**  : calls dsyev routine which computes _all_ eigenvalues and eiegnvectors of a real symmetric matrix
    - **```'evd'```** : calls dsyevd routine which computes _all_ eigenvalues and eiegnvectors for real symmetric matrix using a divide and conquer algorithm
    - **```'evr'```** : calls dsyevr routine which computes a _selected_ eigenvalues and eigenvectors of a real symmetric matrix using the RRR algorithm
    - **```'evx'```** : calls dsyevx routine which computes a _selected_ eigenvalues and eiegnvectors of a real symmetric matrix

  - The default value is set to **```'evr'```** which in general is the recommended choice. **```'evx'```** is faster when only a few of the eigenvalues are desired.

  > **_NOTE:_**  **```subset_by_index```** and **```subset_by_value```** cannot be used together with **```'ev'```** and **```'evd'```** because they compute all eigenvalues.

    ```python
    mlevels.calculate_levels(eig_decomp='lapack', lap_driver='evx')
    ```

<!-- omit in toc -->
#### Eigenvalue Decomposition with ARPACK
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


- **```store_predicted```** - save all computed eigenvalues in file
  
- **```store_info```** - save more information regarding the comutaions in file
  - may be used for debugging purpose
  
- **```store_evecs```** - save the computed eigenvectors in files

- **```sort_output```** - sort the selected eigenvalues before saving them by providing column numbers

- **```sort_predicted```** - sort all computed eigenvalues before saving them by providing column numbers

## Examples

<!-- omit in toc -->
### The minimal working example
An example with the minimal set of nessecary parameters for running a single channel computations with Morse function is presented.

It is done with only a few lines of code.

```python
#!/usr/bin/env python

import diatom

mdata = diatom.MoleculeData()
mdata.molecule = ['58Ni1H']
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
mlevels.calculate_levels()
```
The content of the input potential file:
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
```

<!-- omit in toc -->
### Another example

<!-- omit in toc -->
### A more complex example

# Fitting of the Calculated Energy Levels

The **```Diatomic```** module has several implemented procedures for weighted least-squares fitting of the calculated to the experimental energy levels. In all cases what we want to minimize is the difference between the experimental (observed) energies and the calculated ones. Therefore we define the $\chi^2$ function as:

$$
\chi^2 = \sum_{i=1}^{n} \left[ \frac{E_{i}^{obs} - E_{i}^{cal}(\mathbf{x})}{\sigma_i} \right]^2
$$

where the computed energies are functions of the parameters that we would like to determine by minimizing the $\chi^2$ value. The dependance $E_{i}^{cal}(\mathbf{x})$ is in general nonlinear therefore an iterative procedure will be applied - starting from some trial values of the parameters a corrections will be generated and added to the current values on each iteration that will improve the $\chi^2$ value. The corrections will be found by solving the system:

$$
E_{i}^{obs} - E_{i}^{cal}(\mathbf{x}^{(0)}) - \sum_{j=1}^{m} \frac{\partial{E_{i}^{cal}}}{\partial x_{j}} \bigg\rvert_{x_j=x_{j}^{(0)}} \Delta x_{j} = 0
$$

where we have approximated the dependance $E_{i}^{cal}(\mathbf{x})$ by the first two terms in its Taylor expansion around $\mathbf{x}^{(0)}$. This is a linear system of n equations with m unknowns (usually n > m) in the form $\hat{A}\mathbf{x} - \mathbf{b}$ where the unknown vector x is the vector with the corrections $\Delta x$, the right-hand side vector b is $E_{i}^{obs} - E_{i}^{cal}(\mathbf{x}^{(0)})$, and the coefficient matrix A is formed by the first derivatives of the energies with respect to the parameters. The overall goal of the fit could be summirzied as:

$$
\min_{\mathbf{x}}\;\chi^2 = | \hat{A}\mathbf{x} - \mathbf{b} |^2
$$

As a first step we need to initialize the Fitting object like
```python
fit = diatom.Fitting(mlevels, progress=True)
```

The first parameter is the created **```MoleculeLevels```** object and the second parameter **```progress```** is optional and tells whether a more detailed output to be printed during all the iterations.

## SVD Fit

In general it is not recommnded to solve the above linear system by the method of the normal equations (that uses the matrix inverse) since the matrix is either singular or very close to singular. Sometimes there exist two or more linear combinations of...

Linear Algebra technique for ...

SVD is very special matrix factorization because it can be applied to _any_ matrix, it is unique and guaranteed to exist.

The **```run_svd```** method has only default parameters:

- **```niter```** - the number of iterations. Default is **```niter=1```**

- **```deriv```** - the method used for computation of derivatives. Default is **```deriv='n'```**. When **```deriv='n'```** the derivatives will be computed numerically with the finite diffrence method. Another possibility is **```deriv='a'```** in which case the derivatives of the energies with respect to the parameters will be computed analitically with the Hellman-Feynman theorem:

$$
\frac{\partial E}{\partial a} = \langle \Psi \vert \frac{\partial H}{\partial a}\vert \Psi \rangle
$$

- **```tol```** - the tolerance value. It determines which and how many linear combinations of the fitted parameters will be discarded because the matrix is singular. The singular values which are less than **```tol```** times the largest singular value are treated as zero. The rank of the matrix is determined by the number of the nonzero singular values and **```tol```** governs the effective rank of the matrix that is the number of singular values smaller than some specific value. Default is **```tol=0.1```**

- **```lapack_driver```** - The name of the LAPACK routine used to solve the least-squares problem. The three possible values **```'gelsd'```**, **```'gelsy'```** and **```'gelss'```** correspond to the names of these routines. **```'gelsd'```** uses singular value decomposition and a divide and conquer method, **```'gelsy'```** uses the orthogonal QR factorization of the matrix and **```'gelss'```** uses the singular value decomposition of the matrix. In the most cases **```'gelsd'```** will likely be the most efficient method. Default is **```lapack_driver='gelsd'```**.

- **```step_size```** - used to determine the change in the parameters during the computation of derivatives. Default is **```step_size=1.0e-4```**.

- **```is_weighted```** - whether to apply a weighted least-squares fitting with the method proposed by J. Watson (default is **```False```**) in which case the weights $\sigma_{k}^{-2}$ given by the experimental uncerteinty will be replaced by the expression:

$$
\frac{1}{\sigma_{k}^{2} + 0.3(E_{k}^{exp} - E_{k}^{calc})^{2}}
$$

- **```restart```** - not yet implemented
- **```limit```** - not yet implemented
- **```regular```** - not yet implemented


To find the least-squares solution **```run_svd```** calls <a href="https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.lstsq.html" target="_blank">scipy.linalg.lstsq</a> function. In the most widespread numerical libraries and packages the implementations which use SVD are based on the LAPACK implementation in Fortran - the routine called DGESVD.

## Minuit Fit

It is possible to call the MINUIT

To run the Minuit Fit we should call the method **```run_minuit```** through the created **```Fitting```** object. This method has only optional parameters:

- **```niter```** - The number of iterations. Default is **```niter=0```**.
- **```step_size```** - Default is **```step_size=1.0e-4```**
- **```correlation```** - Default is **```correlation=False```**.
- **```covariance```**.  Default is **```covariance=False```**
- **```uncert```** - Default is **```uncert=False```**

Example

```python
fit.run_minuit(niter=5, step_size=1.0e-3, correlation=True)
```

## Levenberg-Marquard Fit

Not yet implemented

# Computing the Transition Frequencies

## States represented by potentials

## States represented by term values

# Computing the Transition Intensities

## Line strength

### Honl-London Factors

### Frank-Condon Factors

## Einstein A coefficient and Radiative Lifetime

# Plotting

![image info](./hcolormesh.png)

<!-- ## Cross-section -->
