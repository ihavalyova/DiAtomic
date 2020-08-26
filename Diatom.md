
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
  - [Experimental Data](#experimental-data)
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

The **```DiAtom```** module is extensivly tested on Linux platform but works under Windows and MacOS as well.

# Diatomic molecule: basic theoretical concepts

## The total Hamiltonian and basis functions

The total Hamiltonian of a diatomic molecule in the rotating molecule-fixed coordinate system with origin at the center of mass of the molecule can be written as a sum of several terms:

<p align="center"><img src="/tex/0befcf4d6f64295f51648ed6605604bf.svg?invert_in_darkmode&sanitize=true" align=middle width=381.66411525pt height=16.438356pt/></p>

where <img src="/tex/88eff123eb81302eced70d075cfd831d.svg?invert_in_darkmode&sanitize=true" align=middle width=49.06969154999999pt height=24.65753399999998pt/> and <img src="/tex/1565b8417026bdbae4f63e0a7df02e21.svg?invert_in_darkmode&sanitize=true" align=middle width=90.41667525pt height=24.65753399999998pt/> are the vibrational and rotational part of the total nuclear kinetic energy operator in spherical polar coordinates, <img src="/tex/c3a1d7f8e57bd6a18adbaa8585b1b61b.svg?invert_in_darkmode&sanitize=true" align=middle width=40.47569294999999pt height=24.65753399999998pt/> is the kinetic energy of the electrons, 

<!-- <p align="center"><img src="/tex/33361e7842bc9784188a404d33d61dfb.svg?invert_in_darkmode&sanitize=true" align=middle width=107.76820724999999pt height=16.438356pt/></p> -->

## The Scrodinger equation for a single state and coupled system of states

<!-- omit in toc -->
### Single channel approximation for an isolated state

The energy eigenvalues of single isolated state of a diatomic molecule can be obtained by solving the radial Schrodinger equation

<p align="center"><img src="/tex/04ec33ae2ec81c77f6d478424e026fe4.svg?invert_in_darkmode&sanitize=true" align=middle width=428.47916925pt height=40.11819404999999pt/></p>

with internuclear distance labeled with <img src="/tex/1e438235ef9ec72fc51ac5025516017c.svg?invert_in_darkmode&sanitize=true" align=middle width=12.60847334999999pt height=22.465723500000017pt/>, the sum of the second and the third term <img src="/tex/f11fee23e73e65ca5d56841195727865.svg?invert_in_darkmode&sanitize=true" align=middle width=196.35702239999998pt height=26.76175259999998pt/> from the left hand side is called effective potential energy curve, the reduced mass is <img src="/tex/9c774da583d642ab71cf6fc72e1cb939.svg?invert_in_darkmode&sanitize=true" align=middle width=166.2062853pt height=24.65753399999998pt/> where <img src="/tex/6f549764f2f97bec950c14de5352994a.svg?invert_in_darkmode&sanitize=true" align=middle width=22.500061649999992pt height=22.465723500000017pt/> and <img src="/tex/dced8cd0d35e2af2d3499c10d7ee6289.svg?invert_in_darkmode&sanitize=true" align=middle width=22.500061649999992pt height=22.465723500000017pt/> are the masses of the two atoms, J is the rotational quantum number; <img src="/tex/311f29e88a451816c3e621e910343bf4.svg?invert_in_darkmode&sanitize=true" align=middle width=27.481438049999987pt height=22.465723500000017pt/> are the energy eigenvalues of the rovibrational levels and <img src="/tex/6130302efada133fa7ea1cd1ccb12a8f.svg?invert_in_darkmode&sanitize=true" align=middle width=25.14126119999999pt height=22.831056599999986pt/> are their corresponding eigenfunctions.

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

- <img src="/tex/b23332f99af850a48831f80dbf681ed6.svg?invert_in_darkmode&sanitize=true" align=middle width=11.41554479999999pt height=22.465723500000017pt/> and <img src="/tex/9432d83304c1eb0dcb05f092d30a767f.svg?invert_in_darkmode&sanitize=true" align=middle width=11.87217899999999pt height=22.465723500000017pt/> Doubling

## Methods for solution of the Schrodinger equation

Finite-Difference and Fourier Grid Hamiltonain (DVR type method) are the most frequently applied methods for numerical solution of the 1D Schordinger equation for single and coupled channels problems in molecular spectroscopy. In both methods the wavefunction is approximated over an equidistant or non-equidistant grid of points.

<!-- omit in toc -->
#### Uniform grid

In this case the grid points <img src="/tex/6660896e4379722ff79bba94961b201c.svg?invert_in_darkmode&sanitize=true" align=middle width=17.132374049999992pt height=22.465723500000017pt/> in the interval from <img src="/tex/3076bed5ea0f6f7be4dcc194b5375a3b.svg?invert_in_darkmode&sanitize=true" align=middle width=36.92324789999999pt height=22.465723500000017pt/> to <img src="/tex/ea449f9e9a48e2959872aac8fa65e1ca.svg?invert_in_darkmode&sanitize=true" align=middle width=38.73108029999999pt height=22.465723500000017pt/> are determined by:

<p align="center"><img src="/tex/6f277b0df057cb1fe10512647772c30a.svg?invert_in_darkmode&sanitize=true" align=middle width=292.75927724999997pt height=16.438356pt/></p>

where <img src="/tex/83b80e12323c9ced450508d7bd3822a5.svg?invert_in_darkmode&sanitize=true" align=middle width=23.169583799999987pt height=22.465723500000017pt/> is the number of grid points and <img src="/tex/5385180385f31e03adfa5256148b06bb.svg?invert_in_darkmode&sanitize=true" align=middle width=23.66048849999999pt height=22.465723500000017pt/> is the grid step

<p align="center"><img src="/tex/f17deb669038ee0be520dfb23a4a35c0.svg?invert_in_darkmode&sanitize=true" align=middle width=145.76189265pt height=36.09514755pt/></p>

<!-- omit in toc -->
#### Nonuniform grid



<!-- omit in toc -->
  ### Finite-Difference Method (FD)

The second derivative of the wavefunction with respect to the internuclear distance is approximated by five-point central difference schema:

<p align="center"><img src="/tex/ad06c975f5915bd636d67109e3ca03c9.svg?invert_in_darkmode&sanitize=true" align=middle width=465.42534059999997pt height=39.452455349999994pt/></p>

The kinetic energy matrix elements are then computed:

<img src="https://render.githubusercontent.com/render/math?math=
T_{ij} = \frac{\hbar^2}{2\mu \Delta R^2} \times
<p align="center"><img src="/tex/0fe33e26f42f052e3b11185f5494aaa3.svg?invert_in_darkmode&sanitize=true" align=middle width=159.59567084999998pt height=21.14534895pt/></p>

<!-- <p align="center"><img src="/tex/cc7921b82598f66df521a9622aeb6bad.svg?invert_in_darkmode&sanitize=true" align=middle width=555.80065365pt height=42.37444365pt/></p>
V_{ij} = V(R_i) \delta_{ij}
<p align="center"><img src="/tex/d7a5280438df67132698f79caff21763.svg?invert_in_darkmode&sanitize=true" align=middle width=700.2745596pt height=173.33333654999998pt/></p>
  V(R) = T_{e} +D_{e}[1 - e^{\beta(r-r_{e})}]^2
  <p align="center"><img src="/tex/6b57e94b3ae7809bd3bc260fc76f0cc3.svg?invert_in_darkmode&sanitize=true" align=middle width=700.2745991999999pt height=119.90867954999999pt/></p> 
  V_{EMO}(R) = T_{e} + D_{e}[1 - exp(-\beta_{EMO}(R).(R-R_{e}))]^{2} 
  <p align="center"><img src="/tex/ef250ab0d4404e93c40c8058e17a2c10.svg?invert_in_darkmode&sanitize=true" align=middle width=652.33051455pt height=14.611878599999999pt/></p>
  \beta_{EMO}(R) = \sum_{i=0}^{N} \beta_{i} . y(R)^{i}
  <p align="center"><img src="/tex/919acf4ab9ca94d96f4680235b775b5f.svg?invert_in_darkmode&sanitize=true" align=middle width=593.0153229pt height=14.611878599999999pt/></p>
  y(R) = \frac{R^{p} - R_{e}^{p}}{R^{p} + R_{e}^{p}}
  <p align="center"><img src="/tex/c99bf2683851745507e9ec0aa1304b75.svg?invert_in_darkmode&sanitize=true" align=middle width=284.2016628pt height=35.251144499999995pt/></p>
V_{MLR}(R) = T_{e}
<p align="center"><img src="/tex/5556d3573fe4b8c9f75a6d8666f45328.svg?invert_in_darkmode&sanitize=true" align=middle width=1347.7692227999999pt height=8891.3242077pt/></p>
\chi^2 = \sum_{i=1}^{n} \left[ \frac{E_{i}^{obs} - E_{i}^{cal}(\mathbf{x})}{\sigma_i} \right]^2
<p align="center"><img src="/tex/3cde762b99ddde892badc3374179f7aa.svg?invert_in_darkmode&sanitize=true" align=middle width=700.274718pt height=93.51598905pt/></p>
E_{i}^{obs} - E_{i}^{cal}(\mathbf{x}^{(0)}) - \sum_{j=1}^{m} \frac{\partial{E_{i}^{cal}}}{\partial x_{j}} \bigg\rvert_{x_j=x_{j}^{(0)}} \Delta x_{j} = 0
<p align="center"><img src="/tex/623d22c997d14b289d678042f6a91726.svg?invert_in_darkmode&sanitize=true" align=middle width=701.5531330499999pt height=96.0566805pt/></p>
\min_{\mathbf{x}}\;\chi^2 = \frac{1}{\sigma^2} | \hat{A}\mathbf{x} - \mathbf{b} |^2
<p align="center"><img src="/tex/93c3b435ba275cc1e1dcf8a9839c7078.svg?invert_in_darkmode&sanitize=true" align=middle width=898.47997635pt height=634.1552646pt/></p>
\frac{\partial E}{\partial a} = \langle \Psi \vert \frac{\partial H}{\partial a}\vert \Psi \rangle
<p align="center"><img src="/tex/100fd68d56b71b7d7f5a844b4c200841.svg?invert_in_darkmode&sanitize=true" align=middle width=1056.97489095pt height=311.41552859999996pt/></p>
\frac{1}{\sigma_{k}^{2} + 0.3(E_{k}^{exp} - E_{k}^{calc})^{2}}
$$

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
