# SMFA

SMFA is a general program package for performing quantum chemistry calculations on large
molecules, using an energy-based fragmentation approach. The program can calculate
electronic energies, energy gradients and second derivatives; perform geometry
optimization; find first order saddle points (transition states); perform energy optimized
scans along a user-defined path; and evaluate various molecular properties. The program
can use any of the following quantum chemistry packages: GAMESS(US), GAUSSIAN,
NWChem and Q-Chem. In addition, SMFA provides a number of utility programs that, inter
alia, calculate vibrational frequencies and infrared spectra with isotopic substitutions, the
electrostatic potential on the solvent-accessible-surface, and isodesmic and higher order
near-iso-energetic reaction schemes. Calculations of the electronic energy and related
properties can be carried out using a scheme that provides a computation time that is
linearly dependent on the size of the molecule or, if the user has enough processing units
available, in a computation time that is independent of the size of the molecule.

### Table of contents:

* [Requirements](#requirements)
* [Installation](#installation)
* [SMFA Publications](#smfa-publications)
* [Examples](/doc/testcases)
* [Licensing](#licensing)

## Requirements
* fortran compiler
* [cmake](https://cmake.org/)
* [perl](https://www.perl.org/)
* Quantum Chemistry programs
    - [GAMESS](http://www.msg.ameslab.gov/gamess/)
    - [Gaussian](http://gaussian.com/)
    - [NWChem](http://www.nwchem-sw.org/)
    - [QChem](http://www.q-chem.com/)
* Optional: [DALTON](http://daltonprogram.org/)

## Installation

git clone https://github.com/SMFA/SMFA.git      # Clone SMFA source code from GitHub

## SMFA Publications
1. Collins, M. A. Physical Chemistry Chemical Physics 2012, 14, 7744–7751.

## Licensing


## Version
1.0rc1
