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
* [DALTON](http://daltonprogram.org/)
* At least one of the following Quantum Chemistry programs:
    - [GAMESS](http://www.msg.ameslab.gov/gamess/)
    - [Gaussian](http://gaussian.com/)
    - [NWChem](http://www.nwchem-sw.org/)
    - [QChem](http://www.q-chem.com/)


The SMFA_Users_Guide.pdf (Section 2.7) in SMFAPAC/doc contains information
about how to get SMFA to "load" each (at least one) of these quantum chemistry
packages.

Several Perl modules are also required, and these can be installed 
with the following commands (administrator privileges required):

```shell
#> sudo cpan App::cpanminus
#> sudo cpanm Shell
#> sudo apt-get install libncurses5-dev 
#> sudo cpanm Curses
```
One of the optional utility programs in SMFA requires the openbabel program
(see Section 5.2 in SMFA_Users_Guide.pdf), so you will need to install
openbabel if you want to use this feature. You can install openbabel at any
time (the build below does not require it).




## Installation

```shell
#> git clone https://github.com/SMFA/SMFA.git
#> mkdir build
#> cd build
#> cmake ../
#> make install
```

'make install' compiles the binaries and also moves them to the
SMFAPAC/exe directory.  If you need to rebuild, simply remove all
files in the build directory and `make install` again.

The SMFAPAC/bin directory must be on the user's path, which can be achieved by
adding the following to the appropriate rc file: 

For ~/.cshrc:

```shell
># set path = ( $path /path/to/SMFAPAC/bin)
```

For ~/.bashrc:

```shell
># export PATH=$PATH:/path/to/SMFAPAC/bin
```

Ensure /path/to is replaced with the actual path to the directory.




## SMFA Publications

1. The user guide can be found in doc/SMFA_Users_Guide.pdf and 
   contains detailed instructions on how to use the package.
   
2. Collins, M. A. Physical Chemistry Chemical Physics 2012, 14, 7744–7751.

## Licensing


## Version
1.0rc1
