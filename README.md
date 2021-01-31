# tlwall
## Table of contents
 - [Introduction](#introduction)
 - [Installation](#installation)
 - [Tasks performed by tlwall](#tasks-performed-by-tlwall)
 - [Structure of the package](#structure-of-the-package)
 - [Installation code](#installation-code)
 - [References](#references)

## Introduction
"tlwall" is a CERN python code which use the transmission line theory to calculate the resistive wall impedance.
The first version of tlwall has been made in 2011 in matlab by Carlo Zannini, in 2013
started the python developement.
The current tlwall version has undergone a major restyling to follow current CERN guidelines and pep8 recommendations.

The aim of the tlwall is to calculate
- the longitudinal impedance
- the transverse impedance
- the driving and detuning impedance taking care of the form factors (yokoya form factors)
- the space charge impedance (for speed < c)
- the surface impedance

The transmission line equations can be applied recursively to take into account whatever number of layers.
The beam speed and the roughness are taken into account in the formulas.

It is also possible to calculate the impedance of a list of chambers to "build" an entire accelerator.

## Installation
```bash
pip install pypkgexample
```
assuming that the package has been downloaded in the current folder.

## Tasks performed by tlwall


## Structure of the package

```tlwall``` has the structure used for standard python packages. It consists in a folder named after the package "tlwall" (which is also the top level of the git repository) 
that contains the source code, the required code and information for the installation, documentation, unit tests and usage examples. In particular:
 - The **source code** is contained in a subfolder that also has the same name of the python package (tlwall).
 - **Unit tests** is contained in the folder "tests" 
 - **Examples** illustrating the package usage are hosted in the folder "examples"
 - **License** information is contained in the file "LICENSE.txt"
 - This **documentation** is contained in the file "README.md"
 - The **installation process** is defined by the files ["pyproject.toml"](#pyprojecttoml), ["MANIFEST.in"](#manifestin), and  ["setup.py"](#setuppy), which will be described in more detail in the following section.
 
## Installation code

The build and installation process is performed by the pip pacakge installer based on the following files:

### pyproject.toml

The file "pyproject.toml" defines the backend used of the build, which in our case is [setuptools](https://pypi.org/project/setuptools/), and the dependencies that are required to build the package. 


### MANIFEST.in

The file "MANIFEST.in" defines the additional files that need to be copied together with the installed package, together with those that are strictly required for the package to work. In our case we include, "pyproject.toml", this readme file and the license information.

```python
include pyproject.toml

# Include the README
include *.md

# Include the license file
include LICENSE.txt
```

### setup.py


## References
The following resources were used for preparing this package and documentation.



