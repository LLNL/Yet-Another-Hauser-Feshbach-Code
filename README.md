# YAHFC

YAHFC (Yet Another hauser Feshbach Code) computer program framework to model 
low-energy nuclear reactions. The framework allows for reactions with incident 
particles including gammas, protons, neutrons, deuterons, tritons, helium-3,
and alphas. The code can model nuclear reactions or initial populations
specified by the user. The decays principally follow the compound-nucleus 
hypothesis and the statistical decay of Hauser and Feshbach (Phys. Rev. **87**, 366
(1952). The code models nuclear decays with a Monte Carlo process that tracks 
the decay of each state. Data libraries are produced that can be translated into
GNDS (generalized Nuclear Data Structure). In addition, YAHFC can be used as an
event generator, with individual decays written to disk for off-line analysis. 
Both serial and MPI versions are provided.

# Installation

YAHFC code system can can be installed from the GitLab repository. 

## YAHFC enviornment variables  

YAHFC executes in a unix enevironment and requires the following environment 
variables pointing to directories accessed during execution:

YAHFC_DIR = _The upper level directory where YAHFC is installed_

YAHFC_DATA = $YAHFC_DIR/Data

YAHFC_BIN_DIR = _The directory where the executable YAHFC.x is to be installed_

$YAHFC_BIN_DIR should also be placed in the users path.

## Building YAHFC

The source code is in the directory $YAHFC_DIR/Src and a makefile is provided to build the executable YAHFC.x. The makefile should be edited to set the FORTRAN 90 compiler (F90COMPILER) and the compiler flags (F90CFLAGS) and link flags (LFLAGS). The default compiler set in the makefile is gfortran for a serial build and mpifort for an mpi build. Typing "make YAHFC" will build the executable YAHFC.x and copy it to the directory $YAHFC_BIN_DIR. Typing "make YAHFC-MPI" will build the MPI executable YAHFC-MPI.x and copy it to the $YAHFC_BIN_DIR directory.

## FRESCOX

YAHFC also requires the installation of [Frescox](https://github.com/LLNL/Frescox). The executable frescox needs to be placed in a directory within the users path, for example $YAHFC_BIN_DIR above. 

# Documentation

A manual in latex and pdf format is provided in the directory $YAHFC_DIR/Documentation.

# Examples

Several example script files are provided in the subdirectories found $YAHFC_DIR/YAHFC-Examples.

# Libraries

A set of data files are created in the run directory. These can be translated into GNDS using the Python program parseYahfc.py in the directory $YAHFC_DIR/GNDS. The GNDS parser requires that FUDGE (v4.3), which will be released when GNDS-2.0 is adopted.

# External Packages

YAHFC makes use of data files from the Reference Input Parameter Library
(RIPL-3, Nucl. Data Sheets **110**, 3107 (2009)) contained in the directory $YAHFC_DIR/Data. In addition, YAHFC includes two 
subroutines from the EISPACK library that have been modified and updated 
to FORTRAN 90. 

# Author

YAHFC was created by Erich Ormand at LLNL. 

# License

YAHFC is distributed un the terms of the MIT license.

SPDX-License-identifier: MIT

LLNL-CODE-XXXXXX
