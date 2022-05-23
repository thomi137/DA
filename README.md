# DA
Diploma thesis in physics ETH Zurich

## Purpose
This repository comprises the source code used to simulate kicked Bose-Einstein condensates trapped in a lattice with harmonic potential wells. I specifically is used to examine behaviour of the governing Gross-Pitaevskii equation on specific Kick strengths. The non-linearity implies the validity of the KAM theory of classical chaos theory.

At the time it was thought that this gives rise to chaotic behaviour in the quantum realm

## Installation Requirements
In order to be able to compile the file in this repository, you will need to have the following installed:
    
- A working C++ compiler with the standard library
- A working Fortran compiler if you are installing fftw from scratch.
- The Boost extension to the standard library. 
  This can be obtained via your pakage manager, 
  homebrew if you are on mac or simply [here](https://www.boost.org/users/download/)
- Lapack and BLAS: Installation and usage docs founr [here](http://www.netlib.org/lapack/)
- The [Fastest Fourier Transform in the West] (https://www.fftw.org/) Library

**Note:** Some of the files make use of a library called IETL (the Iterative Eigensolver Template Library).
To the best of my knowledge, it is no longer under active development or has been merged into some other project.
To solve the Gross-Pitaevskii equation however, it is not needed, so most of the programs here will compile.

Also, I am in the process of converting some of the makefiles to cmake, which should make it easier to compile modern C/C++