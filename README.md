
# RNAblueprint library written in C++

## Introduction

The RNAblueprint library still needs some documentation! Write here something about the theoretical background Show how to cite the software
Dependencies

### Required:

 * GNU Automake
 * Boost Graph Library
 * Boost Programm Options
 * C++ Standard Library

### Optional:

 * SWIG for interfaces
 * Python for interface
 * Perl for interface
 * Doxygen for documentation
 * LaTeX for PDF documentation
 * libGMP for multiprecision integers
 * openmp for parallel computation
 * Boost Unit Test Framework

## Installation

### Just call these commands:

```bash
./autogen.sh
./configure
make
make install
```

### Most important configure options are:

    --prefix Specify an installation path prefix
    --with-boost Specify the installation directory of the boost library
    --disable-swig Disable all SWIG scripting interfaces
    --disable-program  Disable RNAblueprint program compilation
    --enable-libGMP Enable the calculation of big numbers with multiprecision
    --disable-openmp Disable the usage of parallel computation

TIP: You might want call ./configure –help for all install options!

## Documentation

Documentation is done using Doxygen. Call 'make doc' for a offline version in HTML and PDF.
There is also a online version available here: [http://ribonets.github.io/RNAblueprint/](http://ribonets.github.io/RNAblueprint/)

## How to cite

This is a early release for testing purposes only and a publication for this software package is in progress.
We will update this information as soon as a preprint is available online!
