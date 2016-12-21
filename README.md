
# RNAblueprint library

## Introduction

The RNAblueprint library solves the problem of stochastically sampling RNA/DNA sequences
compatible to multiple structural constraints.
It only creates sequences that fulfill all base pairs specified in any given input structure.
Furthermore, it is possible to specify sequence constraints in IUPAC notation.
Solutions are sampled uniformly from the whole solution space, therefore it is guaranteed,
that there is no bias towards certain sequences.

The library is written in C++ with SWIG scripting interfaces for Python and Perl.
Please cite the software as specified at the bottom of the page!

## Dependencies

### Required:

 * GNU Automake
 * Boost Graph Library
 * C++ Standard Library

### Optional:

 * Boost Programm Options (default: on)
 * SWIG for interfaces (default: on)
 * Python for interface (default: on)
 * Perl for interface (default: on)
 * ExtUtils::Embed module for perl interface (default: on)
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

In case of a local installation, please do not forget to adopt your path variables such as
`PATH`, `LD_LIBRARY_PATH`, `CPLUS_INCLUDE_PATH`, `PYTHONPATH`, `PERL5LIB`

### Most important configure options are:

    --prefix Specify an installation path prefix
    --with-boost Specify the installation directory of the boost library
    --disable-program  Disable RNAblueprint program compilation
    --disable-swig Disable all SWIG scripting interfaces
    --enable-libGMP Enable the calculation of big numbers with multiprecision
    --disable-openmp Disable the usage of parallel computation

TIP: You might want call `./configure --help` for all install options!

## Documentation

Library documentation is done using Doxygen. Call `make doxygen-doc` for a offline version in HTML and PDF.
There is also a online version available here: [http://ribonets.github.io/RNAblueprint/](http://ribonets.github.io/RNAblueprint/)

Documentation of the RNAblueprint program is provided by calling `RNAblueprint --help`.

## Testing

Unit tests are available for many functions of the library. Please call `make check` to run these tests!

## How to cite

Stefan Hammer, Birgit Tschiatschek, Christoph Flamm, Ivo L. Hofacker and Sven Findeiß. Submitted 2016. "Flexible multiple target nucleic acid sequence design" Journal: Pages. doi: DOI.

Publication for this software package is in progress.
We will update this information as soon as a preprint is available online!
