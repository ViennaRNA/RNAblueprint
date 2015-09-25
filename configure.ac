# RNAdesign: configure.ac

# initial information about the project
AC_INIT([RNAdesign],[0.1],[s.hammer@univie.ac.at],[RNAdesign],[http://www.tbi.univie.ac.at/~jango/])

# automake initialisation (mandatory)
AC_PREREQ([2.59])
AM_INIT_AUTOMAKE([1.9 subdir-objects])

# check if the source folder is correct
AC_CONFIG_SRCDIR([src/main.cc])

# include m4 local macros
AC_CONFIG_MACRO_DIR([m4])

# init doxygen
DX_HTML_FEATURE(ON)
DX_CHM_FEATURE(OFF)
DX_CHI_FEATURE(OFF)
DX_MAN_FEATURE(ON)
DX_RTF_FEATURE(OFF)
DX_XML_FEATURE(OFF)
DX_PDF_FEATURE(ON)
DX_PS_FEATURE(OFF)
DX_INIT_DOXYGEN($PACKAGE_NAME, doxygen.cfg, doc)

# check for C++ compiler
# AC_PROG_INSTALL
AC_PROG_CXX
AC_PROG_CC
AC_PROG_CPP
AC_PROG_LIBTOOL
AC_SUBST(LIBTOOL_DEPS)

# use the C++ compiler for the following checks
AC_LANG([C++])

# C++11 dependent program!
AX_CXX_COMPILE_STDCXX_11

# Require Boost Library
# https://raw.githubusercontent.com/tsuna/boost.m4/master/README
##########################################################################

BOOST_REQUIRE([1.55])
#BOOST_SYSTEM
#BOOST_FILESYSTEM

# Check lib boost_program_options
BOOST_PROGRAM_OPTIONS([s])
BOOST_GRAPH([s])
BOOST_TEST([s])
#AC_CHECK_LIB([boost_program_options],[_init])

# debug compilation support
##########################################################################

AC_MSG_CHECKING([whether to build with debug information])
AC_ARG_ENABLE([debug],
    [AS_HELP_STRING([--enable-debug],
        [enable debug data generation (def=no)])],
    [debugit="$enableval"],
    [debugit=no])
AC_MSG_RESULT([$debugit])

if test x"$debugit" = x"yes"; then
    AC_DEFINE([DEBUG],[],[Debug Mode])
    AM_CXXFLAGS="$AM_CXXFLAGS -g -Wall -Werror -Wno-uninitialized -O0"
else
    AC_DEFINE([NDEBUG],[],[No-debug Mode])
    AM_CXXFLAGS="$AM_CXXFLAGS -O3"
fi

# adding the RNA library and its dependencies (static name 'libRNA.a')
##########################################################################

# Check for OpenMP and prepare flags
AC_OPENMP
AC_SUBST(OPENMP_CFLAGS)

# Checks for libraries: math.
AC_CHECK_LIB([m], [exp],
	[LIBS="-lm $LIBS"],
	AC_MSG_ERROR([libMath not found!]))
	
# Checks for libraries: gmp.
AC_CHECK_LIB([gmp], [exp],
	[LIBS="-lgmp $LIBS"],
	AC_MSG_ERROR([libGMP not found!]))

# introduce the optional configure parameter for a non-standard install prefix of RNAlib
AC_ARG_WITH([RNAlib],
    [AS_HELP_STRING([--with-RNAlib=prefix],
        [try this for a non-standard install prefix of the RNA library])],
    [RNAlibPATHSET=1],
    [RNAlibPATHSET=0])
	
# if optional parameter used, extend path flags for compliler and linker
if test $RNAlibPATHSET = 1 ; then
    # extend the compiler and linker flags according to the path set
    AM_CXXFLAGS="$AM_CXXFLAGS -I$with_RNAlib/include"
    AM_LDFLAGS="$AM_LDFLAGS -L$with_RNAlib/lib64 -L$with_RNAlib/lib"
fi

# Checks for ViennaRNA library.
# LIBS="-lRNA $LIBS"
OLD_LDFLAGS=$LDFLAGS
LDFLAGS="-lgomp $LDFLAGS"
AC_CHECK_LIB([RNA],[fold],
	[LIBS="-lRNA $LIBS"],
	AC_MSG_ERROR([libRNA not found!]))
#LDFLAGS=$OLD_LDFLAGS

# Checks for header files.
##########################################################################
OLD_CXXFLAGS=$CXXFLAGS
CPPFLAGS="-std=c++11 $CPPFLAGS"

AC_HEADER_STDC
AC_CHECK_HEADERS([iostream])
AC_CHECK_HEADERS([sstream])
AC_CHECK_HEADERS([fstream])
AC_CHECK_HEADERS([iomanip])

AC_CHECK_HEADERS([string])
AC_CHECK_HEADERS([vector])
AC_CHECK_HEADERS([array])
AC_CHECK_HEADERS([set])
AC_CHECK_HEADERS([unordered_map])

AC_CHECK_HEADERS([functional])
AC_CHECK_HEADERS([limits])
AC_CHECK_HEADERS([chrono])
AC_CHECK_HEADERS([random])

AC_CHECK_HEADERS([boost/lexical_cast.hpp])
AC_CHECK_HEADERS([boost/config.hpp])
AC_CHECK_HEADERS([boost/functional/hash.hpp])
AC_CHECK_HEADERS([boost/iterator/reverse_iterator.hpp])
# AC_CHECK_HEADERS([boost/graph/ear_decomposition.hpp])

AC_CHECK_HEADERS([ViennaRNA/fold.h])
CPPFLAGS=$OLD_CPPFLAGS
# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_SIZE_T

# SWIG Interface
AX_PKG_SWIG(1.3.17, [], [ AC_MSG_ERROR([SWIG is required to build..]) ])
AX_SWIG_ENABLE_CXX
AX_SWIG_MULTI_MODULE_SUPPORT

# Add -fPIC for SWIG interface
AM_CXXFLAGS="-fPIC $AM_CXXFLAGS"

# Python/Perl support
AM_PATH_PYTHON(2.7, [], [ AC_MSG_ERROR([Python is required to build..]) ])
AC_PATH_PROGS(PERL, perl5 perl)

# distribute additional compiler and linker flags
# --> set these variables instead of CXXFLAGS or LDFLAGS
AC_SUBST([AM_CXXFLAGS])
AC_SUBST([AM_LDFLAGS])
AC_SUBST([LIBS])

# files to generate via autotools (.am or .in source files)
AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([lib/Makefile])
AC_CONFIG_FILES([src/Makefile])
AC_CONFIG_FILES([tests/Makefile])
AC_CONFIG_FILES([interfaces/Makefile])
AC_CONFIG_FILES([interfaces/Python/Makefile])
AC_CONFIG_FILES([interfaces/Python/setup.py])
AC_CONFIG_FILES([interfaces/Perl/Makefile])
AC_CONFIG_FILES([interfaces/Perl/Makefile.PL])

# generate the final Makefile etc.
AC_OUTPUT
