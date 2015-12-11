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
##########################################################################
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
##########################################################################
# AC_PROG_INSTALL
AC_PROG_CXX
AC_PROG_CC
AC_PROG_CPP
AC_PROG_LIBTOOL
LT_INIT([pic-only, static])
AC_SUBST(LIBTOOL_DEPS)

# use the C++ compiler for the following checks
AC_LANG([C++])

# C++11 dependent program!
AX_CXX_COMPILE_STDCXX_11([noext], [mandatory])

# Require Boost Library
# https://raw.githubusercontent.com/tsuna/boost.m4/master/README
##########################################################################

BOOST_REQUIRE([1.55])
#BOOST_SYSTEM
#BOOST_FILESYSTEM

# Check lib boost_program_options, graph and regex
BOOST_GRAPH([s])

# build with unit tests
##########################################################################tests
AC_ARG_ENABLE([UNIT_TESTS],
    [AS_HELP_STRING([--disable-tests],
        [Disable support for Boost UNIT tests.])])

AS_IF([test "x$enable_tests" != "xno"], [
    AC_MSG_NOTICE([Building with Boost UNIT tests])
    BOOST_TEST([s])
    [enable_tests_status="yes"]
], [
    [enable_tests_status="no"]
])

AM_CONDITIONAL([ENABLE_TESTS],[test "x$enable_tests_status" = x"yes"])

# build with RNAdesign binary
##########################################################################tests
AC_ARG_ENABLE([RNADESIGN_PROGRAM],
    [AS_HELP_STRING([--disable-program],
        [Disable building the RNAdesign program.])])

AS_IF([test "x$enable_program" != "xno"], [
    AC_MSG_NOTICE([Building with RNAdesign program])
    BOOST_PROGRAM_OPTIONS([s])
    [enable_program_status="yes"]
], [
    [enable_program_status="no"]
])

AM_CONDITIONAL([ENABLE_PROGRAM],[test "x$enable_program_status" = x"yes"])

# debug compilation support
##########################################################################

AC_MSG_CHECKING([whether to build with debug information])
AC_ARG_ENABLE([debug],
    [AS_HELP_STRING([--enable-debug],
        [enable debug data generation (def=no)])],
    [enable_debug="$enableval"],
    [enable_debug=no])
AC_MSG_RESULT([$enable_debug])

if test x"$enable_debug" = x"yes"; then
    AC_DEFINE([DEBUG],[],[Debug Mode])
    AM_CXXFLAGS="$AM_CXXFLAGS -g -Wall -Werror -Wno-uninitialized -O0"
else
    AC_DEFINE([NDEBUG],[],[No-debug Mode])
    AM_CXXFLAGS="$AM_CXXFLAGS -O3"
fi

# Check for OpenMP enabled and prepare flags
##########################################################################
AC_OPENMP
AC_SUBST(OPENMP_CFLAGS)

# Checks for libraries: gmp enabled?
##########################################################################
AC_ARG_ENABLE([libGMP],
    [AS_HELP_STRING([--enable-libGMP],
        [Enable this to calculate everything using boost multiprecision integers])])

AS_IF([test x"$enable_libGMP" = x"yes"], [
    AC_CHECK_LIB([gmp], [__gmp_printf],
        [LIBS="-lgmp $LIBS"],
        [enable_libgmp_status="libGMP not found!"]
        AC_MSG_ERROR([libGMP not found!]))
    AC_DEFINE([LIBGMP],[1],[Using libGMP])
    [enable_libgmp_status="yes"]
], [
    [enable_libgmp_status="no"]
])
AM_CONDITIONAL([ENABLE_LIBGMP], [test x"$enable_libgmp_status" = x"yes"])

# Checks for header files.
##########################################################################
OLD_CXXFLAGS=$CXXFLAGS
CPPFLAGS="-std=c++11 $CPPFLAGS"

AC_CHECK_HEADERS([algorithm])
AC_CHECK_HEADERS([array])
AC_CHECK_HEADERS([chrono])
AC_CHECK_HEADERS([cctype])
AC_CHECK_HEADERS([cstddef])
AC_CHECK_HEADERS([exception])
AC_CHECK_HEADERS([fstream])
AC_CHECK_HEADERS([functional])
AC_CHECK_HEADERS([iomanip])
AC_CHECK_HEADERS([iostream])
AC_CHECK_HEADERS([iterator])
AC_CHECK_HEADERS([limits])
AC_CHECK_HEADERS([list])
AC_CHECK_HEADERS([random])
AC_CHECK_HEADERS([set])
AC_CHECK_HEADERS([sstream])
AC_CHECK_HEADERS([string])
AC_CHECK_HEADERS([unordered_map])
AC_CHECK_HEADERS([unordered_set])
AC_CHECK_HEADERS([utility])
AC_CHECK_HEADERS([vector])

AC_CHECK_HEADERS([boost/config.hpp])
AC_CHECK_HEADERS([boost/functional/hash.hpp])
#AC_CHECK_HEADERS([boost/graph/adjacency_list.hpp])
#AC_CHECK_HEADERS([boost/graph/biconnected_components.hpp])
#AC_CHECK_HEADERS([boost/graph/bipartite.hpp])
#AC_CHECK_HEADERS([boost/graph/breadth_first_search.hpp])
#AC_CHECK_HEADERS([boost/graph/connected_components.hpp])
#AC_CHECK_HEADERS([boost/graph/graphml.hpp])
#AC_CHECK_HEADERS([boost/graph/iteration_macros.hpp])
#AC_CHECK_HEADERS([boost/graph/random_spanning_tree.hpp])
#AC_CHECK_HEADERS([boost/graph/undirected_dfs.hpp])
AC_CHECK_HEADERS([boost/lexical_cast.hpp])
#AC_CHECK_HEADERS([boost/program_options.hpp])
AC_CHECK_HEADERS([boost/property_map/property_map.hpp])
AC_CHECK_HEADERS([boost/property_map/shared_array_property_map.hpp])
AC_CHECK_HEADERS([boost/property_map/vector_property_map.hpp])
#AC_CHECK_HEADERS([boost/test/unit_test.hpp])

AM_COND_IF([ENABLE_LIBGMP], [
    AC_CHECK_HEADERS([boost/multiprecision/gmp.hpp])
    AC_CHECK_HEADERS([boost/multiprecision/random.hpp])
])

CPPFLAGS=$OLD_CPPFLAGS
# Checks for typedefs, structures, and compiler characteristics.
##########################################################################
AC_TYPE_SIZE_T

# SWIG Interface
##########################################################################
AC_ARG_ENABLE([SWIG],
    [AS_HELP_STRING([--disable-swig],
        [Disable all SWIG interfaces.])])

AS_IF([test "x$enable_swig" != "xno"], [
    AX_PKG_SWIG(3.0.7, [
        [enable_swig_status="yes"]
    ], [ 
        AC_MSG_ERROR([SWIG (>= 3.0.7) is required to build.])
        [enable_swig_status="SWIG (>= 3.0.7) is required to build."]
    ])
], [
    [enable_swig_status="no"]
])
AM_CONDITIONAL([ENABLE_SWIG], [test x"$enable_swig_status" = x"yes"])

AM_COND_IF([ENABLE_SWIG], [
    AX_SWIG_ENABLE_CXX
    AX_SWIG_MULTI_MODULE_SUPPORT

    # Python support
    ######################################################################
    AC_ARG_ENABLE([PYTHON],
        [AS_HELP_STRING([--disable-python],
            [Disable the Python interface.])])

    AS_IF([test "x$enable_python" != "xno"], [
        AC_MSG_NOTICE([Building with Python bindings])
        AM_PATH_PYTHON(2.7, [], [ 
            AC_MSG_ERROR([Python is required to build.])
            [enable_python_status="Python is required to build."]
        ])
        AX_PYTHON_DEVEL
        AX_SWIG_PYTHON
        [enable_python_status="yes"]
    ], [
        [enable_python_status="no"]
    ])
    
    # Perl5 support
    ######################################################################
    AC_ARG_ENABLE([PERL5],
        [AS_HELP_STRING([--disable-perl],
            [Disable the Perl5 interface.])])

    AS_IF([test "x$enable_perl" != "xno"], [
        AC_MSG_NOTICE([Building with Perl bindings])
        AX_PERL_EXT([$prefix])
        if test x"$PERL" = x; then
            AC_MSG_ERROR([Perl is required to build.])
            [enable_python_status="Perl is required to build."]
        fi
        AX_PERL_EXT_FLAGS([PERLXS_CFLAGS], [PERLXS_LDFLAGS])
        AC_SUBST([PERLXS_CFLAGS])
        AC_SUBST([PERLXS_LDFLAGS])
        [enable_perl_status="yes"]
    ], [
        [enable_perl_status="no"]
    ])
])

AM_CONDITIONAL([ENABLE_PYTHON],[test "x$enable_python_status" = x"yes"])
AM_CONDITIONAL([ENABLE_PERL5],[test "x$enable_perl_status" = x"yes"])

# distribute additional compiler and linker flags
# set these variables instead of CXXFLAGS or LDFLAGS
##########################################################################
AC_SUBST([AM_CXXFLAGS])
AC_SUBST([AM_LDFLAGS])
AC_SUBST([LIBS])

# files to generate via autotools (.am or .in source files)
##########################################################################
AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([lib/Makefile])
AC_CONFIG_FILES([src/Makefile])
AC_CONFIG_FILES([tests/Makefile])
AC_CONFIG_FILES([interfaces/Makefile])
AC_CONFIG_FILES([interfaces/Python/Makefile])
AC_CONFIG_FILES([interfaces/Perl/Makefile])

# generate the final Makefile etc. and display summary
##########################################################################
AC_OUTPUT

AS_IF([test x"$DX_FLAG_doc" = x"1"], [doxygen_doc="yes"],[doxygen_doc="no"])
AS_IF([test x"$DX_FLAG_html" = x"1"], [doxygen_html="yes"],[doxygen_html="no"])
AS_IF([test x"$DX_FLAG_pdf" = x"1"], [doxygen_pdf="yes"],[doxygen_pdf="no"])
AS_IF([test x"$DX_FLAG_man" = x"1"], [doxygen_man="yes"],[doxygen_man="no"])
AS_IF([test x"$DX_FLAG_ps" = x"1"], [doxygen_ps="yes"],[doxygen_ps="no"])

AC_MSG_NOTICE(
[
------------------------------------------------
Configure successful with the following options:

Options:
  Prefix path:          $prefix
  Custom BOOST path:    $with_boost
  RNAdesign program:    $enable_program_status
  Debug information:    $enable_debug
  Use GMP integers:     $enable_libgmp_status

SWIG enabled:           $enable_swig_status
  Perl interface:       $enable_perl_status
  Python interface:     $enable_python_status

Documentation enabled:  $doxygen_doc
  (HTML):               $doxygen_html
  (PDF):                $doxygen_pdf
  (MAN):                $doxygen_man
  (PS):                 $doxygen_ps

Unit Tests enabled:     $enable_tests_status
])
