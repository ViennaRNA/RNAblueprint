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
        AC_MSG_ERROR([libGMP not found!]))
    AC_DEFINE([LIBGMP],[1],[Using libGMP])
    [libgmpit=true]
])
AM_CONDITIONAL([ENABLE_LIBGMP], [test x"$libgmpit" = x"true"])

# Checks for header files.
##########################################################################
OLD_CXXFLAGS=$CXXFLAGS
CPPFLAGS="-std=c++11 $CPPFLAGS"

AC_CHECK_HEADERS([algorithm])
AC_CHECK_HEADERS([array])
AC_CHECK_HEADERS([chrono])
AC_CHECK_HEADERS([cctype])
AC_CHECK_HEADERS([exception])
AC_CHECK_HEADERS([fstream])
AC_CHECK_HEADERS([functional])
AC_CHECK_HEADERS([iomanip])
AC_CHECK_HEADERS([iostream])
AC_CHECK_HEADERS([iterator])
AC_CHECK_HEADERS([limits])
AC_CHECK_HEADERS([list])
AC_CHECK_HEADERS([random])
AC_CHECK_HEADERS([regex])
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
AX_PKG_SWIG(3.0.7, [], [ AC_MSG_ERROR([SWIG (>= 3.0.7) is required to build.]) ])
AX_SWIG_ENABLE_CXX
AX_SWIG_MULTI_MODULE_SUPPORT

# distribute additional compiler and linker flags
# set these variables instead of CXXFLAGS or LDFLAGS
AC_SUBST([AM_CXXFLAGS])
AC_SUBST([AM_LDFLAGS])
AC_SUBST([LIBS])

# Python support
##########################################################################
AC_ARG_ENABLE([PYTHON],
    [AS_HELP_STRING([--disable-python],
        [Disable the Python interface.])])

AS_IF([test "x$enable_python" != "xno"], [
    AC_MSG_NOTICE([Building with Python bindings])
    AM_PATH_PYTHON(2.7, [], [ AC_MSG_ERROR([Python is required to build.]) ])
    AX_PYTHON_DEVEL
    AX_SWIG_PYTHON
])

AM_CONDITIONAL([ENABLE_PYTHON],[test "x$enable_python" != "xno"])

# Perl5 support
##########################################################################
AC_ARG_ENABLE([PERL5],
    [AS_HELP_STRING([--disable-perl],
        [Disable the Perl5 interface.])])

AS_IF([test "x$enable_perl" != "xno"], [
    AC_MSG_NOTICE([Building with Perl bindings])
    AX_PERL_EXT([$prefix])
    if test x"$PERL" = x; then
        AC_MSG_ERROR([Perl is required to build.])
    fi
    AX_PERL_EXT_FLAGS([PERLXS_CFLAGS], [PERLXS_LDFLAGS])
    AC_SUBST([PERLXS_CFLAGS])
    AC_SUBST([PERLXS_LDFLAGS])
])

AM_CONDITIONAL([ENABLE_PERL5],[test "x$enable_perl" != "xno"])

# files to generate via autotools (.am or .in source files)
##########################################################################
AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([lib/Makefile])
AC_CONFIG_FILES([src/Makefile])
AC_CONFIG_FILES([tests/Makefile])
AC_CONFIG_FILES([interfaces/Makefile])
AM_COND_IF([ENABLE_PYTHON], [
    AC_CONFIG_FILES([interfaces/Python/Makefile])
])
AM_COND_IF([ENABLE_PERL5], [
    AC_CONFIG_FILES([interfaces/Perl/Makefile])
])

# generate the final Makefile etc.
##########################################################################
AC_OUTPUT
