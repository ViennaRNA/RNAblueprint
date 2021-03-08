# RNAblueprint: configure.ac

# initial information about the project
AC_INIT([RNAblueprint],[1.3.1],[sven@bioinf.uni-leipzig.de],[RNAblueprint],[https://github.com/ViennaRNA/RNAblueprint])

# automake initialisation (mandatory)
AC_PREREQ([2.59])
AM_SILENT_RULES([yes])
AM_INIT_AUTOMAKE([1.9 subdir-objects])

# check if the source folder is correct
AC_CONFIG_SRCDIR([src/main.cc])

# include m4 local macros
AC_CONFIG_MACRO_DIR([m4])

# init doxygen
##########################################################################
DX_DOXYGEN_FEATURE(OFF)
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
AC_ARG_ENABLE([tests],
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

# build with RNAblueprint binary
##########################################################################tests
AC_ARG_ENABLE([program],
    [AS_HELP_STRING([--disable-program],
        [Disable building the RNAblueprint program.])])

AS_IF([test "x$enable_program" != "xno"], [
    AC_MSG_NOTICE([Building with RNAblueprint program])
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
#AC_OPENMP
#AC_SUBST(OPENMP_CXXFLAGS)

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
AC_ARG_ENABLE([swig],
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
    # Python2 support
    ######################################################################
    AX_REQUIRE_DEFINED([AX_PYTHON2_DEVEL])
    AC_ARG_ENABLE([python2],
        [AS_HELP_STRING([--disable-python2],
            [Disable the Python 2.x interface.])])

    AS_IF([test "x$enable_python2" != "xno"], [
        AC_MSG_NOTICE([Building with Python 2.x bindings])
        ## check for python2 config
        AX_PYTHON2_DEVEL

        if test "x$python2_enabled_but_failed" != "x"
        then
            with_python2="no"
        else
            AC_SUBST(PYTHON2DIR,$python2dir)
            AC_SUBST(PKGPYTHON2DIR,$pkgpython2dir)
            AC_SUBST(PYEXEC2DIR,$py2execdir)
            AC_SUBST(PKGPYEXEC2DIR,$pkgpy2execdir)

            AC_DEFINE([WITH_PYTHON2_INTERFACE], [1], [Create the python2 interface to RNAblueprint])
            AC_SUBST([PYTHON2_INTERFACE], [Python2])
            AC_CONFIG_FILES([interfaces/Python2/Makefile interfaces/Python2/version.i])
        fi
        [enable_python2_status="yes"]
    ], [
        [enable_python2_status="no"]
    ])

    # Python3 support
    ######################################################################
    AX_REQUIRE_DEFINED([AX_PYTHON3_DEVEL])
    AC_ARG_ENABLE([python],
        [AS_HELP_STRING([--disable-python],
            [Disable the Python 3.x interface.])])

    AS_IF([test "x$enable_python" != "xno"], [
        AC_MSG_NOTICE([Building with Python 3.x bindings])
        ## check for python3 config
        AX_PYTHON3_DEVEL

        if test "x$python3_enabled_but_failed" != "x"
        then
            with_python="no"
        else
            AC_DEFINE([WITH_PYTHON3_INTERFACE], [1], [Create the Python3 interface to RNAlib])
            AC_SUBST([PYTHON3_INTERFACE], [Python3])
        fi

        AC_CONFIG_FILES([interfaces/Python3/Makefile interfaces/Python3/version.i])
        [enable_python_status="yes"]
    ], [
        [enable_python_status="no"]
    ])

    # Perl5 support
    ######################################################################
    AC_ARG_ENABLE([perl],
        [AS_HELP_STRING([--disable-perl],
            [Disable the Perl5 interface.])])

    AS_IF([test "x$enable_perl" != "xno"], [
        AC_MSG_NOTICE([Building with Perl bindings])
        AX_PERL_EXT
        if test "x$PERL" = "x"; then
            AC_MSG_ERROR([Perl is required to build.])
            [enable_perl_status="Perl is required to build."]
        fi
    ])

    AS_IF([test "x$enable_perl" != "xno"], [
        # Compose the correct installation path for perl modules
        #
        # here we actually have to account for INSTALLDIRS env variable, which can be
        #
        # site    = where the local systems administrator installs packages to
        # vendor  = where system packages are installed to, or
        # core    = where perl core packages are installed
        #
        # The default selection is 'site', but upon packaging for a specific distribution
        # we might want the user to set this to 'vendor'
        #
        AS_IF([ test "x$INSTALLDIRS" == "xvendor" ],[
          PERL_ARCH_RELATIVE_INSTALL_DIR=`echo ${PERL_EXT_VENDORARCH} | sed "s,${PERL_EXT_VENDORPREFIX},,"`
          PERL_LIB_RELATIVE_INSTALL_DIR=`echo ${PERL_EXT_VENDORLIB} | sed "s,${PERL_EXT_VENDORPREFIX},,"`
          ],[
          PERL_ARCH_RELATIVE_INSTALL_DIR=`echo ${PERL_EXT_SITEARCH} | sed "s,${PERL_EXT_SITEPREFIX},,"`
          PERL_LIB_RELATIVE_INSTALL_DIR=`echo ${PERL_EXT_SITELIB} | sed "s,${PERL_EXT_SITEPREFIX},,"`
        ])
        AC_SUBST(PERL_ARCH_RELATIVE_INSTALL_DIR)
        AC_SUBST(PERL_LIB_RELATIVE_INSTALL_DIR)

        AC_DEFINE([WITH_PERL_INTERFACE], [1], [Create the perl interface to RNAlib])
        AC_SUBST([PERL_INTERFACE], [Perl])
        AC_CONFIG_FILES([interfaces/Perl/Makefile interfaces/Perl/version.i])

        [enable_perl_status="yes"]
    ], [
        [enable_perl_status="no"]
    ])
])

# Configure supported platforms
##########################################################################
case $host in
  *-apple-darwin*)
    AC_DEFINE([OSX], [], [Define if building on Apple OSX.]) ;;
esac

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
  RNAblueprint program: $enable_program_status
  Debug information:    $enable_debug
  Use GMP integers:     $enable_libgmp_status

SWIG enabled:           $enable_swig_status
  Perl interface:       $enable_perl_status
  Python2 interface:    $enable_python2_status
  Python3 interface:    $enable_python_status

Documentation enabled:  $doxygen_doc
  (HTML):               $doxygen_html
  (PDF):                $doxygen_pdf
  (MAN):                $doxygen_man
  (PS):                 $doxygen_ps

Unit Tests enabled:     $enable_tests_status
])
