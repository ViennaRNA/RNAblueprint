# RNAdesign: configure.ac

# initial information about the project
AC_INIT([RNAdesign],[0.1],[s.hammer@univie.ac.at])

# check if the source folder is correct
AC_CONFIG_SRCDIR([src/main.cc])

# include m4 local macros
AC_CONFIG_MACRO_DIRS([m4])

# check for C++ compiler
AC_PROG_CXX
AC_PROG_CC

# automake initialisation (mandatory) and check for minimal automake API version 1.9
AM_INIT_AUTOMAKE([1.9])

# use the C++ compiler for the following checks
AC_LANG([C++])

# C++11 dependent program!
AX_CXX_COMPILE_STDCXX_11([noext],[mandatory])
CPPFLAGS="-std=c++11 $CPPFLAGS"

# add standard library
LIBS="-lstdc++ $LIBS"

# add boost_program_options
LIBS="-lboost_program_options $LIBS"

# Checks for header files.
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
AC_CHECK_HEADERS([boost/graph/subgraph.hpp])
AC_CHECK_HEADERS([boost/graph/adjacency_list.hpp])
AC_CHECK_HEADERS([boost/graph/connected_components.hpp])
AC_CHECK_HEADERS([boost/graph/biconnected_components.hpp])
AC_CHECK_HEADERS([boost/graph/breadth_first_search.hpp])
AC_CHECK_HEADERS([boost/graph/random_spanning_tree.hpp])
AC_CHECK_HEADERS([boost/functional/hash.hpp])
AC_CHECK_HEADERS([boost/iterator/reverse_iterator.hpp])
AC_CHECK_HEADERS([boost/program_options.hpp])
AC_CHECK_HEADERS([boost/graph/graphml.hpp])
AC_CHECK_HEADERS([boost/graph/graph_utility.hpp])
AC_CHECK_HEADERS([boost/graph/iteration_macros.hpp])

# for test framework
AC_CHECK_HEADERS([boost/test/included/unit_test.hpp])

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_SIZE_T

##########################################################################
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

##########################################################################
# adding the RNA library (with static name 'libRNA.a')
##########################################################################

# adding the lib to the files to link
LIBS="-lRNA -lm $LIBS"

# introduce the optional configure parameter for a non-standard install prefix of XXX
AC_ARG_WITH([RNAlib],
    [AS_HELP_STRING([--with-RNAlib=prefix],
        [try this for a non-standard install prefix of the RNA library])],
    [RNAlibPATHSET=1],
    [RNAlibPATHSET=0])
	
# if optional parameter used, extend path flags for compliler and linker
if test $RNAlibPATHSET = 1 ; then
    # extend the compiler and linker flags according to the path set
    AM_CXXFLAGS="$AM_CXXFLAGS -I$with_RNAlib/include"
    AM_LDFLAGS="$AM_LDFLAGS -L$with_RNAlib/lib"
fi

##########################################################################

# distribute additional compiler and linker flags
# --> set these variables instead of CXXFLAGS or LDFLAGS
AC_SUBST([AM_CXXFLAGS])
AC_SUBST([AM_LDFLAGS])
AC_SUBST([LIBS])

# files to generate via autotools (.am or .in source files)
AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([src/Makefile])

# generate the final Makefile etc.
AC_OUTPUT
