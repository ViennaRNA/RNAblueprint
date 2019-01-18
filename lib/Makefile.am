# this is: lib/Makefile.am

# the library names to build
lib_LTLIBRARIES = libRNAblueprint.la

# where to install the headers on the system
libRNAblueprint_ladir = $(includedir)/RNAblueprint

# ldflags to avoid version information
libRNAblueprint_la_LDFLAGS = -avoid-version

# additional include pathes necessary to compile the C++ library
AM_CPPFLAGS = $(BOOST_CPPFLAGS) @AM_CXXFLAGS@

# the list of header files that belong to the library
libRNAblueprint_la_HEADERS = \
            RNAblueprint.h \
			common.h \
			decompose.h \
			printgraph.h \
			graphcommon.h \
			parsestruct.h \
			pathcoloring.h \
			pairing_matrix.h \
			dependency_graph.h \
			probability_matrix.h \
			uninduced_subgraph.hpp \
			ear-decomposition.hpp

# the sources to add to the library and to add to the source distribution
libRNAblueprint_la_SOURCES =	\
			$(libRNAblueprint_a_HEADERS) \
			RNAblueprint.cc \
			common.cc \
			decompose.cc \
			printgraph.cc \
			parsestruct.cc \
			pathcoloring.cc \
			pairing_matrix.cc \
			dependency_graph.cc \
			probability_matrix.cc
