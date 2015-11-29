# this is: lib/Makefile.am

# the library names to build
lib_LTLIBRARIES = libRNAdesign.la

# where to install the headers on the system
libRNAdesign_ladir = $(includedir)/RNAdesign

# ldflags to avoid version information
libRNAdesign_la_LDFLAGS = -avoid-version

# additional include pathes necessary to compile the C++ library
AM_CPPFLAGS = $(BOOST_CPPFLAGS) @AM_CXXFLAGS@

# the list of header files that belong to the library
libRNAdesign_la_HEADERS = \
            RNAdesign.h \
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
libRNAdesign_la_SOURCES =	\
			$(libRNAdesign_a_HEADERS) \
			RNAdesign.cc \
			common.cc \
			decompose.cc \
			printgraph.cc \
			parsestruct.cc \
			pathcoloring.cc \
			pairing_matrix.cc \
			dependency_graph.cc \
			probability_matrix.cc
