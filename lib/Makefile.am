# this is: lib/Makefile.am

# the library names to build
lib_LIBRARIES = libRNAdesign.a

# where to install the headers on the system
libRNAdesign_adir = $(includedir)/RNAdesign

# additional include pathes necessary to compile the C++ library
AM_CPPFLAGS = $(BOOST_CPPFLAGS) @AM_CXXFLAGS@

# the list of header files that belong to the library
libRNAdesign_a_HEADERS = RNAdesign.h

noinst_HEADERS = \
			common.h \
			decompose.h \
			printgraph.h \
			graphcommon.h \
			parsestruct.h \
			pathcoloring.h \
			pairing_matrix.h \
			dependency_graph.h \
			probability_matrix.h \
			ear-decomposition.hpp

# the sources to add to the library and to add to the source distribution
libRNAdesign_a_SOURCES =	\
			$(libRNAdesign_a_HEADERS) \
			RNAdesign.cc \
			common.cc \
			decompose.cc \
			printgraph.cc \
			graphcommon.cc \
			parsestruct.cc \
			pathcoloring.cc \
			pairing_matrix.cc \
			dependency_graph.cc \
			probability_matrix.cc
