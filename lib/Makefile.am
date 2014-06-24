# this is: lib/Makefile.am

# the library names to build
lib_LIBRARIES = libRNAdesign.a

# where to install the headers on the system
libRNAdesign_adir = $(includedir)/RNAdesign

# additional include pathes necessary to compile the C++ library
AM_CPPFLAGS = $(BOOST_CPPFLAGS) @AM_CXXFLAGS@

# the list of header files that belong to the library
libRNAdesign_a_HEADERS = \
      common.h \
			probability_matrix.h \
			graphcoloring.h \
			pathcoloring.h \
			decompose.h \
			graphcommon.h \
			parsestruct.h \
			printgraph.h \
			dependency_graph.h

# the sources to add to the library and to add to the source distribution
libRNAdesign_a_SOURCES =	\
      $(libRNAdesign_a_HEADERS) \
			common.cc \
			graphcommon.cc \
			printgraph.cc \
			decompose.cc \
			parsestruct.cc \
			pathcoloring.cc \
			probability_matrix.cc \
			graphcoloring.cc \
			dependency_graph.cc
