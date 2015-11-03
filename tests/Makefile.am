# this is: tests/Makefile.am

# the program to build (the names of the final binaries)
check_PROGRAMS = testRNAdesign

# additional include pathes necessary to compile the C++ programs
testRNAdesign_CPPFLAGS = -I$(top_srcdir)/lib $(BOOST_CPPFLAGS) @AM_CXXFLAGS@
testRNAdesign_LDFLAGS = $(BOOST_THREADS_LDFLAGS) $(BOOST_PROGRAM_OPTIONS_LDFLAGS)
testRNAdesign_LDADD = $(LDADD) $(top_builddir)/lib/libRNAdesign.la $(BOOST_THREAD_LIBS) $(BOOST_PROGRAM_OPTIONS_LIBS) $(AM_LDFLAGS)

# list of sources for the 'testRNAdesign' binary
testRNAdesign_SOURCES =	main.cc \
			test_uninduced_subgraph.cc \
			test_subgraph.cc \
			test_pairing_matrix.cc \
			test_probability_matrix.cc \
			test_decompose.cc \
			test_random_gen.cc \
			test_common.cc \
			test_graphcommon.cc \
			test_parsestruct.cc \
			test_pathcoloring.cc \
			test_statistics.cc \
			test_dependency_graph.cc
noinst_HEADERS = test_common.h

TESTS = testRNAdesign
