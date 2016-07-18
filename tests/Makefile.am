# this is: tests/Makefile.am
if ENABLE_TESTS

# the program to build (the names of the final binaries)
check_PROGRAMS = testRNAblueprint

# additional include pathes necessary to compile the C++ programs
testRNAblueprint_CPPFLAGS = -I$(top_srcdir)/lib $(BOOST_CPPFLAGS) @AM_CXXFLAGS@
testRNAblueprint_LDFLAGS = $(BOOST_GRAPH_LDFLAGS) $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
testRNAblueprint_LDADD = $(LDADD) $(top_builddir)/lib/libRNAblueprint.la $(BOOST_GRAPH_LIBS) $(BOOST_UNIT_TEST_FRAMEWORK_LIBS) $(AM_LDFLAGS)

# list of sources for the 'testRNAblueprint' binary
testRNAblueprint_SOURCES =	main.cc \
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
			test_dependency_graph.cc \
			test_RNAblueprint.cc
noinst_HEADERS = test_common.h

TESTS = testRNAblueprint

endif
