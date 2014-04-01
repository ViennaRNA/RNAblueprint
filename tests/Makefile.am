# this is: tests/Makefile.am

# the program to build (the names of the final binaries)
check_PROGRAMS = testRNAdesign

# list of sources for the 'testRNAdesign' binary
testRNAdesign_SOURCES =	main.cc \
			test_decompose.cc \
			test_random_gen.cc \
			test_common.cc \
			test_parsestruct.cc \
			test_graphcoloring.cc \
			test_pathcoloring.cc \
			$(top_builddir)/src/common.h \
                        $(top_builddir)/src/probability_matrix.h \
                        $(top_builddir)/src/graphcoloring.h \
                        $(top_builddir)/src/pathcoloring.h \
                        $(top_builddir)/src/decompose.h \
                        $(top_builddir)/src/graphcommon.h \
                        $(top_builddir)/src/parsestruct.h \
                        $(top_builddir)/src/printgraph.h \
                        $(top_builddir)/src/dependency_graph.h

nodist_include_HEADERS = test_common.h

testRNAdesign_CPPFLAGS = $(BOOST_CPPFLAGS) -I$(top_srcdir)/src
testRNAdesign_LDFLAGS = $(BOOST_THREADS_LDFLAGS) $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
testRNAdesign_LDADD = $(LDADD) $(BOOST_THREAD_LIBS) $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS = testRNAdesign
