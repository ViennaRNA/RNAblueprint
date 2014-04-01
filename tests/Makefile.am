# this is: tests/Makefile.am

# the program to build (the names of the final binaries)
check_PROGRAMS = testRNAdesign

# list of sources for the 'RNAdesign' binary
testRNAdesign_SOURCES =	main.cc \
			test_decompose.cc \
			test_random_gen.cc \
			test_common.cc \
			test_parsestruct.cc \
			test_graphcoloring.cc \
			test_pathcoloring.cc \
			$(top_srcdir)/src/common.cc \
                        $(top_srcdir)/src/graphcommon.cc \
                        $(top_srcdir)/src/printgraph.cc \
                        $(top_srcdir)/src/decompose.cc \
                        $(top_srcdir)/src/parsestruct.cc \
                        $(top_srcdir)/src/pathcoloring.cc \
                        $(top_srcdir)/src/probability_matrix.cc \
                        $(top_srcdir)/src/graphcoloring.cc \
                        $(top_srcdir)/src/dependency_graph.cc

nodist_include_HEADERS = test_common.h \
			$(top_srcdir)/src/common.h \
                        $(top_srcdir)/src/probability_matrix.h \
                        $(top_srcdir)/src/graphcoloring.h \
                        $(top_srcdir)/src/pathcoloring.h \
                        $(top_srcdir)/src/decompose.h \
                        $(top_srcdir)/src/graphcommon.h \
                        $(top_srcdir)/src/parsestruct.h \
                        $(top_srcdir)/src/printgraph.h \
                        $(top_srcdir)/src/dependency_graph.h

testRNAdesign_CPPFLAGS = $(BOOST_CPPFLAGS) -I$(top_srcdir)/src
testRNAdesign_LDFLAGS = $(BOOST_THREADS_LDFLAGS) $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
testRNAdesign_LDADD = $(LDADD) $(BOOST_THREAD_LIBS) $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS = testRNAdesign
