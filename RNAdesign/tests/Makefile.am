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
			test_pathcoloring.cc
nodist_include_HEADERS = test_common.h

testRNAdesign_LDADD = $(LDADD) -lboost_unit_test_framework
testRNAdesign_CPPFLAGS = -I$(top_srcdir)/src/

TESTS = testRNAdesign
