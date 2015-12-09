# this is: src/Makefile.am

# the program to build (the names of the final binaries)
bin_PROGRAMS = RNAdesign

# additional include pathes necessary to compile the C++ programs
RNAdesign_CPPFLAGS = -I$(top_srcdir)/lib $(BOOST_CPPFLAGS) @AM_CXXFLAGS@
RNAdesign_LDFLAGS = $(BOOST_GRAPH_LDFLAGS) $(BOOST_PROGRAM_OPTIONS_LDFLAGS) $(BOOST_REGEX_LDFLAGS)
RNAdesign_LDADD = $(LDADD) $(top_builddir)/lib/libRNAdesign.la $(BOOST_GRAPH_LIBS) $(BOOST_PROGRAM_OPTIONS_LIBS) $(BOOST_REGEX_LIBS) $(AM_LDFLAGS)

# list of sources for the 'RNAdesign' binary
RNAdesign_SOURCES =	main.cc
			
noinst_HEADERS = main.h
