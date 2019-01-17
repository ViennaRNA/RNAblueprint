# this is: src/Makefile.am
if ENABLE_PROGRAM

# the program to build (the names of the final binaries)
bin_PROGRAMS = RNAblueprint

# additional include pathes necessary to compile the C++ programs
RNAblueprint_CPPFLAGS = -I$(top_srcdir)/lib $(BOOST_CPPFLAGS) @AM_CXXFLAGS@
RNAblueprint_LDFLAGS =$(LDFLAGS) $(BOOST_GRAPH_LDFLAGS) $(BOOST_PROGRAM_OPTIONS_LDFLAGS)
RNAblueprint_LDADD = $(LDADD) $(top_builddir)/lib/libRNAblueprint.la $(BOOST_GRAPH_LIBS) $(BOOST_PROGRAM_OPTIONS_LIBS) $(AM_LDFLAGS)

# list of sources for the 'RNAblueprint' binary
RNAblueprint_SOURCES =	main.cc
			
noinst_HEADERS = main.h

endif
