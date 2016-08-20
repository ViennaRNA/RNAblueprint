# this is: interfaces/Perl/Makefile.am
if ENABLE_PERL5
perlextdir = $(PERL_EXT_LIB)
perlext_LTLIBRARIES = RNAblueprint.la

RNAblueprint_la_SOURCES = $(top_srcdir)/interfaces/RNAblueprint.i
nodist_RNAblueprint_la_SOURCES = RNAblueprint_wrap.cxx $(SWIG_SOURCES)
nodist_perlext_SCRIPTS = RNAblueprint.pm

RNAblueprint_la_CPPFLAGS = -I$(top_srcdir)/lib $(PERL_EXT_CPPFLAGS) $(SWIG_PERL_CPPFLAGS) $(PERLXS_CFLAGS) $(BOOST_CPPFLAGS)
RNAblueprint_la_LIBADD = $(LDADD) $(top_srcdir)/lib/libRNAblueprint.la $(PERLXS_LIBS) $(BOOST_GRAPH_LIBS)
RNAblueprint_la_LDFLAGS = -avoid-version -module -shared -export-dynamic $(PERLXS_LDFLAGS) $(BOOST_GRAPH_LDFLAGS)

RNAblueprint_wrap.cxx: $(top_srcdir)/interfaces/RNAblueprint.i
	$(SWIG) -I${top_srcdir}/lib -includeall -o RNAblueprint_wrap.cxx -perl ${top_srcdir}/interfaces/RNAblueprint.i

RNAblueprint.pm: RNAblueprint_wrap.cxx

clean-local:
	rm -rf RNAblueprint_wrap.cxx RNAblueprint.pm

endif
