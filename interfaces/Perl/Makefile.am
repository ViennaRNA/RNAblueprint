# this is: interfaces/Perl/Makefile.am
if ENABLE_PERL5
perlextdir = $(PERL_EXT_LIB)
perlext_LTLIBRARIES = RNAdesign.la

RNAdesign_la_SOURCES = $(top_srcdir)/interfaces/RNAdesign.i
nodist_RNAdesign_la_SOURCES = RNAdesign_wrap.cxx $(SWIG_SOURCES)
nodist_perlext_SCRIPTS = RNAdesign.pm

RNAdesign_la_CPPFLAGS = -I$(top_srcdir)/lib $(PERL_EXT_CPPFLAGS) $(SWIG_PERL_CPPFLAGS) $(PERLXS_CFLAGS)
RNAdesign_la_LIBADD = $(LDADD) $(top_srcdir)/lib/libRNAdesign.la
RNAdesign_la_LDFLAGS = -avoid-version -module -shared -export-dynamic $(PERLXS_LDFLAGS)

RNAdesign_wrap.cxx: $(top_srcdir)/interfaces/RNAdesign.i
	$(SWIG) -I${top_srcdir}/lib -includeall -o RNAdesign_wrap.cxx -perl ${top_srcdir}/interfaces/RNAdesign.i

RNAdesign.pm: RNAdesign_wrap.cxx

clean-local:
	rm -rf RNAdesign_wrap.cxx RNAdesign.pm

endif
