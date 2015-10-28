# this is: interfaces/Perl/Makefile.am

if ENABLE_PERL5

perlextdir = $(PERL_EXT_LIB)

perlext_LTLIBRARIES = RNAdesign.la

RNAdesign_la_CPPFLAGS = -I$(top_srcdir)/lib $(PERL_EXT_CPPFLAGS)
# link to the library
RNAdesign_la_LIBADD = ../lib/RNAdesign.la
# the source code for our extensions module

nodist_RNAdesign_la_SOURCES = RNAdesign_wrap.cxx
# -module forces libtool to generate a dynamically loadable module
RNAdesign_la_LDFLAGS = -avoid-version -module -shared -export-dynamic

nodist_perlext_SCRIPTS = RNAdesign.pm

# this line assumes that the swig interface file RNAdesign.i is located in the interfaces subdirectory
RNAdesign_wrap.cxx: $(top_srcdir)/interfaces/RNAdesign.i
  $(SWIG) -I${top_srcdir}/lib -includeall -o RNAdesign_wrap.cxx -perl ${top_srcdir}/interfaces/RNAdesign.i

RNAdesign.pm: RNAdesign_wrap.cxx

clean-local:
  rm -rf RNAdesign_wrap.cxx RNAdesign.pm

endif
