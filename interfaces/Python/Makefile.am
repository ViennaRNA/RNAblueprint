# this is: interfaces/Python/Makefile.am
if ENABLE_PYTHON

pkgpyexec_LTLIBRARIES = _RNAdesign.la
_RNAdesign_la_CPPFLAGS = -I$(top_srcdir)/lib $(PYTHON_CPPFLAGS) $(SWIG_PYTHON_CPPFLAGS)
# link to the library
_RNAdesign_la_LIBADD = ../lib/RNAdesign.la
# the source code for our extensions module
nodist__RNAdesign_la_SOURCES = RNAdesign_wrap.cxx $(SWIG_SOURCES)
# -module forces libtool to generate a dynamically loadable module
_RNAdesign_la_LDFLAGS = -avoid-version -module -shared -export-dynamic

nodist_pkgpython_PYTHON = RNAdesign.py

# this line assumes that the swig interface file RNAdesign.i is located in the interfaces subdirectory
RNAdesign_wrap.cxx: $(top_srcdir)/interfaces/RNAdesign.i
  $(SWIG) -I${top_srcdir}/lib -includeall -o RNAdesign_wrap.cxx -python ${top_srcdir}/interfaces/RNAdesign.i

RNAdesign.py: RNAdesign_wrap.cxx

clean-local:
  rm -rf RNAdesign_wrap.cxx RNAdesign.py RNAdesign.pyc

endif
