# this is: interfaces/Python/Makefile.am
if ENABLE_PYTHON

pyexec_LTLIBRARIES = _RNAdesign.la

_RNAdesign_la_SOURCES = $(top_srcdir)/interfaces/RNAdesign.i
nodist__RNAdesign_la_SOURCES = RNAdesign_wrap.cxx $(SWIG_SOURCES)
nodist_python_PYTHON = RNAdesign.py

_RNAdesign_la_CPPFLAGS = -I$(top_srcdir)/lib $(PYTHON_CPPFLAGS) $(SWIG_PYTHON_CPPFLAGS) $(BOOST_CPPFLAGS)
_RNAdesign_la_LIBADD = $(LDADD) $(top_srcdir)/lib/libRNAdesign.la $(PYTHON_LIBS) $(PYTHON_EXTRA_LIBS) $(BOOST_GRAPH_LIBS)
_RNAdesign_la_LDFLAGS = -avoid-version -module -shared -export-dynamic $(PYTHON_EXTRA_LDFLAGS) $(BOOST_GRAPH_LDFLAGS)

RNAdesign_wrap.cxx: $(top_srcdir)/interfaces/RNAdesign.i
	$(SWIG) -I${top_srcdir}/lib -includeall -o RNAdesign_wrap.cxx -python ${top_srcdir}/interfaces/RNAdesign.i

RNAdesign.py: RNAdesign_wrap.cxx

clean-local:
	rm -rf RNAdesign_wrap.cxx RNAdesign.py RNAdesign.pyc

endif
