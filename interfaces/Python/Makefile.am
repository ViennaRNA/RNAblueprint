# this is: interfaces/Python/Makefile.am
if ENABLE_PYTHON

pyexec_LTLIBRARIES = _RNAblueprint.la

_RNAblueprint_la_SOURCES = $(top_srcdir)/interfaces/RNAblueprint.i
nodist__RNAblueprint_la_SOURCES = RNAblueprint_wrap.cxx $(SWIG_SOURCES)
nodist_python_PYTHON = RNAblueprint.py

_RNAblueprint_la_CPPFLAGS = -I$(top_srcdir)/lib $(PYTHON_CPPFLAGS) $(SWIG_PYTHON_CPPFLAGS) $(BOOST_CPPFLAGS)
_RNAblueprint_la_LIBADD = $(LDADD) $(top_srcdir)/lib/libRNAblueprint.la $(PYTHON_LIBS) $(PYTHON_EXTRA_LIBS) $(BOOST_GRAPH_LIBS)
_RNAblueprint_la_LDFLAGS = -avoid-version -module -shared -export-dynamic $(PYTHON_EXTRA_LDFLAGS) $(BOOST_GRAPH_LDFLAGS)

RNAblueprint_wrap.cxx: $(top_srcdir)/interfaces/RNAblueprint.i
	$(SWIG) -I${top_srcdir}/lib -includeall -o RNAblueprint_wrap.cxx -python ${top_srcdir}/interfaces/RNAblueprint.i

RNAblueprint.py: RNAblueprint_wrap.cxx

clean-local:
	rm -rf RNAblueprint_wrap.cxx RNAblueprint.py RNAblueprint.pyc

endif
