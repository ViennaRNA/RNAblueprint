include ../generic.mk

LANG_SPECIFIC_INTERFACE_FILES = \
  $(builddir)/version.i

INTERFACE_FILES = $(SWIG_src) \
                  $(LANG_SPECIFIC_INTERFACE_FILES)

pkgpyexecdir = $(py3execdir)/RNAblueprint
pkgpycachedir = $(pkgpyexecdir)/__pycache__

pkgpyexec_LTLIBRARIES = _RNAblueprint.la
pkgpyexec_DATA =  RNAblueprint/__init__.py
pkgpycache_DATA = RNAblueprint/__pycache__/__init__.@PYTHON3_CACHE_TAG@.pyc \
                  RNAblueprint/__pycache__/__init__.@PYTHON3_CACHE_OPT1_EXT@

_RNAblueprint_la_SOURCES = $(INTERFACE_FILES) \
                  $(SWIG_wrapper)

EXTRA_DIST = \
  $(SWIG_wrapper) \
  $(SWIG_module_name).py \
  $(LANG_SPECIFIC_INTERFACE_FILES)

_RNAblueprint_la_CPPFLAGS = -I$(PYTHON3_INC) -I$(top_srcdir)/lib $(BOOST_CPPFLAGS)
_RNAblueprint_la_CXXFLAGS = -Wno-unused-label -Wno-unused-variable

_RNAblueprint_la_LIBADD = \
    $(LDADD) \
    $(top_builddir)/lib/libRNAblueprint_conv.la \
    $(BOOST_GRAPH_LIBS) \
    -lstdc++

_RNAblueprint_la_LDFLAGS = \
                  -avoid-version \
                  -module \
                  -shared \
                  -export-dynamic \
                  -shrext \
                  $(BOOST_GRAPH_LDFLAGS) \
                  $(PYTHON3_SO) \
                  $(PYTHON3_LDFLAGS)

##
## SWIG interface wrapper code generation
##

if ENABLE_SWIG

$(SWIG_wrapper): $(SWIG_src) $(INTERFACE_FILES)
	$(swig_verbose)$(SWIG) -I$(top_srcdir)/src \
                  -I$(srcdir) -I$(builddir) \
                  -includeall \
                  -DPY3 -py3 \
                  -python -c++ -keyword -w511 \
                  -o $(SWIG_wrapper) \
                  $(SWIG_main_src)

$(SWIG_module_name).py: $(SWIG_wrapper)

endif

# We rename RNAblueprint.py to __init__.py so that it can be installed
# in an "RNAblueprint" directory and be treated as a package.
RNAblueprint/__init__.py: $(SWIG_module_name).py
	$(AM_V_GEN)( test -d RNAblueprint || $(MKDIR_P) RNAblueprint ) && \
  cp `test -f RNAblueprint.py || echo '$(srcdir)/'`RNAblueprint.py RNAblueprint/__init__.py

# We "import _RNAblueprint" first so that if we fail to import the glue library
# we don't generate a broken .pyc or .pyo.
RNAblueprint/__pycache__/__init__.@PYTHON3_CACHE_TAG@.pyc: RNAblueprint/__init__.py RNAblueprint/_RNAblueprint$(PYTHON3_SO)
	$(AM_V_GEN)( PYTHONPATH="RNAblueprint:$$PYTHONPATH" $(PYTHON3) -c "import _RNAblueprint" ; \
  PYTHONPATH=".:$$PYTHONPATH" $(PYTHON3) -c "import RNAblueprint" )

RNAblueprint/__pycache__/__init__.@PYTHON3_CACHE_OPT1_EXT@: RNAblueprint/__init__.py RNAblueprint/_RNAblueprint$(PYTHON3_SO)
	$(AM_V_GEN)( PYTHONPATH="RNAblueprint:$$PYTHONPATH" $(PYTHON3) -O -c "import _RNAblueprint" ; \
  PYTHONPATH=".:$$PYTHONPATH" $(PYTHON3) -O -c "import RNAblueprint" )

RNAblueprint/_RNAblueprint$(PYTHON3_SO): _RNAblueprint.la
	$(AM_V_GEN)( test -d RNAblueprint || $(MKDIR_P) RNAblueprint ) && \
  ( $(LIBTOOL) --config > libtoolconfig.tmp ) && \
  ( . ./libtoolconfig.tmp; cp $$objdir/_RNAblueprint$(PYTHON3_SO) RNAblueprint ) && \
  rm -f libtoolconfig.tmp

CLEANFILES = \
    RNAblueprint/_RNAblueprint$(PYTHON3_SO) \
    RNAblueprint/__init__.py \
    RNAblueprint/__pycache__/__init__.@PYTHON3_CACHE_TAG@.pyc \
    RNAblueprint/__pycache__/__init__.@PYTHON3_CACHE_OPT1_EXT@\
    $(srcdir)/*.pyc

all-local:  RNAblueprint/_RNAblueprint$(PYTHON3_SO)

clean-local:
	-rm -rf RNAblueprint

# Remove the .la file - RNAblueprint.la is never linked against (it's a module)
# and Python doesn't use libltdl.  Note that the library gets installed by
# install-data, so that's where we need to hook.
install-data-hook:
	rm -f $(DESTDIR)$(pkgpyexecdir)/_RNAblueprint.la

# Because we don't install the .la file, "make uninstall" doesn't work and
# we need to remove the file ourselves.
uninstall-local:
	eval `grep '^dlname=' $(pkgpyexec_LTLIBRARIES)` ; \
	rm -f $(DESTDIR)$(pkgpyexecdir)/"$$dlname"

## remove all generated content
maintainer-clean-local:
	-rm -f $(SWIG_wrapper) RNAblueprint.py

distclean-local:
	-rm -f $(SWIG_wrapper) RNAblueprint.py
