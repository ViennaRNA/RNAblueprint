include ../generic.mk


PERL_MODULE = $(SWIG_module_name).pm

perlextdir = $(prefix)$(PERL_ARCH_RELATIVE_INSTALL_DIR)/auto/RNAblueprint
perlext_LTLIBRARIES = RNAblueprint.la

LANG_SPECIFIC_INTERFACE_FILES = \
  $(builddir)/version.i

INTERFACE_FILES = $(SWIG_src) \
                  $(LANG_SPECIFIC_INTERFACE_FILES)

RNAblueprint_la_SOURCES =  $(INTERFACE_FILES) \
                  $(SWIG_wrapper)

perlmoddir = $(prefix)$(PERL_LIB_RELATIVE_INSTALL_DIR)
nobase_dist_perlmod_DATA = \
                  $(PERL_MODULE)

RNAblueprint_la_CPPFLAGS = -I$(PERL_EXT_INC) -I$(top_srcdir)/lib $(BOOST_CPPFLAGS)
RNAblueprint_la_CXXFLAGS = -Wno-unused-label -Wno-unused-variable

RNAblueprint_la_LIBADD = \
                  $(LDADD) \
                  $(BOOST_GRAPH_LIBS) \
                  $(top_builddir)/lib/libRNAblueprint_conv.la \
                 -lstdc++

RNAblueprint_la_LDFLAGS = \
                  -avoid-version \
                  -module \
                  -shared \
                  -export-dynamic \
                  -shrext \
                  $(BOOST_GRAPH_LDFLAGS) \
                  $(PERL_SO)

EXTRA_DIST = \
  $(SWIG_wrapper) \
  $(PERL_MODULE) \
  $(LANG_SPECIFIC_INTERFACE_FILES)

##
## SWIG interface wrapper code generation
##

if ENABLE_SWIG

$(SWIG_wrapper) $(PERL_MODULE): $(INTERFACE_FILES)
	$(swig_verbose)$(SWIG) -I$(top_srcdir)/src \
                  -I$(srcdir) -I$(builddir) \
                  -includeall \
                  -perl5 -proxy -const -c++ \
                  -module RNAblueprint -o $(SWIG_wrapper) \
                  $(SWIG_main_src)

endif

##
## copy shared library for perl-tests
##

auto/RNAblueprint/RNAblueprint$(PERL_SO): RNAblueprint.la
	$(AM_V_GEN)( $(MKDIR_P) auto/RNAblueprint ) && \
  ( \
    $(LIBTOOL) --config > libtoolconfig.tmp && \
    ( . ./libtoolconfig.tmp; cp "$$objdir/RNAblueprint$(PERL_SO)" auto/RNAblueprint ) && \
    rm -f libtoolconfig.tmp \
  )

all-local:  auto/RNAblueprint/RNAblueprint$(PERL_SO)

clean-local:
	-rm -rf auto

## Remove the .la file - RNAblueprint.la is never linked against (it's a module)
## and Perl doesn't use libltdl.  Note that the library gets installed by
## install-data, so that's where we need to hook.
install-data-hook:
	rm -f $(DESTDIR)$(perlextdir)/RNAblueprint.la

## Because we don't install the .la file, "make uninstall" doesn't work and
## we need to remove the file ourselves.
uninstall-local:
	eval `grep '^dlname=' $(perlext_LTLIBRARIES)` ; \
	rm -f $(DESTDIR)$(perlextdir)/"$$dlname"

## remove all generated content
maintainer-clean-local:
	-rm -f $(SWIG_wrapper) $(PERL_MODULE)

distclean-local:
	-rm -f $(SWIG_wrapper) $(PERL_MODULE)
