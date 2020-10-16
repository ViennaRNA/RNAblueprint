swig_verbose = $(swig_verbose_@AM_V@)
swig_verbose_ = $(swig_verbose_@AM_DEFAULT_V@)
swig_verbose_0 = @echo "  SWIG     $@";

SWIG_main_src = $(srcdir)/../RNAblueprint.i

SWIG_module_name = RNAblueprint

SWIG_wrapper = RNAblueprint_wrap.cpp

SWIG_src = $(SWIG_main_src)
