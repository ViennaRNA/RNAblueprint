# this is: interfaces/Perl/Makefile.am

SOURCES = Makefile.PL

EXTRA_DIST = RNAdesign.pm RNAdesign_wrap.cxx

Makefile.perl: Makefile.PL Makefile
	@PERL@ Makefile.PL PREFIX=$(prefix) INSTALL_BASE= 

RNAdesign_wrap.cxx: ../RNAdesign.i
	$(SWIG) -perl -outcurrentdir ../RNAdesign.i

all-local: RNAdesign_wrap.cxx Makefile.perl
	$(MAKE) -f Makefile.perl

install-exec-local: all-local
	$(MAKE) -f Makefile.perl install

uninstall-local:
	$(MAKE) -f Makefile.perl uninstall-custom

mostlyclean-local:
	$(MAKE) -f Makefile.perl clean
	rm -f Makefile.perl Makefile.perl.old
	
# remove all generated content
maintainer-clean-local:
	rm -f RNAdesign.pm RNAdesign_wrap.cxx 

