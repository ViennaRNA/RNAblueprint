# this is: interfaces/Perl/Makefile.am

SOURCES = Makefile.PL

RNAdesign_wrap.cxx: ../RNAdesign.i
	$(SWIG) -perl -outcurrentdir ../RNAdesign.i

Makefile.perl: Makefile.PL Makefile
	@PERL@ Makefile.PL PREFIX=$(prefix) INSTALL_BASE= 

perl_interface: RNAdesign_wrap.cxx Makefile.perl
	$(MAKE) -f Makefile.perl

all-local: perl_interface	

install-exec-local: all-local
	$(MAKE) -f Makefile.perl install

uninstall-local:
	$(MAKE) -f Makefile.perl clean

mostlyclean-local: 
	$(MAKE) -f Makefile.perl clean
	rm -rf RNAdesign.pm RNAdesign_wrap.cxx Makefile.perl Makefile.perl.old
