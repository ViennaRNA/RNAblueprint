# this is: interfaces/Python/Makefile.am

SOURCES = setup.py

EXTRA_DIST = RNAdesign_wrap.cxx RNAdesign.py

RNAdesign_wrap.cxx: ../RNAdesign.i
	$(SWIG) -python -outcurrentdir ../RNAdesign.i

all-local: RNAdesign_wrap.cxx setup.py
	@PYTHON@ setup.py build

install-exec-local: all-local
	@PYTHON@ setup.py install --prefix=$(prefix) --record install.record

uninstall-local:
	cat install.record | xargs rm -rf
	
mostlyclean-local: 
	@PYTHON@ setup.py clean --all
	
# remove all generated content
maintainer-clean-local:
	rm -f RNAdesign.py RNAdesign_wrap.cxx
