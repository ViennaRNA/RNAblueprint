# this is: interfaces/Python/Makefile.am

SOURCES = setup.py

EXTRA_DIST = RNAdesign_wrap.cxx setup.py

swig-python: ../RNAdesign.i
	$(SWIG) -python -outcurrentdir ../RNAdesign.i

RNAdesign_wrap.cxx: swig-python

#RNAdesign.py: swig-python
#RNAdesign.py 

python_interface: RNAdesign_wrap.cxx setup.py
	@PYTHON@ setup.py build

all-local: python_interface	

install-exec-local: all-local
	@PYTHON@ setup.py install --prefix=$(prefix) --record install.record

uninstall-local:
	cat install.record | xargs rm -rf
	
mostlyclean-local: 
	@PYTHON@ setup.py clean --all
	
# remove all generated content
maintainer-clean-local:
	rm -f RNAdesign.py RNAdesign_wrap.cxx
