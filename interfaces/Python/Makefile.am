# this is: interfaces/Python/Makefile.am

SOURCES = setup.py

swig-python: ../RNAdesign.i
	$(SWIG) -python -outcurrentdir ../RNAdesign.i

RNAdesign_wrap.cxx: swig-python

RNAdesign.py: swig-python

python_interface: RNAdesign_wrap.cxx RNAdesign.py setup.py
	@PYTHON@ setup.py build

all-local: python_interface	

install-exec-local: all-local
	@PYTHON@ setup.py install --prefix=$(prefix) --record install.record

uninstall-local:
	cat install.record | xargs rm -rf
	
mostlyclean-local: 
	@PYTHON@ setup.py clean --all
	rm -rf RNAdesign.py RNAdesign_wrap.cxx
