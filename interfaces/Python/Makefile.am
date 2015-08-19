# this is: interfaces/Python/Makefile.am

SOURCES = setup.py

RNAdesign_wrap.cxx: ../RNAdesign.i
	$(SWIG) -python -outcurrentdir ../RNAdesign.i

python_interface: RNAdesign_wrap.cxx setup.py
	@PYTHON@ setup.py build

all-local: python_interface	

install-exec-local: all-local
	@PYTHON@ setup.py install --prefix=$(prefix) --record install.record

mostlyclean-local: 
	@PYTHON@ setup.py clean --all
	rm -rf RNAdesign.py RNAdesign_wrap.cxx
