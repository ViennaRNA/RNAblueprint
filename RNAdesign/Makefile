# vars
CC = g++
LL = g++
EXE = RNAdesign
LIBS = -lstdc++
LDFLAGS = -lboost_program_options
DEFS =
WARN = -Wall
OPTI = -g
STD  = -std=c++11
CCFLAGS = $(WARN) $(STD) $(OPTI) $(DEFS) $(LIBS)
CSRC = common.cc graphcommon.cc printgraph.cc decompose.cc treestatistics.cc parsestruct.cc pathcoloring.cc graphcoloring.cc
COBJ = $(CSRC:.cc=.o)
# sources only for program
PCSRC = main.cc
PCOBJ = $(PCSRC:.cc=.o)

# unit test related variables
TEXE = tests/testRNAdesign
TARGS = --log_level=test_suite
TLDFLAGS = -lboost_unit_test_framework
TCSRC = tests/main.cc tests/test_pathcoloring.cc tests/test_graphcoloring.cc tests/test_decompose.cc tests/test_parsestruct.cc
TCOBJ = $(TCSRC:.cc=.o)

# rules
all: $(EXE)
	@echo RNAdesign has been compiled
$(EXE): $(COBJ) $(PCOBJ)
	$(CC) $(CCFLAGS) $(LDFLAGS) -o $(EXE) $(COBJ) $(PCOBJ)
.cc.o:
	$(CC) -c $(CCFLAGS) $< -o $@
run:	$(EXE)
	./$(EXE)
$(TEXE): $(COBJ) $(TCOBJ)
	$(CC) $(CCFLAGS) $(TLDFLAGS) -o $(TEXE) $(COBJ) $(TCOBJ)
test:	$(TEXE)
	./$(TEXE) $(TARGS)
clean:
	rm -f $(COBJ) $(TCOBJ) $(PCOBJ)
	rm -f $(EXE) $(TEXE)
	rm -f *~ */*~
PHONY: clean
