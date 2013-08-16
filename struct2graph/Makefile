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
CSRC = common.cc printgraph.cc decompose.cc treestatistics.cc parsestruct.cc pathcoloring.cc graphcoloring.cc main.cc
COBJ = $(CSRC:.cc=.o)
# rules
all: $(EXE)
	@echo RNAdesign has been compiled
$(EXE): $(COBJ)
	$(CC) $(CCFLAGS) $(LDFLAGS) -o $(EXE) $(COBJ) 
.cc.o:
	$(CC) -c $(CCFLAGS) $< -o $@
clean:
	rm -f $(COBJ)
	rm -f $(EXE)
	rm -f *~ */*~
PHONY: clean
