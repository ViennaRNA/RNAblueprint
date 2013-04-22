# vars
CC = g++
LL = g++
EXE =  struct2graph
LIBS = libstdc++
DEFS =
WARN = -Wall
OPTI = -g
STD  = -std=c++11
CCFLAGS = $(WARN) $(STD) $(OPTI)
CSRC = struct2graph.cc
COBJ = $(foreach file, $(CSRC:.cc=.o), $(file))
# rules
all: $(EXE)
$(COBJ): %.o: %.cc
	$(CC) -c $(CCFLAGS) $< -o $@
clean:
	rm -f $(COBJ)
	rm -f $(EXE) *~
PHONY: clean
