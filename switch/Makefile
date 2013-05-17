# vars
CC = g++
LL = g++
EXE =  switch
LIBS = libstdc++
LDFLAGS += -lboost_program_options
DEFS =
WARN = -Wall
OPTI = -g -lboost_program_options
STD  = -std=c++11
CCFLAGS = $(WARN) $(STD) $(OPTI)
CSRC = switch.cc
COBJ = $(foreach file, $(CSRC:.cc=.o), $(file))
# rules
all: $(EXE)
$(COBJ): %.o: %.cc
	$(CC) -c $(CCFLAGS) $< -o $@
clean:
	rm -f $(COBJ)
	rm -f $(EXE) *~ */*~
PHONY: clean
