IDIR =
CC=gcc
CXX=g++
CFLAGS=-I$(IDIR)

ODIR=obj
LDIR = 

LIBS=

_DEPS = psort.h
DEPS  = $(patsubstr %,$(IDIR)/%,$(_DEPS))

_OBJ = psort.o 
OBJ  = $(patsubstr %,$(ODIR)/%,$(_OBJ))

$(ODIR)/%.o: %.cpp $(DEPS)
	$(CXX) -c -o $@ $< $(CFLAGS)

all: $(OBJ)
	$(CXX) -o $@

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~

