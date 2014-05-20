IDIR = /usr/local/include
CXX = clang++
CFLAGS = -I$(IDIR)

ODIR = obj
LDIR = /usr/local/lib

LIBS = -lgsl -lgslcblas

_DEPS = 
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = main.o 
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))


$(ODIR)/%.o: %.cpp $(DEPS)
	$(CXX) -std=c++11 -c -o $@ $< $(CFLAGS)

maintarget: $(OBJ)
	$(CXX) -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o
