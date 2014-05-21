IDIR = /opt/NAG/clmi623dgl/include
CXX = clang++
CFLAGS = -I$(IDIR)
LFLAGS = -L$(LDIR)

ODIR = obj
LDIR = /opt/NAG/clmi623dgl/lib

LIBS = -lnagc_nag

_DEPS = 
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = main_nag.o 
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))


$(ODIR)/%.o: %.cpp $(DEPS)
	$(CXX) -std=c++11 -c -g -o $@ $< $(CFLAGS)

maintarget: $(OBJ)
	$(CXX) -o $@ $^ $(CFLAGS) $(LFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o
