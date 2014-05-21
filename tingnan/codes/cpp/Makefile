IDIR = /opt/NAG/cll6a23dgl/include 
# /opt/NAG/clmi623dgl/include 
CXX = clang++
CFLAGS = -I$(IDIR)
LFLAGS = -L$(LDIR)

ODIR = obj
LDIR = /opt/NAG/cll6a23dgl/lib 
# /opt/NAG/clmi623dgl/lib

LIBS = -lnagc_nag -lpthread

_DEPS = 
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = main_nag.o 
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))


$(ODIR)/%.o: %.cpp $(DEPS)
	$(CXX) -std=c++11 -c -O2 -o $@ $< $(CFLAGS)

maintarget: $(OBJ)
	$(CXX) -o $@ $^ $(CFLAGS) $(LFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o
