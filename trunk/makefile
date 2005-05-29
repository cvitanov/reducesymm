# Makefile to produce a static library and modules from nr*.f90 modules.

COMP = pgf90
OPTFULL = -fastsse
TARGET = libnr
OBJ = nrtype.o nr.o nrutil.o 
MOD = nrtype.mod nr.mod nrutil.mod
INSTALLDIR = /home/vasimos/lib64/nr

$(TARGET): $(OBJ)
	ar cr $(TARGET).a $(OBJ)
	ranlib $(TARGET).a

%.o: %.f90
	$(COMP) -c $(OPTFULL) -o $@  $< 

#$(MOD): $(OBJ)

$(OBJ): nrtype.f90


install:	
	cp *.mod $(INSTALLDIR)
	cp $(TARGET).a $(INSTALLDIR)
clean:
	rm -f *.o *.mod *.a