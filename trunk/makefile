# Makefile to produce a static library and modules from nr*.f90 modules.
# by vasimos

COMP = pgf90
OPTFULL = -fastsse
TARGET = libnr
OBJ = nrtype.o nr.o nrutil.o 
MOD = nrtype.mod nr.mod nrutil.mod
# Where to find user libraries:
userLibPATH = /home/vasimos/lib64
# Where to find user modules:
userModPATH = $(userLibPATH)/modules

$(TARGET): $(OBJ)
	ar cr $(TARGET).a $(OBJ)
	ranlib $(TARGET).a

%.o: %.f90
	$(COMP) -c $(OPTFULL) -o $@  $< 

#$(MOD): $(OBJ)

$(OBJ): nrtype.f90


install:
	cp *.mod $(userModPATH)
	cp $(TARGET).a $(userLibPATH)
clean:
	rm -f *.o *.mod *.a