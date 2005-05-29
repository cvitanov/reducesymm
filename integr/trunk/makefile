# Makefile to produce a static library and modules related to integration of ODE's.
# by vasimos

# Choose compiler:
COMP = pgf90 
OPTFULL = -fastsse
TARGET = libintegr
OBJ = rk4.o rk4driver.o
MOD = ifc_integr.mod
# Where to find user libraries:
userLibPATH = /home/vasimos/lib64
# Where to find user modules:
userModPATH = $(userLibPATH)/modules
nrLIB = nr

$(TARGET): $(OBJ)
	ar cr $(TARGET).a $(OBJ)
	ranlib $(TARGET).a

%.o: %.f90
	$(COMP) -c $(OPTFULL) -o $@ -I$(userModPATH) -L$(userLibPATH)  $< -l$(nrLIB)

$(OBJ): $(MOD)

$(MOD): 
	$(COMP) -c -I$(userModPATH) -L$(userLibPATH) ifc_integr.f90 -l$(nrLIB)


install:
	cp *.mod $(userModPATH)
	cp $(TARGET).a $(userLibPATH)
clean:
	rm -f *.o *.mod *.a