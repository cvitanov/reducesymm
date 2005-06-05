# Makefile to produce a static library and modules related to integration of ODE's.
# by vasimos

# Choose compiler:
COMP = pgf90 
OPTFULL = -fastsse
TARGET = libintegr
OBJ = eulerJ.o integrP.o rk2J.o rk2Jdriver.o rk4.o rk4driver.o rk4P.o rk4Pdriver.o 
targetMod = ifc_integr.mod
# Where to find user libraries:
userLibPATH = /home/vasimos/lib64
# Where to find user modules:
userModPATH = $(userLibPATH)/modules
# Names of user libraries
nrLIB = nr
integrLIB= integr

$(TARGET): $(OBJ)
	ar cr $(TARGET).a $(OBJ)
	ranlib $(TARGET).a

%.o: %.f90
	$(COMP) -c $(OPTFULL) -o $@ -I$(userModPATH) -L$(userLibPATH)  $< -l$(nrLIB)

$(targetMod): $(OBJ)
	$(COMP) -c -I$(userModPATH) ifc_integr.f90 


install:
	cp *.mod $(userModPATH)
	cp $(TARGET).a $(userLibPATH)

uninstall:
	rm -f $(userModPATH)/$(targetMod)
	rm -f $(userLibPATH)/$(TARGET).a

clean:
	rm -f *.o *.mod *.a