# Makefile to produce a static library and modules related to integration of ODE's.
# by vasimos

# Choose compiler:
COMP = pgf90 
OPTFULL = -fastsse
TARGET = libnewton
OBJ = newtonPO.o setLHM.o  setRHS.o
srcMod = ifc_newton.f90 
targetMod = ifc_newton.mod
# Where to find user libraries:
userLibPATH = /home/vasimos/lib64
# Where to find user modules:
userModPATH = $(userLibPATH)/modules
# Names of user libraries
nrLIB = nr
luLIB = lu 
utilLIB = util

$(TARGET): $(OBJ) $(targetMod) uninstall
	ar cr $(TARGET).a $(OBJ)
	ranlib $(TARGET).a

%.o: %.f90 
	$(COMP) -c $(OPTFULL) -o $@ -I$(userModPATH) -L$(userLibPATH)  $< -l$(luLIB) -l$(utilLIB) -l$(nrLIB)

$(OBJ): $(targetMod)

$(targetMod): $(srcMod)
	$(COMP) -c $(OPTFULL) -I$(userModPATH) $(srcMod)  


install:
	cp *.mod $(userModPATH)
	cp $(TARGET).a $(userLibPATH)

uninstall:
	rm -f $(userModPATH)/$(targetMod)
	rm -f $(userLibPATH)/$(TARGET).a

clean:
	rm -f *.o *.mod *.a
