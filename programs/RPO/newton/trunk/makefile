# Makefile to produce a static library and modules related to integration of ODE's.
# by vasimos

# Choose compiler:
COMP = pgf90 
OPTFULL = -fastsse
TARGET = librpo
OBJ = newtonRPOetdrk4.o setLHM.o  setRHS.o
srcMod = ifc_rpo.f90 
targetMod = ifc_rpo.mod
# Where to find user libraries:
userLibPATH = /home/vasimos/lib64
# Where to find user modules:
userModPATH = $(userLibPATH)/modules
lapackModPATH = /usr/pgi/linux86-64/5.2/lib/lapack95_modules
lapackLibPATH = /usr/pgi/linux86-64/5.2/lib/
# Names of user libraries
nrLIB = nr
utilLIB = util
lapackLIB = lapack95
acmlLIB = acml

$(TARGET): $(OBJ) $(targetMod) uninstall
	ar cr $(TARGET).a $(OBJ)
	ranlib $(TARGET).a

%.o: %.f90 
	$(COMP) -c $(OPTFULL) -o $@ -I$(lapackModPATH)  -L$(lapackLibPATH) -I$(userModPATH) -L$(userLibPATH)  $< -l$(lapackLIB) -l$(acmlLIB) -l$(utilLIB) -l$(nrLIB)

$(OBJ): $(targetMod)

$(targetMod): $(srcMod)
	$(COMP) -c $(OPTFULL) -I$(lapackModPATH) -I$(userModPATH) $(srcMod)  


install:
	cp *.mod $(userModPATH)
	cp $(TARGET).a $(userLibPATH)

uninstall:
	rm -f $(userModPATH)/$(targetMod)
	rm -f $(userLibPATH)/$(TARGET).a

clean:
	rm -f *.o *.mod *.a
