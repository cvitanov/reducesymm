# Makefile to produce a static library and modules related to integration of ODE's.
# by $(USER)

# Choose compiler:
COMP = 
TARGET = librpo
OBJ = ifc_rpo.o newtonRPOetdrk4.o setLHM.o  setRHS.o
srcMod = ifc_rpo.f90 
targetMod = ifc_rpo.mod
# Where to find user libraries:
userLibPATH = /home/$(USER)/lib-$(COMP)
# Where to find user modules:
userModPATH = $(userLibPATH)/modules
lapackModPATH = /usr/lib64/lapack95_modules
lapackLibPATH = /usr/lib64/
# Names of user libraries
nrLIB = nr
utilLIB = util
lapackLIB = lapack95
acmlLIB = acml

$(TARGET): uninstall $(OBJ) $(targetMod) 
	ar cr $(TARGET).a $(OBJ)
	ranlib $(TARGET).a

%.o: %.f90 
	$(COMP) -c $(CFLAGS) -o $@ -I$(lapackModPATH)  -L$(lapack95path) -L$(lapackLibPATH) -I$(userModPATH) -L$(userLibPATH)  $< -l$(lapackLIB) -l$(acmlLIB) -l$(utilLIB) -l$(nrLIB)

$(OBJ): $(targetMod)

$(targetMod): $(srcMod)
	$(COMP) -c $(CFLAGS) -I$(lapackModPATH) -I$(userModPATH) $(srcMod) -l$(lapackLIB) 
 


install:
	cp *.mod $(userModPATH)
	cp $(TARGET).a $(userLibPATH)

uninstall:
	rm -f $(userModPATH)/$(targetMod)
	rm -f $(userLibPATH)/$(TARGET).a

clean:
	rm -f *.o *.mod *.a
