# Makefile to produce a static library and modules related to integration of ODE's.
# by $(USER)

# Choose compiler:
TARGET = libpo
OBJ = ifc_po.o newtonPOetdrk4.o setLHM.o  setRHS.o
srcMod = ifc_po.f90 
targetMod = ifc_po.mod
# Where to find user libraries:
userLibPATH = /home/$(USER)/lib-$(COMP)
# Where to find user modules:
userModPATH = $(userLibPATH)/modules
lapackModPATH = /usr/lib/lapack95_modules-$(COMP)
lapackLibPATH = /usr/lib/
# Names of user libraries
nrLIB = nr
utilLIB = util
lapack95 = lapack95-$(COMP)
lapack77 =lapack 
blast = blast
acmlLIB = acml

$(TARGET): uninstall $(OBJ) $(targetMod)
	ar cr $(TARGET).a $(OBJ)
	ranlib $(TARGET).a

%.o: %.f90
	$(COMP) -c $(CFLAGS) -o $@ -I$(lapackModPATH)  -L$(lapack95path)/$(COMP) -L$(lapackLibPATH) -I$(userModPATH) -L$(userLibPATH)  $< -l$(lapack95) -l$(lapack77) -l-l$(blast) -l$(utilLIB) -l$(nrLIB)

$(OBJ): $(targetMod)

$(targetMod): $(srcMod)
	$(COMP) -c $(CFLAGS) -I$(lapackModPATH) -I$(userModPATH) $(srcMod)

install:
	cp *.mod $(userModPATH)
	cp $(TARGET).a $(userLibPATH)

uninstall:
	rm -f $(userModPATH)/$(targetMod)
	rm -f $(userLibPATH)/$(TARGET).a

clean:
	rm -f *.o *.mod *.a
