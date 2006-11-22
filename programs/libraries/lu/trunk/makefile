# Makefile to produce a static library and modules related to integration of ODE's.
# by $(USER)

# Choose compiler:
TARGET = liblu
OBJ = lubksb.o ludcmp.o
srcMod = ifc_lu.f90 
targetMod = ifc_lu.mod
# Where to find user libraries:
userLibPATH = /home/$(USER)/lib-$(COMP)
# Where to find user modules:
userModPATH = $(userLibPATH)/modules
# Names of user libraries
nrLIB = nr
luLIB= lu 

$(TARGET): $(OBJ) $(targetMod) uninstall
	ar cr $(TARGET).a $(OBJ)
	ranlib $(TARGET).a

%.o: %.f90 
	$(COMP) -c $(CFLAGS) -o $@ -I$(userModPATH) -L$(userLibPATH)  $< -l$(nrLIB)

$(OBJ): $(targetMod)

$(targetMod): $(srcMod)
	$(COMP) -c $(CFLAGS) -I$(userModPATH) $(srcMod)  


install:
	cp *.mod $(userModPATH)
	cp $(TARGET).a $(userLibPATH)

uninstall:
	rm -f $(userModPATH)/$(targetMod)
	rm -f $(userLibPATH)/$(TARGET).a

clean:
	rm -f *.o *.mod *.a
