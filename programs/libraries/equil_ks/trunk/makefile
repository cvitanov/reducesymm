# Makefile to produce a static library and modules related to integration of ODE's.
# by vasimos

TARGET = libequil_ks
OBJ = ifc_equil_ks.o ksFJ.o SetNlin_KS.o
srcMod = ifc_equil_ks.f90 
targetMod = ifc_equil_ks.mod
# Where to find user libraries:
userLibPATH = /home/vasimos/lib-$(COMP)
# Where to find user modules:
userModPATH = $(userLibPATH)/modules
# Names of user libraries
nrLIB = nr
integrLIB= integr

$(TARGET): uninstall $(OBJ) $(srcMod) 
	make $(targetMod)
	ar cr $(TARGET).a $(OBJ)
	ranlib $(TARGET).a

%.o: %.f90 
	$(COMP) -c $(CFLAGS) -o $@ -I$(userModPATH)  $<

$(OBJ): $(srcMod)

$(targetMod): $(srcMod)
	$(COMP) -c $(CFLAGS) -I$(userModPATH) $(srcMod)  


install:
	cp *.mod $(userModPATH)
	cp $(TARGET).a $(userLibPATH)

uninstall:
	cd $(userModPATH)/;\
	rm -f $(targetMod)
	rm -f $(userLibPATH)/$(TARGET).a

clean:
	rm -f *.o *.mod *.a
