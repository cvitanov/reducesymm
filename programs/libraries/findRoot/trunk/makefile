# Choose compiler:
OPTFULL = -O2
TARGET = libfindRoot
OBJ = ifc_findRoot.o rtbis.o zbrac.o zbracB.o
srcMod = ifc_findRoot.f90 
targetMod = ifc_findRoot.mod
# Where to find user libraries:
userLibPATH = /home/$(USER)/lib-$(COMP)
# Where to find user modules:
userModPATH = $(userLibPATH)/modules
# Names of user libraries
nrLIB = nr

$(TARGET): uninstall $(OBJ) $(srcMod) 
	make $(targetMod)
	ar cr $(TARGET).a $(OBJ)
	ranlib $(TARGET).a

%.o: %.f90 
	$(COMP) -c $(OPTFULL) -o $@ -I$(userModPATH)  $<

$(OBJ): $(srcMod)

$(targetMod): $(srcMod)
	$(COMP) -c $(OPTFULL) -I$(userModPATH) $(srcMod)  


install:
	cp *.mod $(userModPATH)
	cp $(TARGET).a $(userLibPATH)

uninstall:
	cd $(userModPATH)/;\
	rm -f $(targetMod)
	rm -f $(userLibPATH)/$(TARGET).a

clean:
	rm -f *.o *.mod *.a
