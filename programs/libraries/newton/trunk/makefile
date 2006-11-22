##
TARGET = libnewt
OBJ = ifc_newt.o mnewt.o mnewtPO.o mnewtRPO.o mnewtRPOdamp.o mnewtTW.o mnewt_c.o 
srcMod = ifc_newt.f90 
targetMod = ifc_newt.mod
# Where to find user libraries:
userLibPATH = /home/$(USER)/lib-$(COMP)
# Where to find user modules:
userModPATH = $(userLibPATH)/modules
lapackModPATH = $(lapack95path)/$(COMP)/lapack95_modules
lapackLibPATH = /usr/lib$(local_arch)/
# Names of user libraries
nrLIB = nr
lapack95 = lapack95-$(COMP)
lapack = lapack77
blast = blast
acmlLIB = acml

$(TARGET): uninstall $(OBJ) $(targetMod)
	ar cr $(TARGET).a $(OBJ)
	ranlib $(TARGET).a

%.o: %.f90
	$(COMP) -c $(CFLAGS) -o $@ -I$(lapackModPATH)  -L$(lapack95path)/$(COMP) -L$(lapackLibPATH) -I$(userModPATH) -L$(userLibPATH)  $< -l$(lapack95) -l$(lapack) -l$(blast) -l$(nrLIB)

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
