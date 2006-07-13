arch=64
TARGET = libnewt
OBJ = ifc_newt.o mnewt.o mnewtPO.o mnewtRPO.o mnewtRPOdamp.o mnewtTW.o mnewt_c.o 
srcMod = ifc_newt.f90 
targetMod = ifc_newt.mod
# Where to find user libraries:
userLibPATH = /home/vasimos/lib-$(COMP)
# Where to find user modules:
userModPATH = $(userLibPATH)/modules
lapackModPATH = /usr/lib$(arch)/lapack95/$(COMP)/lapack95_modules
lapackLibPATH = /usr/lib$(arch)/
# Names of user libraries
nrLIB = nr
lapack95 = lapack95-$(COMP)
lapack77 = lapack
blast = blast
acmlLIB = acml

$(TARGET): uninstall $(OBJ) $(targetMod)
	ar cr $(TARGET).a $(OBJ)
	ranlib $(TARGET).a

%.o: %.f90
	$(COMP) -c $(CFLAGS) -o $@ -I$(lapackModPATH)  -L$(lapackLibPATH) -I$(userModPATH) -L$(userLibPATH)  $< -l$(lapack95) -l$(lapack77) -l$(blast) -l$(nrLIB)

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
