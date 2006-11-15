# Makefile to produce a static library and modules for rpo shooting in KSe.
# by vasimos

TARGET = libks
OBJ = ks.o Rc.o Rr.o SetLin_KS.o SetNlin_KS.o SetAndiag_KS.o SetAndiag_KS_diff.o
srcMod = ks.f90 
targetMod = ks.mod
# Where to find user libraries:
userLibPATH = /home/vasimos/lib-$(COMP)
# Where to find user modules:
userModPATH = $(userLibPATH)/modules
# Where to find shared libraries: 
sharedLibPATH = /usr/lib64
# Where to find lapack:
lapackLibPATH = /usr/lib64
lapackModPATH = $(lapackLibPATH)/lapack95_modules
# Names of user libraries
nrLIB = nr
integrLIB= integr

$(TARGET): uninstall $(OBJ) $(srcMod) 
	make $(targetMod)
	ar cr $(TARGET).a $(OBJ)
	ranlib $(TARGET).a

%.o: %.f90 
	$(COMP) -c $(CFLAGS) -o $@ -I$(userModPATH) -I$(lapackModPATH) $<

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