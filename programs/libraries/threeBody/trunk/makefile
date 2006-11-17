OPTFULL = -O2
TARGET = libthreeBody
OBJ = CartToKSMcG_TBL.o findIC.o findSep.o ifc_threeBody.o KSMcGtoCart_TBL.o TBLint.o TBLint_sep.o TBrkqs.o vectorize.o
srcMod = ifc_threeBody.f90 
targetMod = ifc_threeBody.mod
# Where to find user libraries:
userLibPATH = /home/$(USER)/lib-$(COMP)
# Where to find user modules:
userModPATH = $(userLibPATH)/modules
# Names of user libraries
nrLIB = nr
integrLIB= integr

$(TARGET): uninstall $(OBJ) $(srcMod) 
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
	rm -f $(userModPATH)/$(targetMod)
	rm -f $(userLibPATH)/$(TARGET).a

clean:
	rm -f *.o *.mod *.a
