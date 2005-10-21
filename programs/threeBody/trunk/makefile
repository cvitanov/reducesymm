# by vasimos

# Choose compiler:
COMP = pgf90 
OPTFULL = -fastsse
TARGET = libthreeBody
OBJ = CartToKSMcG_TBL.o ifc_threeBody.o KSMcGtoCart_TBL.o TBLint.o TBrkqs.o
srcMod = ifc_threeBody.f90 
targetMod = ifc_threeBody.mod
# Where to find user libraries:
userLibPATH = /home/vasimos/lib
# Where to find user modules:
userModPATH = $(userLibPATH)/modules
# Names of user libraries
nrLIB = nr
integrLIB= integr

$(TARGET): uninstall $(OBJ) $(targetMod) 
	ar cr $(TARGET).a $(OBJ)
	ranlib $(TARGET).a

%.o: %.f90 
	$(COMP) -c $(OPTFULL) -o $@ -I$(userModPATH)  $<

$(OBJ): $(targetMod)

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
