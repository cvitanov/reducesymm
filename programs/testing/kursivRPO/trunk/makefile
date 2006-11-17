# command line compiler options can be passed in CFLAGS

# Choose compiler
COMP=
OBJECTS = RPOshootingTest.o parameters.o SetANdiag_KS.o SetLin_KS.o  SetNlin_KS.o 
# Set the outputfile
TARGET = ks 	
# Where to find user libraries:
userLibPATH = /home/$(USER)/lib-$(COMP)
# Where to find user modules:
userModPATH = $(userLibPATH)/modules
# Where to find shared libraries: 
sharedLibPATH = /usr/lib64
# Where to find lapack:
lapackModPATH = /usr/pgi/linux86-64/5.2/lib/lapack95_modules
lapackLibPATH = /usr/pgi/linux86-64/5.2/lib/
paramMod = parameters.mod
nrLIB = nr
integrLIB = integr
shootingLIB = rpo 
utilLIB = util
lapackLIB = lapack95
acmlLIB = acml

$(TARGET): $(OBJECTS)
	$(COMP) -o $(TARGET) -I$(userModPATH) -I$(lapackModPATH)  -L$(lapackLibPATH) -L$(userLibPATH) $(CFLAGS) $(OBJECTS) $(sharedLibPATH)/libfftw3.so -l$(shootingLIB) -l$(lapackLIB) -l$(acmlLIB)  -l$(integrLIB) -l$(utilLIB) -l$(nrLIB) 

%.o: %.f90 $(paramMod)
	$(COMP) -c -o $@ -I$(userModPATH) -I$(lapackModPATH) $(CFLAGS) $< 

$(paramMod): parameters.f90
	$(COMP) -c -o parameters.o -I$(userModPATH) -I$(lapackModPATH) $(CFLAGS) parameters.f90 


	
clean:
	rm -f *.o $(TARGET) 
	
cleanall:
	rm -f *.o *.mod $(TARGET) 

cleandat: 
	rm -f *.dat	

cleanedit: 
	rm -f *~


	
	
