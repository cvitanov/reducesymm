# command line compiler options can be passed in CFLAGS

# Choose compiler
COMP = 	
OBJECTS = etdrk4Test.o parameters.o SetLin_KS.o  SetNlin_KS.o
# Set the outputfile
TARGET = ks 	
# Where to find user libraries:
userLibPATH = /home/vasimos/lib-$(COMP)
# Where to find user modules:
userModPATH = $(userLibPATH)/modules
# Where to find shared libraries: 
sharedLibPATH = /usr/lib64
paramMod = parameters.mod
nrLIB = nr
integrLIB =integr


$(TARGET): $(OBJECTS)
	$(COMP) -o $(TARGET) -I$(userModPATH) -L$(userLibPATH) $(CFLAGS) $(OBJECTS) $(sharedLibPATH)/libfftw3.so -l$(integrLIB) -l$(nrLIB)

%.o: %.f90 $(paramMod)
	$(COMP) -c -o $@ -I$(userModPATH) $(CFLAGS) $< 

$(paramMod): parameters.f90
	$(COMP) -c -o parameters.o -I$(userModPATH) $(CFLAGS) parameters.f90 


	
clean:
	rm -f *.o $(TARGET) 
	
cleanall:
	rm -f *.o *.mod $(TARGET) 

cleandat: 
	rm -f *.dat	

cleanedit: 
	rm -f *~


	
	
