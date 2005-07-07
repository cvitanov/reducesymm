# command line compiler options can be passed in CFLAGS

# Choose compiler
COMP = pgf90	
OBJECTS = integr.o  KSfield.o parameters.o
# Set the outputfile
TARGET = ks 	
# Where to find user libraries:
userLibPATH = /home/vasimos/lib64
# Where to find user modules:
userModPATH = $(userLibPATH)/modules
paramMod = parameters.mod
nrLIB = nr
integrLIB =integr


$(TARGET): $(OBJECTS)
	$(COMP) -o $(TARGET) -I$(userModPATH) -L$(userLibPATH) $(CFLAGS) $(OBJECTS) -l$(integrLIB) -l$(nrLIB)

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


	
	
