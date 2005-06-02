# command line compiler options can be passed in CFLAGS

# Choose compiler
COMP = pgf90	
OBJECTS = integr.o  oscillator.o
# Set the outputfile
TARGET = osc 	
# Where to find user libraries:
userLibPATH = /home/vasimos/lib64
# Where to find user modules:
userModPATH = $(userLibPATH)/modules
nrLIB = nr
integrLIB =integr


$(TARGET): $(OBJECTS)
	$(COMP) -o $(TARGET) -I$(userModPATH) -L$(userLibPATH) $(CFLAGS) $(OBJECTS) -l$(integrLIB) -l$(nrLIB)

%.o: %.f90 
	$(COMP) -c -o $@ -I$(userModPATH) $(CFLAGS) $< 

clean:
	rm -f *.o $(TARGET) 
	
cleanall:
	rm -f *.o *.mod $(TARGET) 

cleandat: 
	rm -f *.dat	

cleanedit: 
	rm -f *~


	
	