# command line compiler options can be passed in CFLAGS

# Choose compiler
COMP = pgf90	
OBJECTS = main.o roesslerField.o roesslerFieldP.o roesslerVar.o
# Set the outputfile
TARGET = j	
# Where to find user libraries:
userLibPATH = /home/vasimos/lib64
# Where to find user modules:
userModPATH = $(userLibPATH)/modules
paramMod = parameters.mod
paramSRC = parameters.f90
nrLIB = nr
integrLIB = integr
utilLIB = util

$(TARGET): $(OBJECTS)
	$(COMP) -o $(TARGET) -I$(userModPATH) -L$(userLibPATH) $(CFLAGS) $(OBJECTS) -l$(integrLIB) -l$(utilLIB) -l$(nrLIB)

%.o: %.f90 $(paramMod)
	$(COMP) -c -o $@ -I$(userModPATH) $(CFLAGS) $< 

$(paramMod): $(paramSRC)
	$(COMP) -c -I$(userModPATH) $(CFLAGS) $< 
	
clean:
	rm -f *.o $(TARGET) 
	
cleanall:
	rm -f *.o *.mod $(TARGET) 

cleandat: 
	rm -f *.dat	

cleanedit: 
	rm -f *~


	
	