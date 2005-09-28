# command line compiler options can be passed in CFLAGS

# Choose compiler
COMP = pgf90	
OBJECTS = main.o parameters.o SetANdiag_roessler.o SetNlin_roessler.o SetLin_roessler.o
# Set the outputfile
TARGET = j	
# Where to find user libraries:
userLibPATH = /home/vasimos/lib64
# Where to find user modules:
userModPATH = $(userLibPATH)/modules
ll95PATH = /usr/pgi/linux86-64/5.2/lib
ml95PATH = /lib64/lapack95_modules
paramMod = parameters.mod
paramSRC = parameters.f90
nrLIB = nr
integrLIB = integr
utilLIB = util

$(TARGET): $(OBJECTS)
	$(COMP) -o $(TARGET) -I$(userModPATH) -I$(ml95PATH) -L$(userLibPATH) -L$(ll95PATH) $(CFLAGS) $(OBJECTS) -l$(integrLIB) -l$(utilLIB) -l$(nrLIB) -llapack95 -lacml

%.o: %.f90 $(paramMod)
	$(COMP) -c -o $@ -I$(userModPATH) -I$(ml95PATH) $(CFLAGS) $< 

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


	
	
