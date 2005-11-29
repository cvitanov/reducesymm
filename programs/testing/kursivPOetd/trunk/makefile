# command line compiler options can be passed in CFLAGS

# Choose compiler
OBJECTS = RPOshootingTest.o parameters.o SetANdiag_KS.o SetLin_KS.o  SetNlin_KS.o 
# Set the outputfile
TARGET = ks 	
# Where to find user libraries:
userLibPATH = /home/vasimos/lib-$(COMP)
# Where to find user modules:
userModPATH = $(userLibPATH)/modules
# Where to find shared libraries: 
sharedLibPATH = /usr/lib64
# Where to find lapack:
lapackLibPATH = /usr/lib64/lapack95/$(COMP)
lapackModPATH = $(lapackLibPATH)/lapack95_modules
paramMod = parameters.mod
nrLIB = nr
integrLIB = integr
shootingLIB = po 
utilLIB = utils

ifdef COMP
include $(COMP).inc
endif

$(TARGET): $(OBJECTS)
	$(COMP) -o $(TARGET) -I$(userModPATH) -I$(lapackModPATH) -L$(laPATH) -L$(lapackLibPATH) -L$(userLibPATH) $(CFLAGS) $(OBJECTS) $(sharedLibPATH)/libfftw3.so -l$(shootingLIB) -l$(lapack95LIB) -l$(lapackLIB) -l$(blasLIB) -l$(integrLIB) -l$(utilLIB) -l$(nrLIB) 

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


	
	
