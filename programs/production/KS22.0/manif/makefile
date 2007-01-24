# command line compiler options can be passed in CFLAGS

# Choose compiler
OBJECTS = main.o
# Set the outputfile
TARGET = ks
# Where to find user libraries:
userLibPATH = /home/$(USER)/lib-$(COMP)
# Where to find user modules:
userModPATH = $(userLibPATH)/modules
# Where to find shared libraries: 
sharedLibPATH = /usr/lib$(local_arch)
# Where to find lapack:
lapackLibPATH = /usr/lib$(local_arch) 
lapackModPATH = $(lapack95path)/$(COMP)/lapack95_modules
nrLIB = nr
shootingLIB = newt
utilLIB = utils
integrLIB = integr
ksLIB= ks


# ifdef COMP
# include $(COMP).inc
# endif

$(TARGET): $(OBJECTS)
	$(COMP) -o $(TARGET) -I$(userModPATH) -I$(lapackModPATH) -L$(lapackLibPATH) -L$(lapack95path)/$(COMP) -L$(userLibPATH) $(CFLAGS) $(OBJECTS) $(sharedLibPATH)/libfftw3.so -l$(shootingLIB) -l$(ksLIB) -l$(integrLIB) -llapack95 -llapack -lblas -l$(utilLIB) -l$(nrLIB) /usr/lib/libgfortran.so.1

%.o: %.f90
	$(COMP) -c -o $@ -I$(userModPATH) -I$(lapackModPATH) $(CFLAGS) $< 


	
clean:
	rm -f *.o $(TARGET) 
	
cleanall:
	rm -f *.o *.mod $(TARGET) 

cleandat: 
	rm -f *.dat	

cleanedit: 
	rm -f *~


	
	
