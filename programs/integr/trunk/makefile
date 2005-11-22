# Makefile to produce a static library and modules related to integration of ODE's.
# by vasimos

# Choose compiler:
COMP =  
OPTFULL = -g
TARGET = libintegr
OBJ = ifc_integr.o derivsJ.o etdrk4Diag.o etdrk4DiagDriverS.o etdrk4DiagDriverS_a.o etdrk4DiagJ.o etdrk4DiagJif.o etdrk4DiagJDriverS.o etdrk4DiagPrefactors.o eulerJ.o integrPrk4.o odeint.o rk2J.o rk2Jdriver.o rk2Jsdriver.o rk4.o rk4driver.o rk4J.o rk4Jdriver.o rk4P.o rk4Pdriver.o rkck.o rkqs.o
srcMod = ifc_integr.f90 
targetMod = ifc_integr.mod ode_path.mod
# Where to find user libraries:
userLibPATH = /home/vasimos/lib-$(COMP)
# Where to find user modules:
userModPATH = $(userLibPATH)/modules
# Names of user libraries
nrLIB = nr
integrLIB= integr

$(TARGET): uninstall $(OBJ) $(srcMod) 
	make $(targetMod)
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
	cd $(userModPATH)/;\
	rm -f $(targetMod)
	rm -f $(userLibPATH)/$(TARGET).a

clean:
	rm -f *.o *.mod *.a
