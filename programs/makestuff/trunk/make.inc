# Makefile to produce a static library and modules related to integration of ODE's.
# by vasimos

# Choose compiler:
COMP = pgf90 
OPTFULL = -fastsse
# Where to find user libraries:
userLibPATH = /home/vasimos/lib
# Where to find user modules:
userModPATH = $(userLibPATH)/modules
# Names of user libraries
nrLIB = nr
integrLIB= integr

