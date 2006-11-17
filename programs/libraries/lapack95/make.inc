#	Modified: Nov 17 2006 by Evangelos Siminos
#  -- LAPACK95 interface driver routine (version 2.0) --
#     UNI-C, Denmark; Univ. of Tennessee, USA; NAG Ltd., UK
#     August 5, 2000
##	Customized for ifort compiler

FC	 = $(COMP)
FC1      = $(COMP)

OPTS0    = $(CFLAGS)
MODLIB   = -I./../lapack95_modules
OPTS1    = -c $(OPTS0)
OPTS3    = $(OPTS1) $(MODLIB)
OPTL     = -o 
OPTLIB   =

LAPACK_PATH = /usr/lib64

LAPACK95 = ../lapack95.a
LAPACK77 = $(LAPACK_PATH)/lapack.a
TMG77    = /usr/lib64/tmglib.a
BLAS     = $(LAPACK_PATH)/blas.a

LIBS     = $(LAPACK95) $(TMG77) $(LAPACK77) $(BLAS)
SUF      = f90

XX = 'rm' -f $@; \
        'rm' -f $@.res; \
	$(FC) $(OPTS0) -o $@ $(MODLIB) $@.$(SUF) $(OPTLIB) $(LIBS); \
        $@ < $@.dat > $@.res; \
        'rm' -f $@

YY = $(FC) $(OPTS0) -o $@ $(MODLIB) $@.$(SUF) $(OPTLIB) $(LIBS)

.SUFFIXES: .f90 .f .o

.$(SUF).o: 
	$(FC) $(OPTS3) $<

.f.o:
	$(FC1) $(OPTS3) $<

# To install system-wide: 	$> cp -r lapack95_modules /usr/lib64 
#				$> cp lapack95.a /usr/lib64/liblapack95.a
# Compile with something like
# $> ifort -L/usr/lib64 -I/usr/lib64/lapack95_modules driver.f90  -llapack95 -llapack -lblas      

