include $(NEMOLIB)/makedefs

#FC            = f77
#LDFLAGS       = 

MAKEFILE      = Makefile

OBJS          = 

PROGRAM       = bulgerot

SRCS          = bulgerot.f

 
$(PROGRAM):     $(OBJS)  $(MAKEFILE) 
		$(FC) $(LDFLAGS)  $(SRCS)    -o $(PROGRAM)

all:    $(PROGRAM)        

clean: 
	rm -rf *.o $(PROGRAM) sersic1.dat sersic1D.dat bench0

install: all
	cp $(PROGRAM) $(NEMOBIN)
	mknemo runbulgerot

# classic fortran version
bench: all
	time ./bulgerot < bulgerot.in
	-diff sersic0.dat sersic1.dat

# this is the bench with nemo frontend
bench0:
	time runbulgerot bench0 galaxy=sersic:n=1
	diff sersic0.dat bench0/bulgerot.dat

# use this for older compilers who still write out the E format for REAL*8
bench1: sersic1.dat
	-sed s/E/D/g sersic1.dat > sersic1D.dat 
	-diff sersic0.dat sersic1D.dat


