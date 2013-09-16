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


#   sersic1.dat is to compare with the benchmark in sersic0.dat
clean: 
	rm -f *.o $(PROGRAM) sersic1.dat sersic1D.dat

install: all
	cp $(PROGRAM) $(NEMOBIN)
	mknemo runbulgerot

bench: all
	./bulgerot < bulgerot.in
	-diff sersic0.dat sersic1.dat

# use this for older compilers who still write out the E format for REAL*8
bench1: sersic1.dat
	-sed s/E/D/g sersic1.dat > sersic1D.dat 
	-diff sersic0.dat sersic1D.dat


