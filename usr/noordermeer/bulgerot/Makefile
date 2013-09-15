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
	rm -f *.o $(PROGRAM) sersic1.dat

install: all
	cp $(PROGRAM) $(NEMOBIN)
	mknemo runbulgerot

bench: all
	./bulgerot < bulgerot.in
	-diff sersic0.dat sersic1.dat


