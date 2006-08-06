include $(NEMOLIB)/makedefs

#FC            = f77
#LDFLAGS       = 

MAKEFILE      = Makefile

OBJS          = \
		locate.o \
		magalie.o 

PROGRAM       = magalie

SRCS          = \
		superbox.f \
		locate.f 
 
$(PROGRAM):     $(OBJS)  $(MAKEFILE) 
		$(FC) $(LDFLAGS)  $(OBJS)    -o $(PROGRAM).exe

pgf: 
	$(MAKE) "FC = pgf77" "FFLAGS = -O3"  $(PROGRAM)

DEC77: 
	$(MAKE) "FC = f77" "FFLAGS = -O3"  $(PROGRAM)

g77: 
	$(MAKE) "FC = g77" "FFLAGS = $(FFLAGS) -O6 -D GNUF77"   $(PROGRAM)

gfortran: 
	$(MAKE) "FC = gfortran"  $(PROGRAM)

intel:
	$(MAKE) "FC = ifort"     $(PROGRAM)

all:    $(PROGRAM)        


clean: 
	rm -f *.o $(PROGRAM).exe

install: all
	cp $(PROGRAM).exe $(NEMOBIN)

#FFLAGS        = -O6 -m486 
#FFLAGS        = -g -pg -ansi -pedantic -Wall






