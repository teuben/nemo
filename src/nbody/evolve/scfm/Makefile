#	SCFM
include $(NEMOLIB)/makedefs
#	makefile for Lars' scfm version  (Version 1)
#
#   7-jul-97	installed from Lars version 1 in usr tree
#  20-jan-98    installed in src tree
#   4-apr-01    NEMO V3 changes

NBODY = 1024

#FFLAGS = -O
#CFLAGS = -O
#				not actually needed yet
L = $(NEMOLIB)/libnemo.a

ALLSRC = scfm.f second.f \
  README Makefile 
INC = scf.h
MAN =
EXE =

#				second.* may not be needed
SRC = scfm.f
OBJ = scfm.o

#		replace this with second0.o, second.o or whichever
#		if your make(1) utility cannot handle 'ifeq' constructs
#SECOND = second.o
ifeq "$(FC)" "g77"
  SECOND = 
else
  SECOND = second.o
endif

BINS = scfm.exe runscfm

TAR = $(ALLSRC) $(INC) $(MAN) $(EXE) SCFBI SCFPAR

help:
	@echo No help yet

nemo_lib:
	@echo No nemo_lib in `pwd`

nemo_bin:	$(BINS)
	mv $? $(NEMOBIN)
	rm -f *.o

bins: $(BINS)
	mv $(BINS) $(NEMOBIN)

nemo_src:
	-@echo `pwd` scfm
		
		
install:	nemo_bin

.install_src:	$(ALLSRC) $(INC)
	cp $? $(NEMO)/src/nbody/evolve/scfm
	@touch .install_src

.install_man:	$(MAN)
	cp $? $(NEMO)/man/man1
	@touch .install_man

scfm:	scfm.exe

scfm.exe: $(OBJ) $(SECOND)
	$(FC) $(FFLAGS) -o scfm.exe $(OBJ) $(SECOND)

runscfm: runscfm.c
	$(CC) $(CFLAGS) -o runscfm runscfm.c $(NEMO_LIBS)

clean:
	rm -fr $(BINS) core *.o SNAP??? SCFLOG SCFOUT SCFEL??? scfm.tar bench

tar:
	@echo Creating some sort of TAR file with all the essentials
	@echo Plus a little test file. Assume 2 example SCFBI and SCFPAR
	@echo files are present.
	@tar cvf scfm.tar $(TAR)


# 	don't do this, we like to keep a standard dataset
newdata:
	mkplummer - $(NBODY) seed=$(NBODY) |\
	 snapprint - m,x,y,z,vx,vy,vz header=t > SCFBI
#
bench: SCFBI SCFPAR scfm.exe
	mkdir bench; cd bench; cp ../SCFBI ../SCFPAR . ; time ../scfm.exe
