#	NEMO/src/nbody/io
include $(NEMOLIB)/makedefs

#	8-jun-88	Version 1.0	Piet Hut
#	21-nov-89	HOSTTYPE added	PJT
#       5-feb-01        adopted for current NEMO 	PJT

#..............................................................................

#	SRCDIR is the full directory name in nemo:
SRCDIR = $(NEMO)/src/piet/newton0/vanilla

# CFLAGS = -O
TREECFLAGS = -c -DTREE $(CFLAGS)
TREEQCFLAGS = -c -DTREE -DQUADPOLE $(CFLAGS)
REGCFLAGS = -c -DREGULARIZATION $(CFLAGS)
EXTCFLAGS = -c -DEXTRAPOLATION $(CFLAGS)

DBXCFLAGS =  -g
DBXTREECFLAGS = -c -DTREE $(DBXCFLAGS)
DBXTREEQCFLAGS = -c  -DTREE -DQUADPOLE $(DBXCFLAGS)
DBXREGCFLAGS = -c  -DREGULARIZATION $(DBXCFLAGS)
DBXEXTCFLAGS = -c  -DEXTRAPOLATION $(DBXCFLAGS)

L = $(NEMOLIB)/libnemo.a 

#

OBJ=binaryin.o binaryout.o \
        bodyalgebra.o bodyconversion.o create.o \
        diagnose.o differentiate.o integrate.o \
        newton0.o orbit.o out.o save.o soften.o \
        statealgebra.o systemalgebra.o \
        systemconversion.o timestep.o
TREEOBJ=binaryin_t.o binaryout_t.o \
            bodyalgebra_t.o bodyconversion_t.o create_t.o \
            diagnose_t.o differentiate_t.o differentiatetree_t.o \
            integrate_t.o transformtree_t.o newton0_t.o orbit_t.o \
            out_t.o save_t.o soften_t.o \
            statealgebra_t.o systemalgebra_t.o \
            systemconversion_t.o timestep_t.o
TREEQOBJ=binaryin_tq.o binaryout_tq.o \
             bodyalgebra_tq.o bodyconversion_tq.o create_tq.o \
             diagnose_tq.o differentiate_tq.o differentiatetree_tq.o \
             integrate_tq.o transformtree_tq.o newton0_tq.o orbit_tq.o \
             out_tq.o save_tq.o soften_tq.o \
             statealgebra_tq.o systemalgebra_tq.o \
             systemconversion_tq.o timestep_tq.o
REGOBJ=binaryin_r.o binaryout_r.o \
           bodyalgebra_r.o bodyconversion_r.o create_r.o \
           diagnose_r.o differentiate_r.o differentiatereg_r.o \
           integrate_r.o newton0_r.o orbit_r.o \
           out_r.o save_r.o soften_r.o \
           statealgebra_r.o systemalgebra_r.o \
           systemconversion_r.o timestep_r.o transformreg_r.o
EXTOBJ=binaryin_e.o binaryout_e.o \
           bodyalgebra_e.o bodyconversion_e.o create_e.o \
           diagnose_e.o differentiate_e.o differentiateext_e.o \
           integrate_e.o newton0_e.o orbit_e.o \
           out_e.o save_e.o soften_e.o \
           statealgebra_e.o systemalgebra_e.o \
           systemconversion_e.o timestep_e.o transformext_e.o

#

DBXOBJ=dbxbinaryin.o dbxbinaryout.o \
           dbxbodyalgebra.o dbxbodyconversion.o dbxcreate.o \
           dbxdiagnose.o dbxdifferentiate.o dbxintegrate.o \
           dbxnewton0.o dbxorbit.o dbxout.o dbxsave.o dbxsoften.o \
           dbxstatealgebra.o dbxsystemalgebra.o \
           dbxsystemconversion.o dbxtimestep.o
DBXTREEOBJ=dbxbinaryin_t.o dbxbinaryout_t.o \
               dbxbodyalgebra_t.o dbxbodyconversion_t.o dbxcreate_t.o \
               dbxdiagnose_t.o dbxdifferentiate_t.o dbxdifferentiatetree_t.o \
               dbxintegrate_t.o dbxtransformtree_t.o \
               dbxnewton0_t.o dbxorbit_t.o \
               dbxout_t.o dbxsave_t.o dbxsoften_t.o dbxstatealgebra_t.o \
               dbxsystemalgebra_t.o dbxsystemconversion_t.o dbxtimestep_t.o
DBXTREEQOBJ=dbxbinaryin_tq.o dbxbinaryout_tq.o dbxbodyalgebra_tq.o \
                dbxbodyconversion_tq.o dbxcreate_tq.o dbxdiagnose_tq.o \
                dbxdifferentiate_tq.o dbxdifferentiatetree_tq.o \
                dbxintegrate_tq.o dbxtransformtree_tq.o \
                dbxnewton0_tq.o dbxorbit_tq.o \
                dbxout_tq.o dbxsave_tq.o dbxsoften_tq.o dbxstatealgebra_tq.o \
                dbxsystemalgebra_tq.o dbxsystemconversion_tq.o dbxtimestep_tq.o
DBXREGOBJ=dbxbinaryin_r.o dbxbinaryout_r.o dbxbodyalgebra_r.o \
              dbxbodyconversion_r.o dbxcreate_r.o \
              dbxdiagnose_r.o dbxdifferentiate_r.o dbxdifferentiatereg_r.o \
              dbxintegrate_r.o dbxnewton0_r.o dbxorbit_r.o \
              dbxout_r.o dbxsave_r.o dbxsoften_r.o dbxstatealgebra_r.o \
              dbxsystemalgebra_r.o dbxsystemconversion_r.o \
              dbxtimestep_r.o dbxtransformreg_r.o
DBXEXTOBJ=dbxbinaryin_e.o dbxbinaryout_e.o dbxbodyalgebra_e.o \
              dbxbodyconversion_e.o dbxcreate_e.o \
              dbxdiagnose_e.o dbxdifferentiate_e.o dbxdifferentiateext_e.o \
              dbxintegrate_e.o dbxnewton0_e.o dbxorbit_e.o \
              dbxout_e.o dbxsave_e.o dbxsoften_e.o dbxstatealgebra_e.o \
              dbxsystemalgebra_e.o dbxsystemconversion_e.o \
              dbxtimestep_e.o dbxtransformext_e.o

#

TOOLS= newton0 newton0tree newton0treeq newton0reg newton0ext
DBXTOOLS= dbxnewton0 dbxnewton0tree dbxnewton0treeq dbxnewton0reg dbxnewton0ext

INC = $(STDINC) $(TREEINC) $(REGINC) $(EXTINC)
STDINC = body.h ctrl.h diag.h newton0.h spec.h state.h syst.h
TREEINC = node.h 
REGINC = 
EXTINC = bodyext.h 

SRC = Makefile $(INC) \
      $(STDSRC) $(TREESRC) $(REGSRC) $(EXTSRC)
STDSRC = binaryin.c binaryout.c \
        bodyalgebra.c bodyconversion.c create.c \
        diagnose.c differentiate.c integrate.c \
        newton0.c orbit.c out.c save.c soften.c \
        statealgebra.c systemalgebra.c \
        systemconversion.c timestep.c
TREESRC = differentiatetree.c transformtree.c
REGSRC = differentiatereg.c transformreg.c
EXTSRC = differentiateext.c transformext.c

MAN1 = newton0.1

#==============================================================================

help:
	@echo This is one of piet's  Makefiles, residing within nemo in
	@echo \" $(SRCDIR) \" .
	@echo
	@echo This Makefile offers the following make options:
	@echo install: invokes all of the following install-... options:
	@echo install_inc: install all header files
	@echo install_src: install all source files
	@echo install_bin: install all tools
	@echo install_man: install all man pages
	@echo cleanup: removes the non-ASCII files

#

tools: $(TOOLS)
	

test: $(TOOLS)
	-@for i in $(TOOLS); do\
	(echo "### $$i";$$i); done

	

#==============================================================================
#
#	WARNING: the following commands should NEVER by issued by
#                hand, but ONLY automatically by nemo from within nemo!
#..............................................................................
#
#	The following commands are normally called from one level up in nemo:

nemo_lib:
	@echo No nemo_lib here

nemo_bin:	$(TOOLS)
	mv $? $(NEMOBIN)
	rm *.o

nemo_src:
	-@for i in $(TOOLS); do \
	echo `pwd` $$i ; done 
		
clean:
	rm -f core *.o $(TOOLS) $(DBXTOOLS)
#
#==============================================================================
#
#	USER area: the following commands can be issued freely by the user:

install:  install_inc install_src install_lib install_bin install_man

install_inc:	.install_inc
install_src:	.install_src
install_lib:	.install_lib
install_bin:	.install_bin
install_man:	.install_man

.install_inc: $(INC)
	/bin/cp $?  $(NEMOPATH)/inc
	@touch .install_inc

.install_src: $(SRC)
	@if [ ! -d $(SRCDIR) ];  \
 	then \
 		mkdir $(SRCDIR); \
 		chmod 777 $(SRCDIR); \
 	fi
	/bin/cp $? $(SRCDIR)
	@touch .install_src

.install_lib:
	@echo No lib here

.install_bin: $(TOOLS)
	/bin/cp $? $(NEMOBIN)
	@touch .install_bin

.install_man:	.install_man1

.install_man1:	$(MAN1)
	/bin/cp $? $(NEMOPATH)/man/man1
	@touch .install_man1

#..............................................................................

cleanup:
	rm *.o *.a $(TOOLS) $(DBXTOOLS)

#..............................................................................
#
newton0: $(OBJ) 
	 $(CC) -o newton0 $(OBJ) $(L) -lm

binaryin.o: binaryin.c $(STDINC)

binaryout.o: binaryout.c $(STDINC)

bodyalgebra.o: bodyalgebra.c $(STDINC)

bodyconversion.o: bodyconversion.c $(STDINC)

create.o: create.c $(STDINC)

diagnose.o: diagnose.c $(STDINC)

differentiate.o: differentiate.c $(STDINC)

integrate.o: integrate.c $(STDINC)

newton0.o: newton0.c $(STDINC)

orbit.o: orbit.c $(STDINC)

out.o: out.c $(STDINC)

save.o: save.c $(STDINC)

soften.o: soften.c $(STDINC)

statealgebra.o: statealgebra.c $(STDINC)

systemalgebra.o: systemalgebra.c $(STDINC)

systemconversion.o: systemconversion.c $(STDINC)

timestep.o: timestep.c $(STDINC)

#..............................................................................
#
dbxnewton0: $(DBXOBJ)
	 $(CC) -o dbxnewton0 $(DBXOBJ) $(L) -lm

dbxbinaryin.o: binaryin.c $(STDINC)
	$(CC) -c $(DBXCFLAGS) -o dbxbinaryin.o binaryin.c

dbxbinaryout.o: binaryout.c $(STDINC)
	$(CC) -c $(DBXCFLAGS) -o dbxbinaryout.o binaryout.c

dbxbodyalgebra.o: bodyalgebra.c $(STDINC)
	$(CC) -c $(DBXCFLAGS) -o dbxbodyalgebra.o bodyalgebra.c 

dbxbodyconversion.o: bodyconversion.c $(STDINC)
	$(CC) -c $(DBXCFLAGS) -o dbxbodyconversion.o bodyconversion.c 

dbxcreate.o: create.c $(STDINC)
	$(CC) -c $(DBXCFLAGS) -o dbxcreate.o create.c 

dbxdiagnose.o: diagnose.c $(STDINC)
	$(CC) -c $(DBXCFLAGS) -o dbxdiagnose.o diagnose.c

dbxdifferentiate.o: differentiate.c $(STDINC)
	$(CC) -c $(DBXCFLAGS) -o dbxdifferentiate.o differentiate.c

dbxintegrate.o: integrate.c $(STDINC)
	$(CC) -c $(DBXCFLAGS) -o dbxintegrate.o integrate.c

dbxnewton0.o: newton0.c $(STDINC)
	$(CC) -c $(DBXCFLAGS) -o dbxnewton0.o newton0.c

dbxorbit.o: orbit.c $(STDINC)
	$(CC) -c $(DBXCFLAGS) -o dbxorbit.o orbit.c

dbxout.o: out.c $(STDINC)
	$(CC) -c $(DBXCFLAGS) -o dbxout.o out.c

dbxsave.o: save.c $(STDINC)
	$(CC) -c $(DBXCFLAGS) -o dbxsave.o save.c

dbxsoften.o: soften.c $(STDINC)
	$(CC) -c $(DBXCFLAGS) -o dbxsoften.o soften.c

dbxstatealgebra.o: statealgebra.c $(STDINC)
	$(CC) -c $(DBXCFLAGS) -o dbxstatealgebra.o statealgebra.c 

dbxsystemalgebra.o: systemalgebra.c $(STDINC)
	$(CC) -c $(DBXCFLAGS) -o dbxsystemalgebra.o systemalgebra.c 

dbxsystemconversion.o: systemconversion.c $(STDINC)
	$(CC) -c $(DBXCFLAGS) -o dbxsystemconversion.o systemconversion.c 

dbxtimestep.o: timestep.c $(STDINC)
	$(CC) -c $(DBXCFLAGS) -o dbxtimestep.o timestep.c

#..............................................................................
#
newton0tree: $(TREEOBJ)
	 $(CC) -o newton0tree $(TREEOBJ) $(L) -lm

binaryin_t.o: binaryin.c $(STDINC) $(TREEINC)
	$(CC) $(TREECFLAGS) -o binaryin_t.o binaryin.c

binaryout_t.o: binaryout.c $(STDINC) $(TREEINC)
	$(CC) $(TREECFLAGS) -o binaryout_t.o binaryout.c

bodyalgebra_t.o: bodyalgebra.c $(STDINC) $(TREEINC)
	$(CC) $(TREECFLAGS) -o bodyalgebra_t.o bodyalgebra.c 

bodyconversion_t.o: bodyconversion.c $(STDINC) $(TREEINC)
	$(CC) $(TREECFLAGS) -o bodyconversion_t.o bodyconversion.c 

create_t.o: create.c $(STDINC) $(TREEINC)
	$(CC) $(TREECFLAGS) -o create_t.o create.c 

diagnose_t.o: diagnose.c $(STDINC) $(TREEINC)
	$(CC) $(TREECFLAGS) -o diagnose_t.o diagnose.c

differentiate_t.o: differentiate.c $(STDINC) $(TREEINC)
	$(CC) $(TREECFLAGS) -o differentiate_t.o differentiate.c

differentiatetree_t.o: differentiatetree.c $(STDINC) $(TREEINC)
	$(CC) $(TREECFLAGS) -o differentiatetree_t.o differentiatetree.c

integrate_t.o: integrate.c $(STDINC) $(TREEINC)
	$(CC) $(TREECFLAGS) -o integrate_t.o integrate.c

transformtree_t.o: transformtree.c $(STDINC) $(TREEINC)
	$(CC) $(TREECFLAGS) -o transformtree_t.o transformtree.c

newton0_t.o: newton0.c $(STDINC) $(TREEINC)
	$(CC) $(TREECFLAGS) -o newton0_t.o newton0.c

orbit_t.o: orbit.c $(STDINC) $(TREEINC)
	$(CC) $(TREECFLAGS) -o orbit_t.o orbit.c

out_t.o: out.c $(STDINC) $(TREEINC)
	$(CC) $(TREECFLAGS) -o out_t.o out.c

save_t.o: save.c $(STDINC) $(TREEINC)
	$(CC) $(TREECFLAGS) -o save_t.o save.c

soften_t.o: soften.c $(STDINC) $(TREEINC)
	$(CC) $(TREECFLAGS) -o soften_t.o soften.c

statealgebra_t.o: statealgebra.c $(STDINC) $(TREEINC)
	$(CC) $(TREECFLAGS) -o statealgebra_t.o statealgebra.c 

systemalgebra_t.o: systemalgebra.c $(STDINC) $(TREEINC)
	$(CC) $(TREECFLAGS) -o systemalgebra_t.o systemalgebra.c 

systemconversion_t.o: systemconversion.c $(STDINC) $(TREEINC)
	$(CC) $(TREECFLAGS) -o systemconversion_t.o systemconversion.c 

timestep_t.o: timestep.c $(STDINC) $(TREEINC)
	$(CC) $(TREECFLAGS) -o timestep_t.o timestep.c

#..............................................................................
#
dbxnewton0tree: $(DBXTREEOBJ)
	 $(CC) -o dbxnewton0tree $(DBXTREEOBJ)  $(L) -lm

dbxbinaryin_t.o: binaryin.c $(STDINC) $(TREEINC)
	$(CC) $(DBXTREECFLAGS) -o dbxbinaryin_t.o binaryin.c

dbxbinaryout_t.o: binaryout.c $(STDINC) $(TREEINC)
	$(CC) $(DBXTREECFLAGS) -o dbxbinaryout_t.o binaryout.c

dbxbodyalgebra_t.o: bodyalgebra.c $(STDINC) $(TREEINC)
	$(CC) $(DBXTREECFLAGS) -o dbxbodyalgebra_t.o bodyalgebra.c 

dbxbodyconversion_t.o: bodyconversion.c $(STDINC) $(TREEINC)
	$(CC) $(DBXTREECFLAGS) -o dbxbodyconversion_t.o bodyconversion.c 

dbxcreate_t.o: create.c $(STDINC) $(TREEINC)
	$(CC) $(DBXTREECFLAGS) -o dbxcreate_t.o create.c 

dbxdiagnose_t.o: diagnose.c $(STDINC) $(TREEINC)
	$(CC) $(DBXTREECFLAGS) -o dbxdiagnose_t.o diagnose.c

dbxdifferentiate_t.o: differentiate.c $(STDINC) $(TREEINC)
	$(CC) $(DBXTREECFLAGS) -o dbxdifferentiate_t.o differentiate.c

dbxdifferentiatetree_t.o: differentiatetree.c $(STDINC) $(TREEINC)
	$(CC) $(DBXTREECFLAGS) -o dbxdifferentiatetree_t.o differentiatetree.c

dbxintegrate_t.o: integrate.c $(STDINC) $(TREEINC)
	$(CC) $(DBXTREECFLAGS) -o dbxintegrate_t.o integrate.c

dbxtransformtree_t.o: transformtree.c $(STDINC) $(TREEINC)
	$(CC) $(DBXTREECFLAGS) -o dbxtransformtree_t.o transformtree.c

dbxnewton0_t.o: newton0.c $(STDINC) $(TREEINC)
	$(CC) $(DBXTREECFLAGS) -o dbxnewton0_t.o newton0.c

dbxorbit_t.o: orbit.c $(STDINC) $(TREEINC)
	$(CC) $(DBXTREECFLAGS) -o dbxorbit_t.o orbit.c

dbxout_t.o: out.c $(STDINC) $(TREEINC)
	$(CC) $(DBXTREECFLAGS) -o dbxout_t.o out.c

dbxsave_t.o: save.c $(STDINC) $(TREEINC)
	$(CC) $(DBXTREECFLAGS) -o dbxsave_t.o save.c

dbxsoften_t.o: soften.c $(STDINC) $(TREEINC)
	$(CC) $(DBXTREECFLAGS) -o dbxsoften_t.o soften.c

dbxstatealgebra_t.o: statealgebra.c $(STDINC) $(TREEINC)
	$(CC) $(DBXTREECFLAGS) -o dbxstatealgebra_t.o statealgebra.c 

dbxsystemalgebra_t.o: systemalgebra.c $(STDINC) $(TREEINC)
	$(CC) $(DBXTREECFLAGS) -o dbxsystemalgebra_t.o systemalgebra.c 

dbxsystemconversion_t.o: systemconversion.c $(STDINC) $(TREEINC)
	$(CC) $(DBXTREECFLAGS) -o dbxsystemconversion_t.o systemconversion.c 

dbxtimestep_t.o: timestep.c $(STDINC) $(TREEINC)
	$(CC) $(DBXTREECFLAGS) -o dbxtimestep_t.o timestep.c

#..............................................................................
#
newton0treeq: $(TREEQOBJ)
	 $(CC) -o newton0treeq $(TREEQOBJ) $(L) -lm

binaryin_tq.o: binaryin.c $(STDINC) $(TREEINC)
	$(CC) $(TREEQCFLAGS) -o binaryin_tq.o binaryin.c

binaryout_tq.o: binaryout.c $(STDINC) $(TREEINC)
	$(CC) $(TREEQCFLAGS) -o binaryout_tq.o binaryout.c

bodyalgebra_tq.o: bodyalgebra.c $(STDINC) $(TREEINC)
	$(CC) $(TREEQCFLAGS) -o bodyalgebra_tq.o bodyalgebra.c 

bodyconversion_tq.o: bodyconversion.c $(STDINC) $(TREEINC)
	$(CC) $(TREEQCFLAGS) -o bodyconversion_tq.o bodyconversion.c 

create_tq.o: create.c $(STDINC) $(TREEINC)
	$(CC) $(TREEQCFLAGS) -o create_tq.o create.c 

diagnose_tq.o: diagnose.c $(STDINC) $(TREEINC)
	$(CC) $(TREEQCFLAGS) -o diagnose_tq.o diagnose.c

differentiate_tq.o: differentiate.c $(STDINC) $(TREEINC)
	$(CC) $(TREEQCFLAGS) -o differentiate_tq.o differentiate.c

differentiatetree_tq.o: differentiatetree.c $(STDINC) $(TREEINC)
	$(CC) $(TREEQCFLAGS) -o differentiatetree_tq.o differentiatetree.c

integrate_tq.o: integrate.c $(STDINC) $(TREEINC)
	$(CC) $(TREEQCFLAGS) -o integrate_tq.o integrate.c

transformtree_tq.o: transformtree.c $(STDINC) $(TREEINC)
	$(CC) $(TREEQCFLAGS) -o transformtree_tq.o transformtree.c

newton0_tq.o: newton0.c $(STDINC) $(TREEINC)
	$(CC) $(TREEQCFLAGS) -o newton0_tq.o newton0.c

orbit_tq.o: orbit.c $(STDINC) $(TREEINC)
	$(CC) $(TREEQCFLAGS) -o orbit_tq.o orbit.c

out_tq.o: out.c $(STDINC) $(TREEINC)
	$(CC) $(TREEQCFLAGS) -o out_tq.o out.c

save_tq.o: save.c $(STDINC) $(TREEINC)
	$(CC) $(TREEQCFLAGS) -o save_tq.o save.c

soften_tq.o: soften.c $(STDINC) $(TREEINC)
	$(CC) $(TREEQCFLAGS) -o soften_tq.o soften.c

statealgebra_tq.o: statealgebra.c $(STDINC) $(TREEINC)
	$(CC) $(TREEQCFLAGS) -o statealgebra_tq.o statealgebra.c 

systemalgebra_tq.o: systemalgebra.c $(STDINC) $(TREEINC)
	$(CC) $(TREEQCFLAGS) -o systemalgebra_tq.o systemalgebra.c 

systemconversion_tq.o: systemconversion.c $(STDINC) $(TREEINC)
	$(CC) $(TREEQCFLAGS) -o systemconversion_tq.o systemconversion.c 

timestep_tq.o: timestep.c $(STDINC) $(TREEINC)
	$(CC) $(TREEQCFLAGS) -o timestep_tq.o timestep.c

#..............................................................................
#
dbxnewton0treeq: $(DBXTREEQOBJ)
	 $(CC) -o dbxnewton0treeq $(DBXTREEQOBJ) $(L) -lm

dbxbinaryin_tq.o: binaryin.c $(STDINC) $(TREEINC)
	$(CC) $(DBXTREEQCFLAGS) -o dbxbinaryin_tq.o binaryin.c

dbxbinaryout_tq.o: binaryout.c $(STDINC) $(TREEINC)
	$(CC) $(DBXTREEQCFLAGS) -o dbxbinaryout_tq.o binaryout.c

dbxbodyalgebra_tq.o: bodyalgebra.c $(STDINC) $(TREEINC)
	$(CC) $(DBXTREEQCFLAGS) -o dbxbodyalgebra_tq.o bodyalgebra.c 

dbxbodyconversion_tq.o: bodyconversion.c $(STDINC) $(TREEINC)
	$(CC) $(DBXTREEQCFLAGS) -o dbxbodyconversion_tq.o bodyconversion.c 

dbxcreate_tq.o: create.c $(STDINC) $(TREEINC)
	$(CC) $(DBXTREEQCFLAGS) -o dbxcreate_tq.o create.c 

dbxdiagnose_tq.o: diagnose.c $(STDINC) $(TREEINC)
	$(CC) $(DBXTREEQCFLAGS) -o dbxdiagnose_tq.o diagnose.c

dbxdifferentiate_tq.o: differentiate.c $(STDINC) $(TREEINC)
	$(CC) $(DBXTREEQCFLAGS) -o dbxdifferentiate_tq.o differentiate.c

dbxdifferentiatetree_tq.o: differentiatetree.c $(STDINC) $(TREEINC)
	$(CC) $(DBXTREEQCFLAGS) -o dbxdifferentiatetree_tq.o differentiatetree.c

dbxintegrate_tq.o: integrate.c $(STDINC) $(TREEINC)
	$(CC) $(DBXTREEQCFLAGS) -o dbxintegrate_tq.o integrate.c

dbxtransformtree_tq.o: transformtree.c $(STDINC) $(TREEINC)
	$(CC) $(DBXTREEQCFLAGS) -o dbxtransformtree_tq.o transformtree.c

dbxnewton0_tq.o: newton0.c $(STDINC) $(TREEINC)
	$(CC) $(DBXTREEQCFLAGS) -o dbxnewton0_tq.o newton0.c

dbxorbit_tq.o: orbit.c $(STDINC) $(TREEINC)
	$(CC) $(DBXTREEQCFLAGS) -o dbxorbit_tq.o orbit.c

dbxout_tq.o: out.c $(STDINC) $(TREEINC)
	$(CC) $(DBXTREEQCFLAGS) -o dbxout_tq.o out.c

dbxsave_tq.o: save.c $(STDINC) $(TREEINC)
	$(CC) $(DBXTREEQCFLAGS) -o dbxsave_tq.o save.c

dbxsoften_tq.o: soften.c $(STDINC) $(TREEINC)
	$(CC) $(DBXTREEQCFLAGS) -o dbxsoften_tq.o soften.c

dbxstatealgebra_tq.o: statealgebra.c $(STDINC) $(TREEINC)
	$(CC) $(DBXTREEQCFLAGS) -o dbxstatealgebra_tq.o statealgebra.c 

dbxsystemalgebra_tq.o: systemalgebra.c $(STDINC) $(TREEINC)
	$(CC) $(DBXTREEQCFLAGS) -o dbxsystemalgebra_tq.o systemalgebra.c 

dbxsystemconversion_tq.o: systemconversion.c $(STDINC) $(TREEINC)
	$(CC) $(DBXTREEQCFLAGS) -o dbxsystemconversion_tq.o systemconversion.c 

dbxtimestep_tq.o: timestep.c $(STDINC) $(TREEINC)
	$(CC) $(DBXTREEQCFLAGS) -o dbxtimestep_tq.o timestep.c

#..............................................................................
#
newton0reg: $(REGOBJ)
	 $(CC) -o newton0reg $(REGOBJ) $(L) -lm

binaryin_r.o: binaryin.c $(STDINC) $(REGINC)
	$(CC) $(REGCFLAGS) -o binaryin_r.o binaryin.c

binaryout_r.o: binaryout.c $(STDINC) $(REGINC)
	$(CC) $(REGCFLAGS) -o binaryout_r.o binaryout.c

bodyalgebra_r.o: bodyalgebra.c $(STDINC) $(REGINC)
	$(CC) $(REGCFLAGS) -o bodyalgebra_r.o bodyalgebra.c 

bodyconversion_r.o: bodyconversion.c $(STDINC) $(REGINC)
	$(CC) $(REGCFLAGS) -o bodyconversion_r.o bodyconversion.c 

create_r.o: create.c $(STDINC) $(REGINC)
	$(CC) $(REGCFLAGS) -o create_r.o create.c 

diagnose_r.o: diagnose.c $(STDINC) $(REGINC)
	$(CC) $(REGCFLAGS) -o diagnose_r.o diagnose.c

differentiate_r.o: differentiate.c $(STDINC) $(REGINC)
	$(CC) $(REGCFLAGS) -o differentiate_r.o differentiate.c

differentiatereg_r.o: differentiatereg.c $(STDINC) $(REGINC)
	$(CC) $(REGCFLAGS) -o differentiatereg_r.o differentiatereg.c

integrate_r.o: integrate.c $(STDINC) $(REGINC)
	$(CC) $(REGCFLAGS) -o integrate_r.o integrate.c

newton0_r.o: newton0.c $(STDINC) $(REGINC)
	$(CC) $(REGCFLAGS) -o newton0_r.o newton0.c

orbit_r.o: orbit.c $(STDINC) $(REGINC)
	$(CC) $(REGCFLAGS) -o orbit_r.o orbit.c

out_r.o: out.c $(STDINC) $(REGINC)
	$(CC) $(REGCFLAGS) -o out_r.o out.c

save_r.o: save.c $(STDINC) $(REGINC)
	$(CC) $(REGCFLAGS) -o save_r.o save.c

soften_r.o: soften.c $(STDINC) $(REGINC)
	$(CC) $(REGCFLAGS) -o soften_r.o soften.c

statealgebra_r.o: statealgebra.c $(STDINC) $(REGINC)
	$(CC) $(REGCFLAGS) -o statealgebra_r.o statealgebra.c 

systemalgebra_r.o: systemalgebra.c $(STDINC) $(REGINC)
	$(CC) $(REGCFLAGS) -o systemalgebra_r.o systemalgebra.c 

systemconversion_r.o: systemconversion.c $(STDINC) $(REGINC)
	$(CC) $(REGCFLAGS) -o systemconversion_r.o systemconversion.c 

timestep_r.o: timestep.c $(STDINC) $(REGINC)
	$(CC) $(REGCFLAGS) -o timestep_r.o timestep.c

transformreg_r.o: transformreg.c $(STDINC) $(REGINC)
	$(CC) $(REGCFLAGS) -o transformreg_r.o transformreg.c

#..............................................................................
#
dbxnewton0reg: $(DBXREGOBJ)
	 $(CC) -o dbxnewton0reg $(DBXREGOBJ) $(L) -lm

dbxbinaryin_r.o: binaryin.c $(STDINC) $(REGINC)
	$(CC) $(DBXREGCFLAGS) -o dbxbinaryin_r.o binaryin.c

dbxbinaryout_r.o: binaryout.c $(STDINC) $(REGINC)
	$(CC) $(DBXREGCFLAGS) -o dbxbinaryout_r.o binaryout.c

dbxbodyalgebra_r.o: bodyalgebra.c $(STDINC) $(REGINC)
	$(CC) $(DBXREGCFLAGS) -o dbxbodyalgebra_r.o bodyalgebra.c 

dbxbodyconversion_r.o: bodyconversion.c $(STDINC) $(REGINC)
	$(CC) $(DBXREGCFLAGS) -o dbxbodyconversion_r.o bodyconversion.c 

dbxcreate_r.o: create.c $(STDINC) $(REGINC)
	$(CC) $(DBXREGCFLAGS) -o dbxcreate_r.o create.c 

dbxdiagnose_r.o: diagnose.c $(STDINC) $(REGINC)
	$(CC) $(DBXREGCFLAGS) -o dbxdiagnose_r.o diagnose.c

dbxdifferentiate_r.o: differentiate.c $(STDINC) $(REGINC)
	$(CC) $(DBXREGCFLAGS) -o dbxdifferentiate_r.o differentiate.c

dbxdifferentiatereg_r.o: differentiatereg.c $(STDINC) $(REGINC)
	$(CC) $(DBXREGCFLAGS) -o dbxdifferentiatereg_r.o differentiatereg.c

dbxintegrate_r.o: integrate.c $(STDINC) $(REGINC)
	$(CC) $(DBXREGCFLAGS) -o dbxintegrate_r.o integrate.c

dbxnewton0_r.o: newton0.c $(STDINC) $(REGINC)
	$(CC) $(DBXREGCFLAGS) -o dbxnewton0_r.o newton0.c

dbxorbit_r.o: orbit.c $(STDINC) $(REGINC)
	$(CC) $(DBXREGCFLAGS) -o dbxorbit_r.o orbit.c

dbxout_r.o: out.c $(STDINC) $(REGINC)
	$(CC) $(DBXREGCFLAGS) -o dbxout_r.o out.c

dbxsave_r.o: save.c $(STDINC) $(REGINC)
	$(CC) $(DBXREGCFLAGS) -o dbxsave_r.o save.c

dbxsoften_r.o: soften.c $(STDINC) $(REGINC)
	$(CC) $(DBXREGCFLAGS) -o dbxsoften_r.o soften.c

dbxstatealgebra_r.o: statealgebra.c $(STDINC) $(REGINC)
	$(CC) $(DBXREGCFLAGS) -o dbxstatealgebra_r.o statealgebra.c 

dbxsystemalgebra_r.o: systemalgebra.c $(STDINC) $(REGINC)
	$(CC) $(DBXREGCFLAGS) -o dbxsystemalgebra_r.o systemalgebra.c 

dbxsystemconversion_r.o: systemconversion.c $(STDINC) $(REGINC)
	$(CC) $(DBXREGCFLAGS) -o dbxsystemconversion_r.o systemconversion.c 

dbxtimestep_r.o: timestep.c $(STDINC) $(REGINC)
	$(CC) $(DBXREGCFLAGS) -o dbxtimestep_r.o timestep.c

dbxtransformreg_r.o: transformreg.c $(STDINC) $(REGINC)
	$(CC) $(DBXREGCFLAGS) -o dbxtransformreg_r.o transformreg.c

#..............................................................................
#
newton0ext: $(EXTOBJ)
	 $(CC) -o newton0ext $(EXTOBJ) $(L) -lm

binaryin_e.o: binaryin.c $(STDINC) $(EXTINC)
	$(CC) $(EXTCFLAGS) -o binaryin_e.o binaryin.c

binaryout_e.o: binaryout.c $(STDINC) $(EXTINC)
	$(CC) $(EXTCFLAGS) -o binaryout_e.o binaryout.c

bodyalgebra_e.o: bodyalgebra.c $(STDINC) $(EXTINC)
	$(CC) $(EXTCFLAGS) -o bodyalgebra_e.o bodyalgebra.c 

bodyconversion_e.o: bodyconversion.c $(STDINC) $(EXTINC)
	$(CC) $(EXTCFLAGS) -o bodyconversion_e.o bodyconversion.c 

create_e.o: create.c $(STDINC) $(EXTINC)
	$(CC) $(EXTCFLAGS) -o create_e.o create.c 

diagnose_e.o: diagnose.c $(STDINC) $(EXTINC)
	$(CC) $(EXTCFLAGS) -o diagnose_e.o diagnose.c

differentiate_e.o: differentiate.c $(STDINC) $(EXTINC)
	$(CC) $(EXTCFLAGS) -o differentiate_e.o differentiate.c

differentiateext_e.o: differentiateext.c $(STDINC) $(EXTINC)
	$(CC) $(EXTCFLAGS) -o differentiateext_e.o differentiateext.c

integrate_e.o: integrate.c $(STDINC) $(EXTINC)
	$(CC) $(EXTCFLAGS) -o integrate_e.o integrate.c

newton0_e.o: newton0.c $(STDINC) $(EXTINC)
	$(CC) $(EXTCFLAGS) -o newton0_e.o newton0.c

orbit_e.o: orbit.c $(STDINC) $(EXTINC)
	$(CC) $(EXTCFLAGS) -o orbit_e.o orbit.c

out_e.o: out.c $(STDINC) $(EXTINC)
	$(CC) $(EXTCFLAGS) -o out_e.o out.c

save_e.o: save.c $(STDINC) $(EXTINC)
	$(CC) $(EXTCFLAGS) -o save_e.o save.c

soften_e.o: soften.c $(STDINC) $(EXTINC)
	$(CC) $(EXTCFLAGS) -o soften_e.o soften.c

statealgebra_e.o: statealgebra.c $(STDINC) $(EXTINC)
	$(CC) $(EXTCFLAGS) -o statealgebra_e.o statealgebra.c 

systemalgebra_e.o: systemalgebra.c $(STDINC) $(EXTINC)
	$(CC) $(EXTCFLAGS) -o systemalgebra_e.o systemalgebra.c 

systemconversion_e.o: systemconversion.c $(STDINC) $(EXTINC)
	$(CC) $(EXTCFLAGS) -o systemconversion_e.o systemconversion.c 

timestep_e.o: timestep.c $(STDINC) $(EXTINC)
	$(CC) $(EXTCFLAGS) -o timestep_e.o timestep.c

transformext_e.o: transformext.c $(STDINC) $(EXTINC)
	$(CC) $(EXTCFLAGS) -o transformext_e.o transformext.c

#..............................................................................
#
dbxnewton0ext: $(DBXEXTOBJ)
	 $(CC) -o dbxnewton0ext $(DBXEXTOBJ) $(L) -lm

dbxbinaryin_e.o: binaryin.c $(STDINC) $(EXTINC)
	$(CC) $(DBXEXTCFLAGS) -o dbxbinaryin_e.o binaryin.c

dbxbinaryout_e.o: binaryout.c $(STDINC) $(EXTINC)
	$(CC) $(DBXEXTCFLAGS) -o dbxbinaryout_e.o binaryout.c

dbxbodyalgebra_e.o: bodyalgebra.c $(STDINC) $(EXTINC)
	$(CC) $(DBXEXTCFLAGS) -o dbxbodyalgebra_e.o bodyalgebra.c 

dbxbodyconversion_e.o: bodyconversion.c $(STDINC) $(EXTINC)
	$(CC) $(DBXEXTCFLAGS) -o dbxbodyconversion_e.o bodyconversion.c 

dbxcreate_e.o: create.c $(STDINC) $(EXTINC)
	$(CC) $(DBXEXTCFLAGS) -o dbxcreate_e.o create.c 

dbxdiagnose_e.o: diagnose.c $(STDINC) $(EXTINC)
	$(CC) $(DBXEXTCFLAGS) -o dbxdiagnose_e.o diagnose.c

dbxdifferentiate_e.o: differentiate.c $(STDINC) $(EXTINC)
	$(CC) $(DBXEXTCFLAGS) -o dbxdifferentiate_e.o differentiate.c

dbxdifferentiateext_e.o: differentiateext.c $(STDINC) $(EXTINC)
	$(CC) $(DBXEXTCFLAGS) -o dbxdifferentiateext_e.o differentiateext.c

dbxintegrate_e.o: integrate.c $(STDINC) $(EXTINC)
	$(CC) $(DBXEXTCFLAGS) -o dbxintegrate_e.o integrate.c

dbxnewton0_e.o: newton0.c $(STDINC) $(EXTINC)
	$(CC) $(DBXEXTCFLAGS) -o dbxnewton0_e.o newton0.c

dbxorbit_e.o: orbit.c $(STDINC) $(EXTINC)
	$(CC) $(DBXEXTCFLAGS) -o dbxorbit_e.o orbit.c

dbxout_e.o: out.c $(STDINC) $(EXTINC)
	$(CC) $(DBXEXTCFLAGS) -o dbxout_e.o out.c

dbxsave_e.o: save.c $(STDINC) $(EXTINC)
	$(CC) $(DBXEXTCFLAGS) -o dbxsave_e.o save.c

dbxsoften_e.o: soften.c $(STDINC) $(EXTINC)
	$(CC) $(DBXEXTCFLAGS) -o dbxsoften_e.o soften.c

dbxstatealgebra_e.o: statealgebra.c $(STDINC) $(EXTINC)
	$(CC) $(DBXEXTCFLAGS) -o dbxstatealgebra_e.o statealgebra.c 

dbxsystemalgebra_e.o: systemalgebra.c $(STDINC) $(EXTINC)
	$(CC) $(DBXEXTCFLAGS) -o dbxsystemalgebra_e.o systemalgebra.c 

dbxsystemconversion_e.o: systemconversion.c $(STDINC) $(EXTINC)
	$(CC) $(DBXEXTCFLAGS) -o dbxsystemconversion_e.o systemconversion.c 

dbxtimestep_e.o: timestep.c $(STDINC) $(EXTINC)
	$(CC) $(DBXEXTCFLAGS) -o dbxtimestep_e.o timestep.c

dbxtransformext_e.o: transformext.c $(STDINC) $(EXTINC)
	$(CC) $(DBXEXTCFLAGS) -o dbxtransformext_e.o transformext.c

#..............................................................................
