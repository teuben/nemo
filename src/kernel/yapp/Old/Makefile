# $NEMOPATH/src/nemo/yapp/Makefile: Yet Another Plotting Package.
include $(NEMOLIB)/makedefs
#		5-March-89	JEB	added yapp_cg_new.c
#	22-nov HOSTTYPE
#		21-may-90	PJT	added yapp_sm
CFLAGS = -g
LDFLAGS =

LIBG = -lcore -lsunwindow -lpixrect
LIBJ = $(NEMOLIB)/libT.a $(NEMOLIB)/libJ.a

LIBSM = /saturn/teuben/sm/callable/dummy.o \
	$(NEMOLIB)/libplotsub.a $(NEMOLIB)/libdevices.a $(NEMOLIB)/libutils.a \
         -lsuntool -lsunwindow -lpixrect -lX11 

# Files exported to Nemo.

OBJFILES = yapp_ps.o yapp_ps_new.o yapp.o yapp_cg.o yapp_mongo.o \
	   yapp_suntools.o yapp_sv.o yapp_cg_new.o

BINFILES = catps pltext_ps plarrow_ps

MAN3FILES = yapp.3
MAN5FILES = yapp.5

# Target provided for helpless people.

help:
	@echo
	@echo object files exported are:
	@echo $(OBJFILES)
	@echo executables:
	@echo $(BINFILES)

# Targets used by Nemo during re-build.

nemo_lib: $(OBJFILES)
	mv $? $(NEMOLIB)

nemo_bin: $(BINFILES)
	mv $? $(NEMOBIN)

clean:
	rm -f $(OBJFILES) $(BINFILES)

# Targets used by Nemo when mantaining this directory.

install: install_man install_lib install_bin

install_lib: $(OBJFILES)
	cp $? $(NEMOLIB)
	@touch install_lib

install_bin: $(BINFILES)
	cp $? $(NEMOBIN)
	@touch install_bin

install_man: .install_man3 .install_man5

.install_man3: $(MAN3FILES) 
	cp $? $(NEMOPATH)/man/man3
	@touch .install_man3

.install_man5: $(MAN5FILES) 
	cp $? $(NEMOPATH)/man/man5
	@touch .install_man5

tidy:
	rm -f *.o yapptest* yapp.ps
#
# Dependencies.

yapp_cg.o: yapp_cg.c my_get_view_surface.c

yapp_cg_new.o: yapp_cg_new.c my_get_view_surface.c

yapp_suntools.o: yapp_suntools.c
	$(CC) $(CFLAGS) -DMOUSE -c yapp_suntools.c

# Tests.

tests: yapptest yapptest_cg yapptest_ps

yapptest: yapp.c
	$(CC) $(CFLAGS) -o yapptest yapp.c -DTESTBED $(LIBG) $(LIBJ) -lm

yapptest_cg: yapp_cg.c
	$(CC) $(CFLAGS) -o yapptest_cg yapp_cg.c -DTESTBED $(LIBG) $(LIBJ) -lm

yapptest_sm: yapp_sm.c
	$(CC) $(CFLAGS) -o yapptest_sm yapp_sm.c -DTESTBED $(LIBSM) -lm

yapptest_cg_new: yapp_cg.c
	$(CC) $(CFLAGS) -o yapptest_cg_new yapp_cg_new.c -DTESTBED \
	   $(LIBG) $(LIBJ) -lm

yapptest_ps: yapp_ps.c
	$(CC) $(CFLAGS) -o yapptest_ps yapp_ps.c -DTESTBED $(LIBJ) -lm

yapptest_ps_new: yapp_ps_new.c
	$(CC) $(CFLAGS) -o yapptest_ps_new yapp_ps_new.c -DTESTBED $(LIBJ) -lm

txttst: txttst.c yapp.o
	$(CC) -o txttst txttst.c yapp.o $(LIBJ) $(LIBG) -lm

txttst_ps: txttst.c yapp_ps.o
	$(CC) -o txttst_ps txttst.c yapp_ps_new.o $(LIBJ) -lm

# Utilities

catps: catps.c
	$(CC) -o catps catps.c $(LIBJ)

pltext: pltext.c
	$(CC) -o pltext pltext.c $(NEMOLIB)/yapp.o $(LIBJ) $(LIBG) -lm

pltext_ps: pltext.c
	$(CC) -o pltext_ps pltext.c $(NEMOLIB)/yapp_ps_new.o $(LIBJ) -lm

plarrow: plarrow.c
	$(CC) -o plarrow plarrow.c $(NEMOLIB)/yapp.o $(LIBJ) $(LIBG) -lm

plarrow_ps: plarrow.c
	$(CC) -o plarrow_ps plarrow.c $(NEMOLIB)/yapp_ps_new.o $(LIBJ) -lm

# the end

yapptest_pg:	yapp_pgplot.c
	$(CC) $(CFLAGS) -DTESTBED -o yapptest_pg yapp_pgplot.c \
	  $(LIBJ) $(MIRLIB)/libpgnew.a \
	  -lsuntool -lcgi77 -lcgi -lsunwindow -lpixrect \
	 -lF77 -lI77 -lU77 -lm
