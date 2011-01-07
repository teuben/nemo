# 	NEMO/src/kernel/yapp	YAPP plotting API
include $(NEMOLIB)/makedefs
LOCAL_INC = -DYAPP_$(YAPP) 
#
#	 5-Mar-89	JEB	added yapp_cg_new.c
#	22-nov-89	HOSTTYPE
#	26-oct-90	new Nemo dir.str. - put all Yapp's together, some
#			don't work though
#	9-jul-91	added untested YAPP_GT from Josh

# Standard Yapp version for library
YAPP =
# Extra YAPPO to be linked before libnemo
YAPPO = $(NEMOLIB)/yapp_ps.o
# Usual extra libs which can be used in linking
EL =

L = $(NEMOLIB)/libnemo.a
.PRECIOUS: $(L)

# old peculiar things to this directory:
BSDFORLIBS = -lF77 -lI77 -lU77	# FORLIBS itself now imported from NEMORC
LIBG = -Bdynamic -lcore -lsunwindow -lpixrect
#LIBX = -lX11 -lsocket -lnls
LIBX = $(YAPP_XLIB)

# Various predefined EL macros to easy making various yapptest_xxx's
LIBSUNVIEW  = -Bdynamic -lsuntool -lsunwindow -lpixrect
LIBSUNTOOLS = -Bdynamic -lcore -lsunwindow -lpixrect
LIBMONGO    = /usr/local/lib/libmongo.a $(LIBSUNVIEW) $(FORLIBS)
#LIBPGPLOT   = $(NEMOLIB)/libpgplot.a $(LIBSUNVIEW) -lX11 $(FORLIBS)
# linux:
LIBPGPLOT   = $(NEMOLIB)/libpgplot.a $(LIBX) $(FORLIBS)
#LIBPLPLOT   = -lplplotftk -ltk -litcl -ltcl $(LIBX) -lieee -lm -lf2c
LIBPLPLOT   = -lplplot -ltk -ltcl  $(LIBX)
# -lieee -lm
# LIBSM       = /saturn/teuben/sm/callable/dummy.o 
#  this is an SM without sunview, only X11 (check your sm/src/options.h)
LIBSM       =  -L$(NEMOLIB) -lplotsub -ldevices -lutils $(LIBX)
#  pick SGI real thing, or the VOGL emulator 
#LIBGL       = -lgl  $(LIBX)
LIBGL       = -lvogl  $(LIBX)
#LIBMESA	    = -lMesaGL $(LIBX)
LIBMESA	    = -lGL $(LIBX)

# Files exported to $NEMOLIB 

OBJFILES = yapp.o
YAPPOBJS = yapp_ps.o yapp_null.o yapp_pgplot.o yapp_plplot.o \
	yapp_cpgplot.o yapp_meta.o
BINFILES = catps catpgps yapp
# wow: the following are TRUE original NEMO V1 programs !!!    can't be in BINFILES yet
# pltext_ps plarrow_ps


TESTFILES = yapptest \
	   yapptest_mongo yapptest_pgplot yapptest_core yapptest_suntools \
	   yapptest_gl \
	   yapptest_ps yapptest_ps_new yapptest_psg \
	   yapptest_cg yapptest_cg_new \
	   yapptest_sv yapptest_sm yapptest_gt yapptest_plplot \
	   yapptest_x

# Target provided for helpless people.

help:
	@echo "YAPP: "
	@echo "Important macros used: "
	@echo "YAPP  = $(YAPP)"
	@echo "YAPPO = $(YAPPO)"
	@echo ""
	@echo object files exported are:
	@echo $(OBJFILES) $(YAPPOBJS)
	@echo executables:
	@echo $(BINFILES)
	@echo test:
	@echo $(TESTFILES)

# Targets used by Nemo during re-build.

nemo_lib: $(OBJFILES)
	mv $? $(NEMOLIB)

nemo_bin: $(BINFILES)
	mv $? $(NEMOBIN)

bins: $(BINFILES)
	mv $(BINFILES) $(NEMOBIN)

nemo_src:
	-@for i in $(BINFILES); do \
	echo `pwd` $$i ; done 

clean:
	rm -f $(OBJFILES) $(BINFILES) $(TESTFILES) *.o *.a core

cleanlib:

tests: $(TESTFILES)

# Targets used by Nemo when mantaining this directory.

install: .install_yapp

#	We will not, for now, allow yapp.o to be in $(L)
.install_lib: $(OBJFILES)
	ar ruv $(L) $?
	$(RANLIB) $(L)
	@touch .install_lib

.install_yapp: $(YAPPOBJS)
	cp $? $(NEMOLIB)

.install_bin: $(BINFILES)
	cp $? $(NEMOBIN)
	@touch .install_bin

# Rules for This Yapp Directory (now in makedefs)
#.c.o:
#	$(CC) $(CFLAGS) -DYAPP_$(YAPP) -c $<

# Dependencies.

yapp_cg.o: yapp_cg.c my_get_view_surface.c

yapp_cg_new.o: yapp_cg_new.c my_get_view_surface.c

yapp_suntools.o: yapp_suntools.c
	$(CC) $(CFLAGS) -DMOUSE -c yapp_suntools.c

# Tests.

#OK
yapptest: yapp.c
	$(CC) $(CFLAGS) -DYAPP_$(YAPP) -o yapptest yapp.c -DTESTBED \
		$(NEMO_LIBS) $(EL) -lm
#OK
yapptest_cg: yapp_cg.c
	$(CC) $(CFLAGS) -o yapptest_cg yapp_cg.c -DTESTBED \
		$(NEMO_LIBS) $(LIBSUNTOOLS) $(EL) -lm
#OK
yapptest_mongo: yapp_mongo.c testyapp.c
	$(CC) $(CFLAGS) -o yapptest_mongo yapp_mongo.c testyapp.c \
		$(NEMO_LIBS) $(LIBMONGO) $(EL) -lm
#OK
yapptest_suntools: yapp_suntools.c testyapp.c
	$(CC) $(CFLAGS) -o yapptest_suntools yapp_suntools.c testyapp.c \
		$(NEMO_LIBS) $(LIBSUNTOOLS) $(EL) -lm
#OK
yapptest_pgplot: yapp_pgplot.c testyapp.c
	$(CC) $(CFLAGS) -o yapptest_pgplot yapp_pgplot.c testyapp.c \
		$(NEMO_LIBS) $(LIBPGPLOT) $(EL) -lm

yapptest_plplot: yapp_plplot.c testyapp.c
	$(CC) $(CFLAGS) -o yapptest_plplot yapp_plplot.c testyapp.c \
		$(NEMO_LIBS) $(LIBPLPLOT) $(EL) -lm
#OK???
yapptest_cg_new: yapp_cg.c testyapp.c
	$(CC) $(CFLAGS) -o yapptest_cg_new yapp_cg_new.c -DTESTBED \
	   	$(NEMO_LIBS) $(LIBSUNTOOLS) $(EL) -lm
#OK
yapptest_core: yapp_core.c testyapp.c
	$(CC) $(CFLAGS) -o yapptest_core yapp_core.c testyapp.c \
	   	$(NEMO_LIBS) $(LIBSUNTOOLS) $(EL) -lm
#OK
yapptest_ps: yapp_ps.c testyapp.c
	$(CC) $(CFLAGS) -o yapptest_ps yapp_ps.c testyapp.c \
		$(NEMO_LIBS) $(EL) -lm

#OK
yapptest_psg: yapp_psg.c testyapp.c
	$(CC) $(CFLAGS) -o yapptest_psg yapp_psg.c testyapp.c \
		$(NEMO_LIBS) $(EL) -lm

#OK
yapptest_x: testyapp.c
	$(CC) $(CFLAGS) -o yapptest_x testyapp.c \
		$(NEMO_LIBS) -lyapp_x $(YAPP_X11) $(EL) -lm

yapptest_sv: yapp_sv.c testyapp.c
	$(CC) $(CFLAGS) -o yapptest_sv yapp_sv.c testyapp.c \
		$(NEMO_LIBS) $(LIBSUNVIEW) -lm

yapptest_gl: yapp_gl.c testyapp.c
	$(CC) $(CFLAGS) -o yapptest_gl yapp_gl.c testyapp.c \
		$(NEMO_LIBS) $(LIBGL) -lm

yapptest_mesa: yapp_mesa.c testyapp.c
	$(CC) -DMESA $(CFLAGS) -o yapptest_mesa yapp_mesa.c testyapp.c \
		$(NEMO_LIBS) $(LIBMESA) -lm

yapptest_sm: yapp_sm.c testyapp.c
	$(CC) $(CFLAGS) -o yapptest_sm yapp_sm.c testyapp.c \
		$(NEMO_LIBS) $(LIBSM) -lm

yapptest_gt: yapp_gt.c testyapp.c
	$(CC) $(CFLAGS) -o yapptest_gt yapp_gt.c testyapp.c \
		$(NEMO_LIBS) $(LIBGT) -lm

#OK
yapptest_ps_new: yapp_ps_new.c testyapp.c
	$(CC) $(CFLAGS) -o yapptest_ps_new yapp_ps_new.c testyapp.c \
		$(NEMO_LIBS) $(EL) -lm

txttst: txttst.c
	$(CC) $(CFLAGS) -o txttst txttst.c $(YAPPO) $(NEMO_LIBS) $(EL) -lm

txttst_ps: txttst.c yapp_ps.o
	$(CC) $(CFLAGS) -o txttst_ps txttst.c yapp_ps_new.o $(NEMO_LIBS) -lm

# Utilities

catps: catps.c
	$(CC) $(CFLAGS) -o catps catps.c $(NEMO_LIBS) -lm

catpgps: catpgps.c
	$(CC) $(CFLAGS) -o catpgps catpgps.c $(NEMO_LIBS) -lm

pltext: pltext.c
	$(CC) -o pltext pltext.c  $(YAPPO) $(NEMO_LIBS) $(YAPPLIB) -lm

pltext_ps: pltext.c
	$(CC) -o pltext_ps pltext.c $(NEMOLIB)/yapp_ps_new.o $(NEMO_LIBS) -lm

plarrow: plarrow.c
	$(CC) -o plarrow plarrow.c $(YAPPO) $(NEMO_LIBS) $(YAPPLIB) -lm

plarrow_ps: plarrow.c
	$(CC) -o plarrow_ps plarrow.c $(NEMOLIB)/yapp_ps_new.o $(NEMO_LIBS) -lm


yapp: yappmain.c
	$(CC) $(CFLAGS) -o yapp yappmain.c $(NEMO_LIBS) $(YAPPLIB) -lm

# the end
