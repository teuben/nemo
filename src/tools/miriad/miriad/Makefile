#  Makefile for 'miriad' shell
include $(NEMOLIB)/makedefs

#  No special CC is needed for miriad shell yet
CC = /bin/cc
CFLAGS = -g
#  assign $(TCL)  and  $(LIBTCL) to something sensible to get TCL
#
TCL = -DTCL -L$(MIRLIB) -I$(MIRINC) 
LIBTCL = -ltcl
#  GNU READLINE stuff
GNUFLAGS = -DREADLINE
GNULIB   = -L$(NEMOLIB) -lreadline -ltermcap

#	Miriad src install thingos
SUBMIT = miriad.c
SUBMITDIR = $(MIR)/src/tools
#	Nemo src install thingos
SRC = Makefile miriad.c README
SRCDIR = $(NEMO)/src/tools/miriad/miriad

help:
	@echo Makefile for Miriad shell 
	@echo ""
	@echo "***FLAGS which can be be used, and their defaults: "
	@echo TCL=$(TCL)
	@echo LIBTCL=$(LIBTCL)
	@echo CFLAGS=$(CFLAGS)
	@echo SUBMIT=$(SUBMIT)
	@echo ""
	@echo Targets:
	@echo "  miriad		compile a clean miriad, without TCL"
	@echo "  tcl		compile miriad with TCL assumed"
	@echo "  cleanup"
	@echo "  clean"
	@echo "  install	copy stuff into Nemo"
	@echo "  submit         submit stuff into Miriad"


miriad:	miriad.c
	$(CC) $(CFLAGS) -o miriad miriad.c

tcl:	miriad.c
	$(CC) $(TCL) $(CFLAGS) -o miriad miriad.c $(LIBTCL)

gnu:	miriad.c
	$(CC) $(GNUFLAGS) -o miriad miriad.c $(GNULIB)

clean:  cleanup

cleanup:
	rm -f flint.log core miriad *.o *.a

install: .install_src .install_bin 

.install_src: $(SRC)
	cp $? $(SRCDIR)
	@touch .install_src

.install_bin: miriad
	cp $? $(NEMOBIN)
	@touch .install_bin

#	Submit for miriad
submit: .submit

.submit: $(SUBMIT)
	mir.submit $? $(SUBMITDIR)
	@touch .submit
