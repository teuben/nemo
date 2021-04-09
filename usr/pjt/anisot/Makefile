#	local Makefile for anisot
#
include $(NEMOLIB)/makedefs

#CFLAGS = -g

L = $(NEMOLIB)/libnemo.a

#LIBG =  -lcore -lsunwindow -lpixrect
#LIBM =  $(MONGOPATH)/libmongo.a -lsuntool -lsunwindow -lpixrect
#LIBCGI = -lcgi -lsunwindow -lpixrect

#YAPP = $(NEMOPATH)/lib/yapp.o $(LIBG)
#YAPP = $(NEMOPATH)/lib/yapp_mongo.o $(LIBM) -lF77 -lI77 -lU77 -lm
#YAPPPS =  $(NEMOPATH)/lib/yapp_ps.o
#YAPPCGI =  $(NEMOPATH)/lib/yapp_cgi.o $(LIBCGI)

# what's in this directory, or supposed to be there

BIN =  anisot

SRC =  anisot.c

HDR =

MAN1 = anisot.1

MAN3 = 

MAN5 =

# and where should the sources go (directory is created if not present)

SRCDIR = $(NEMOPATH)/usr/pjt/nbody/


help:
	@echo Makefile for $(BIN)
	@echo
	@echo make nemo: installs all binaries and removes them and other junk
	@echo make all: makes all binaries
	@echo make install: installs new binaries and src, man, and Makefile
tidy:
	rm -f *.o test* tmp*

clean:
	rm -f *.o test* tmp* $(BIN)

#	install by nemo (called through nemo's master Makefile)

nemo_lib:
	@echo no nemo_lib here

nemo_bin:	$(BIN)
	mv $(BIN) $(NEMOBIN)

#	install options: by user

install: install_bin install_src install_man

install_bin: $(BIN)
	cp $? $(NEMOBIN)
	@touch install_bin

install_man: $(MAN1)
	cp $? $(NEMO)/man/man1
	@touch install_man

install_src: $(SRC) Makefile
	-@if [ ! -d $(SRCDIR) ]; \
 	then \
                echo creating $(SRCDIR); \
 		mkdir $(SRCDIR); \
 		chmod 755 $(SRCDIR); \
 	fi
	cp $? $(SRCDIR)
	@touch install_src

# now come the individual programs in this directory

a:	anisot.c
	$(CC) $(CFLAGS) -o a anisot.c $(L) -lm
