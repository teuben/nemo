# This makefile is supposed to be standalone and can be used within
# an existing NEMO distribution to extract the relevant files for
# uNEMO to be built. 
# Within NEMO it normally lives in $NEMO/src/uNEMO.
#
#	31-may-2001	Created			Peter Teuben

CFLAGS = -g


INC = config.h filefn.h	filestruct.h history.h snapshot.h  strlib.h \
	extstring.h filesecret.h getparam.h options.h stdinc.h version.h

LOCAL_SRC = getparam_fake.c simple_read.c simple_write.c

SRC = allocate.c bswap.c convert.c dprintf.c error.c extstring.c filefn.c \
	filesecret.c history.c strlib.c stropen.c date_id.c


OBJ = allocate.o bswap.o convert.o dprintf.o error.o extstring.o filefn.o \
	filesecret.o getparam_fake.o history.o strlib.o stropen.o

BIN = simple_read simple_write

LIB = libunemo.a

TARFILES = README.unemo $(INC) $(SRC) Makefile simple_read.c simple_write.c 


help:
	@echo there is no help here .... yet

clean:
	rm -f $(LIB)  *.o core  $(BIN) junk.dat

clean_all:
	@rm -f $(SRC) $(INC)

lib:	$(OBJ)
	ar ruv $(LIB) $?

bin:	$(BIN)

all:	lib bin

tar:
	tar cf unemo.tar $(TARFILES)

.c.o:
	$(CC) -I. $(CFLAGS) -c $*.c 

simple_read:
	$(CC) -I. $(CFLAGS) -o simple_read simple_read.c libunemo.a -lm

simple_write:
	$(CC) -I. $(CFLAGS) -o simple_write simple_write.c libunemo.a -lm


### Difference the current version with the NEMO version

DIFF = diff

NEMO_INC = extstring.h filefn.h filestruct.h getparam.h \
	   history.h stdinc.h options.h strlib.h version.h

NEMO_ISS = snapshot.h

NEMO_CORES = allocate.c bswap.c date_id.c error.c  filefn.c strlib.c

NEMO_IO = convert.c dprintf.c extstring.c filesecret.c filesecret.h \
	stropen.c  history.c

NEMO_LIB = config.h

diff: diff_inc diff_iss diff_cores diff_io diff_lib

diff_inc:
	-@for i in $(NEMO_INC); do\
	($(DIFF) $$i $(NEMOINC)); done

diff_iss:
	-@for i in $(NEMO_ISS); do\
	($(DIFF) $$i $(NEMOINC)/snapshot); done

diff_cores:
	-@for i in $(NEMO_CORES); do\
	($(DIFF) $$i $(NEMO)/src/kernel/cores); done

diff_io:
	-@for i in $(NEMO_IO); do\
	($(DIFF) $$i $(NEMO)/src/kernel/io); done

diff_lib:
	-@for i in $(NEMO_LIB); do\
	($(DIFF) $$i $(NEMOLIB)); done
	
####  Copy a new version from $NEMO

CP = cp

cp:	cp_inc cp_iss cp_cores cp_io cp_lib

cp_inc:
	-@for i in $(NEMO_INC); do\
	($(CP) $(NEMOINC)/$$i .); done

cp_iss:
	-@for i in $(NEMO_ISS); do\
	($(CP) $(NEMOINC)/snapshot/$$i .); done

cp_cores:
	-@for i in $(NEMO_CORES); do\
	($(CP) $(NEMO)/src/kernel/cores/$$i .); done

cp_io:
	-@for i in $(NEMO_IO); do\
	($(CP) $(NEMO)/src/kernel/io/$$i .); done

cp_lib:
	-@for i in $(NEMO_LIB); do\
	($(CP) $(NEMOLIB)/$$i .); done




