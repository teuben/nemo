#	NEMO/src/nbody/io_nemo   - for NEMO V3
include $(NEMOLIB)/makedefs

#--------------------------------------------------------
# IO_NEMO Makefile
#--------------------------------------------------------

TARGET_BIN = snapmask_s snapmask_d io_nemo_test

#--------------------------------------------------------
# Chemin de librairies
LIBPATH = -L$(NEMOLIB) -L.


L = $(NEMOLIB)/libnemo.a

#--------------------------------------------------------
# Librairies
# NEMO
LIBNEMO = -lnemo

OBJFILES = io_nemo.o parameters.o check_file.o io_init.o \
	    io_get_put.o toolsnemo.o

#--------------------------------------------------------

lib : libio_nemo.a

clean :
	rm -f *.o core

cleanall : 
	/bin/rm -f *.o $(TARGET_BIN) *.a plum*

test  : $(TARGET_BIN)

install:	.install_lib

.install_lib: $(OBJFILES)
	ar ruv $(L) $?
	ranlib $(L)


# dependencies

io_nemo.o    : maxio.h
parameters.o : flags_data.h
check_file.o : file_data.h maxio.h
io_init.o    : flags_data.h file_data.h maxio.h
io_get_put.o : snapshot_data.h flags_data.h file_data.h maxio.h


#--------------------------------------------------------
snapmask_s : snapmask_s.o
	$(CC)  -o $@ $@.o -L$(NEMOLIB) -lnemo -lm

snapmask_s.o : snapmask2.c
	$(CC) $(CFLAGS) -DSINGLEPREC -c snapmask2.c
	mv snapmask2.o snapmask_s.o

snapmask_d : snapmask_d.o
	$(CC)  -o $@ $@.o -L$(NEMOLIB) $(NEMO_LIBS)

snapmask_d.o : snapmask2.c
	$(CC) $(CFLAGS) -c snapmask2.c
	mv snapmask2.o snapmask_d.o
#--------------------------------------------------------
io_nemo_test : io_nemo_test.o 
	$(CC) -o $@ $@.o $(NEMO_LIBS)


#
