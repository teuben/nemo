
include $(NEMOLIB)/makedefs

#SRC = fitsverify.c ftverify.c fvrf_data.c fvrf_file.c fvrf_head.c fvrf_key.c fvrf_misc.c
SRC = ftverify.c fvrf_data.c fvrf_file.c fvrf_head.c fvrf_key.c fvrf_misc.c


BINS = 	fitsverify


help:
	
all:	$(BINS)

bins: 	all
	mv $(BINS) $(NEMOBIN)
	
fitsverify:	$(SRC)
	$(CC) $(CFLAGS) -DSTANDALONE -o fitsverify $(SRC) $(CFITSIO_INC) $(CFITSIO_LIB) -lm
#	$(CC) $(CFLAGS) -DSTANDALONE -o fitsverify $(SRC) -I$(NEMOINC) -L$(NEMOLIB) -lcfitsio -lm

clean:
	rm -f *.o $(BINS)


