include $(NEMOLIB)/makedefs

# sun doesn't work yet, but you'd also need these :-0
SUNFLAGS = -lsocket -lnsl -lm -ldl

tcppipe:	tcppipe_nemo.c
	$(CC) $(CFLAGS) tcppipe_nemo.c -o tcppipe $(NEMOLIB)/libnemo.a -lm

clean:
	rm -f tcppipe *.o core

install:	tcppipe
	mv tcppipe $(NEMOBIN)
