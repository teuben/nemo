#	config utilities
include $(NEMOLIB)/makedefs

help:
	more config.c

config:	config.c
	$(CC) $(CFLAGS) -o config config.c

clean:
	rm -f config core *.o

