

ULIB = $(NEMOLIB)
UINC = $(NEMOINC)

CC =		cc
CFLAGS =	-O4 
STATICFLAGS =	-Bstatic

# CC = gcc
# CFLAGS = 	-fpic -O2 -ansi -traditional \
# 		-W -Wunused -Wcomment -Wtrigraphs -Wformat -Wuninitialized \
# 		-Wparentheses -Wshadow -Wpointer-arith -Wcast-qual \
# 		-Wconversion -Wswitch -Wid-clash-32
# STATICFLAGS =	-static


LDFLAGS   = 	# -assert pure-text
LINTFLAGS = 	-abhuxz

# don't insert a line break after the -d
SHARFLAGS =	-v -l40 -M -d__KRILL_ARE_YUMMY_AND_CRUNCHY

VERSION =	0.1
SHELL =		/bin/sh

OBJS = 	merge_sort.o \
	quick_sort.o \
	heap_sort.o \
	shell_sort.o \
	bubble_sort.o \
	bogo_sort.o \
	insert_sort.o

flogger: flogger.o libsort.a
	$(CC) $(CFLAGS) $(STATICFLAGS) -o flogger flogger.o -L. -lsort

flogger.dynamic: flogger.o libsort.so.$(VERSION)
	$(CC) $(CFLAGS) -o flogger.dynamic flogger.o -L. -lsort

libsort.a: $(OBJS)
	ar ru libsort.a $(OBJS)
	ranlib libsort.a

libsort.so.$(VERSION): $(OBJS)
	#ld $(LDFLAGS) $(OBJS) -o libsort.so.$(VERSION)

$(OBJS): sorting.h

install: libsort.a sorting.h
	cp libsort.a $(ULIB)
	cp sorting.h $(UINC)


install_dyn:
	cp libsort.so.$(VERSION) $(ULIB)
	chmod 755 $(ULIB)/libsort.so.$(VERSION)
	#/usr/etc/ldconfig

lint:	lint.out

lint.out: *.c sorting.h
	lint $(LINTFLAGS) *.c > lint.out 2>&1

clean:
	-rm -f core shar0? $(OBJS) \
	libsort.so.$(VERSION) libsort.a flogger.o \
	flogger flogger.static flogger.dynamic lint.out

shar:	
	shar $(SHARFLAGS) -o shar \
	README INSTALLATION TODO OPTIMIZING Makefile bubble_sort.c \
	heap_sort.c merge_sort.c quick_sort.c shell_sort.c insert_sort.c \
	sorting.h flogger.c bogo_sort.c

