prefix = /usr/local
exec_prefix = $(prefix)
bindir = $(exec_prefix)/bin

CC = gcc
INSTALL = install
CFLAGS = -Wformat
INSTALLFLAGS = -s

programs = tcppipe

.PHONY: all clean install uninstall ed

all: $(programs)

clean:
	rm -f $(programs)

install: $(programs)
	$(INSTALL) $(INSTALLFLAGS) $(programs) $(bindir)

uninstall:
	for prog in $(programs); do rm -f $(bindir)/$$prog; done

ed:
	vi tcppipe.c
	make
