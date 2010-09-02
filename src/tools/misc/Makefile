# 	$NEMO/src/tools/misc
#
#	dec-90	moved small things from $NEMO/src/nemo/util into here

include $(NEMOLIB)/makedefs

#
CFLAGS = -O
LDFLAGS = -O

#UTILITIES = crc hd changed redir age  accudate
UTILITIES = hd redir accudate

help:
	@echo "misc NEMO stuff ($(UTILITIES))"
	@echo No more specific help here

# old entries
nemo_lib:
	@echo no nemo_lib here

nemo_bin: $(UTILITIES)
	mv $? $(NEMOBIN)
	rm -f *.o *.a core $(UTILITIES)

bins: $(UTILITIES)
	mv $? $(NEMOBIN)

nemo_src:
	-@for i in $(UTILITIES) ; do \
	echo `pwd` $$i ; done 
		

all:	$(UTILITIES)
		
clean:
	rm -f $(UTILITIES) *.o *.a core 

####				install by USER (nemo)

install: install_bin
	@echo installing utilities: $(UTILITIES) ....

install_bin: $(UTILITIES)
	mv $? $(NEMOBIN)

tidy:
	rm -f *.o *.BAK *.CKP core

crc: crc.c
	$(CC) $(CFLAGS) -o crc crc.c 

hd: hd.c
	$(CC) $(CFLAGS) -o hd hd.c

changed:	changed.c
	$(CC) $(CFLAGS) -o changed changed.c

redir:	redir.c
	$(CC) $(CFLAGS) -o redir redir.c

age:	age.c
	$(CC) $(CFLAGS) -o age age.c

accudate:	accudate.c
	$(CC) $(CFLAGS) -o accudate accudate.c


