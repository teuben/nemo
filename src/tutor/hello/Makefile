#
include $(NEMOLIB)/makedefs
#   Example Makefile for hello.c within NEMO environment
#CFLAGS = -g
#LDFLAGS = -g
BIN = hello hello++ bye hallo
SRC = hello.c hallo.f
MAN = hello.1
CCC = g++ -I$(NEMOINC) -L$(NEMOLIB)
SRCDIR =  $(NEMO)/src/tutor/hello
L = $(NEMOLIB)/libnemo.a

#	catch default
help:
	@echo Example Makefile for 'hello'

#	required nemo install targets (for system installation)
nemo_lib:
	@echo Directory `pwd`: no nemo_lib 

nemo_bin:	$(BIN)
	mv $? $(NEMOBIN)

nemo_src:  $(BIN)
	-@do i in $(BIN) ;\
	(echo `pwd` $$i);\
	enddo

#	default rules

.cc:
	$(CXX) $(CFLAGS) -o $@ $*.cc -L$(NEMOLIB) -lnemo++ $(L) $(EL) -lm

#	optional user install targets (for user installation)
install:	install_src install_bin install_man

install_src:	$(SRC)
	cp $? $(SRCDIR)
	@touch install_src

install_bin:	$(BIN)
	cp $? $(NEMOBIN)
	@touch install_bin

install_man:	$(MAN)
	cp $? $(NEMO)/man/man1
	@touch install_man

clean:
	rm -f *.o *.a core $(BIN)
#	non-standard worker targets
