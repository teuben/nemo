# supposed to work with STARLAB 2.1, but only in real=double mode
include $(NEMOLIB)/makedefs
# does not work in real=float

help:
	@echo These are programs that convert data between NEMO and STARLAB
	@echo The Makefile assumes the STARLAB_* variables are present,
	@echo most notably:
	@echo " STARLAB_PATH=$(STARLAB_PATH)"
	@echo " STARLAB_CCC=$(STARLAB_CCC) "
	@echo " STARLAB_CCC_FLAGS=$(STARLAB_CCC_FLAGS) "
	@echo "If any of these are not present, STARLAB is probably not loaded"


#CFLAGS = -O

BIN = stod dtos
OBJ = stod_subs.o nemomain.o 
INC = 
SRC = stod_subs.c nemomain.C stod.C $(INC) Makefile
SRCDIR = $(NEMO)/src/nbody/io/starlab

LIBNEMO = -L$(NEMOLIB) -lnemo
#LIBSTAR = -L$(STARLAB_PATH)/lib  -ltdyn -ldyn -lnode -lstd -lsstar
LIBSTAR = -L$(STARLAB_INSTALL_PATH)/lib/starlab  -ltdyn -ldyn -lnode -lstd 

STARLAB_CCC_FLAGS =  -I$(STARLAB_INSTALL_PATH)/include/starlab

.SUFFIXES: .C

#	nemo:
#.c.o:
#	$(CC) $(CFLAGS) -c $*.c

# 	starlab:
.C.o:
	$(CXX) $(STARLAB_CCC_FLAGS) $(CFLAGS) \
		 -c $*.C

nemomain.o: nemomain.C
	$(CXX) $(CFLAGS) -c $*.C


#                 -I$(NEMO)/inc  -I$(NEMOLIB) 

clean:
	rm -f $(BIN) $(OBJ) *.o core

all:	$(BIN)

stod:	symlink $(OBJ) stod.o
	$(CXX) $(STARLAB_CCC_FLAGS) $(CFLAGS) \
               -o stod $(OBJ) stod.o \
		$(LIBNEMO) $(LIBSTAR) -lm

dtos:	symlink $(OBJ) dtos.o
	$(CXX) -o dtos $(OBJ) dtos.o \
		$(LIBNEMO) $(LIBSTAR) -lm

symlink:  $(NEMOINC)/nemo_stdinc.h

$(NEMOINC)/nemo_stdinc.h:
	@if [ ! -e $(NEMOINC)/nemo_stdinc.h ]; \
	then \
	  echo Create symlink for nemo_stdinc.h for starlab+nemo;\
	  ln -s $(NEMOINC)/stdinc.h $(NEMOINC)/nemo_stdinc.h;\
	else \
	  echo OK; \
	fi

pjt:	.pjt

.pjt:	$(SRC)
	cp $? $(SRCDIR)
	touch .pjt

submit:
