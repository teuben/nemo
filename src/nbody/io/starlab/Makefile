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


CFLAGS = -O

BIN = stod dtos
OBJ = stod_subs.o nemomain.o 
INC = 
SRC = stod_subs.c nemomain.C stod.C $(INC) Makefile
SRCDIR = $(NEMO)/src/nbody/io/starlab

LIBNEMO = -L$(NEMOLIB) -lnemo
LIBSTAR = -L$(STARLAB_PATH)/lib  -ldyn -lnode -lstd 

.SUFFIXES: .C

#	nemo:
.c.o:
	$(CC) $(CFLAGS) -c $*.c

# 	starlab:
.C.o:
	$(STARLAB_CCC) $(STARLAB_CCC_FLAGS) $(CFLAGS) \
                        -I$(NEMO)/inc  -I$(NEMOLIB) -I$(STARLAB_PATH)/inc -c $*.C

clean:
	rm -f $(BIN) $(OBJ) *.o core
	
all:	$(BIN)
	
stod:	$(OBJ) stod.o
	$(STARLAB_CCC) $(STARLAB_CCC_FLAGS) $(CFLAGS) \
               -o stod $(OBJ) stod.o \
		$(LIBNEMO) $(LIBSTAR) -lm

dtos:	$(OBJ) dtos.o
	$(STARLAB_CCC) -o dtos $(OBJ) dtos.o \
		$(LIBNEMO) $(LIBSTAR) -lm

pjt:	.pjt

.pjt:	$(SRC)
	cp $? $(SRCDIR)
	touch .pjt

submit:
