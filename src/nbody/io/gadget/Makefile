# -*-makefile-*-

#include $(NEMOLIB)/makedefs

# Compilation otions
CPP      = g++
CPPFLAGS = -I$(NEMOINC) -I$(NEMOLIB) -Wall -O2
LNEMO    = -L$(NEMOLIB)  -lio_nemo -lnemo++ -lnemo
LDL      = -ldl

# Files OBJ
FOBJ = compile/obj
# Files BIN
FBIN = compile/bin
# FILES SOURCES
FSRC = ./

# targets
bin : 	dirs $(FBIN)/gadget2nemo $(FBIN)/nemo2gadget
gad :   dirs $(FBIN)/gadget2nemo
nem :   dirs $(FBIN)/nemo2gadget
install : bin
	\cp -p $(FBIN)/gadget2nemo $(FBIN)/nemo2gadget ${NEMOBIN}

#--
#
OBJTOOLS = $(FOBJ)/gadget_endian_tools.o
$(OBJTOOLS) : $(FSRC)/gadget_endian_tools.cc $(FSRC)/gadget_endian_tools.h    \
					     $(FSRC)/gadget_data_structure.h
	$(CPP) 	$(CPPFLAGS)  -o $@ -c  $(FSRC)/gadget_endian_tools.cc
#--
OBJ1    = $(FOBJ)/gadget2nemo.o 

$(FBIN)/gadget2nemo :  $(OBJ1) $(OBJTOOLS)
	$(CPP)  -o $@ $(OBJ1) $(OBJTOOLS) $(LNEMO) ${LDL} -lstdc++ -lm 

$(OBJ1) : $(FSRC)/gadget2nemo.cc  $(FSRC)/gadget_endian_tools.h $(FSRC)/gadget_data_structure.h 
	$(CPP) 	$(CPPFLAGS)  -o $@ -c  $(FSRC)/gadget2nemo.cc
#--
OBJ2    = $(FOBJ)/nemo2gadget.o 

$(FBIN)/nemo2gadget :  $(OBJ2) $(OBJTOOLS)
	$(CPP)  -o $@ $(OBJ2) $(OBJTOOLS) $(LNEMO) ${LDL} -lstdc++ -lm 

$(OBJ2) : $(FSRC)/nemo2gadget.cc $(FSRC)/gadget_endian_tools.h $(FSRC)/gadget_data_structure.h 
	$(CPP) 	$(CPPFLAGS)  -o $@ -c  $(FSRC)/nemo2gadget.cc



# targets
dirs : 
	@mkdir -p ${FOBJ} ${FBIN}

tar:
		(cd ../..; tar czhvf gadget-tools/gadget-tools.tar.gz `ls gadget-tools/src/*.{cc,h,txt}  gadget-tools/src/Makefile`)
#
