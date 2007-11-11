# -*- makefile -*-                                                             
#-----------------------------------------------------------------------------+
#                                                                             |
# Makefile                                                                    |
#                                                                             |
#                                                                             |
# Copyright Jean-Charles LAMBERT - 2006                                       |
# e-mail:   Jean-Charles.Lambert@oamp.fr                                      |
# address:  Dynamique des galaxies                                            |
#           Laboratoire d'Astrophysique de Marseille                          |
#           2, place Le Verrier                                               |
#           13248 Marseille Cedex 4, France                                   |
#           CNRS U.M.R 6110                                                   |
#                                                                             |
#-----------------------------------------------------------------------------+

#------------------------------------------------------------------------------
# gyrfalcon location and options                                               

# check if proprietary version of the code exist
PROP := $(shell if test -d ${NEMO}/usr/dehnen/falcON.P; then echo 1; fi)

# according to Proprietary and Public version of falcON      
# you must populate GYRFALCON with falcON compiling tree path
# and the rights compilation flags                           
ifeq ($(PROP),1) 
    VERSION        := "PROPRIETARY"
    GYRFALCON      := $(NEMO)/usr/dehnen/falcON.P
    FULLPROPER     := -DfalcON_INDI -DfalcON_NEMO -DfalcON_PROPER -DfalcON_SSE -DNEMO
else
    VERSION        := "PUBLIC"
    GYRFALCON      := $(NEMO)/usr/dehnen/falcON
    FULLPROPER     := -DfalcON_INDI -DfalcON_NEMO
endif

ifeq ($(GYRFALCON),"")
$(error variable GYRFALCON inside Makefile must be populate with falcON compiling tree path)
endif
#------------------------------------------------------------------------------
# macros: directories & names                                                  
#
SRC            := src/
OBJ            := obj/
MANIP          := shared-lib/
INC            := -I./ -I./inc -I$(GYRFALCON)/inc -I$(GYRFALCON)/utils/inc -I$(GYRFALCON)/../utils/inc -I$(NEMOINC) -I$(NEMOLIB)
MNAME          := glnemo

#------------------------------------------------------------------------------
# C compiler                                                                   
FC  = g77
CPP = /usr/bin/g++-3.3.6
CC  = gcc
CPPFLAGS = -g $(FULLPROPER) -rdynamic

NET_ENGINE = $(OBJ)MessageBuffer.o  $(OBJ)ChunkBuffer.o          $(OBJ)Socket.o        $(OBJ)SocketServer.o
MTS_ENGINE = $(OBJ)generic_thread.o $(OBJ)master_server_thread.o $(OBJ)server_thread.o
GLNOBJ     = $(NET_ENGINE) $(MTS_ENGINE) $(OBJ)glnemo.o
#------------------------------------------------------------------------------
# list of manip targets                                                        
ommanips: dirs                        \
	  $(MANIP)glnemo.so

#------------------------------------------------------------------------------
# manip compiling rules                                                        

$(MANIP)glnemo.so    :	$(GLNOBJ)
			@echo "****"
			@echo "**** Linking " $(VERSION) " version of the manipulator"
			@echo "****"
			$(CPP) $(GLNOBJ) -o $@ $(CPPFLAGS)  $(INC) -lpthread -shared

$(OBJ)glnemo.o        :	$(SRC)glnemo.cc $(SRC)master_server_thread.h 
			$(CPP) -fPIC -c $< -o $@ $(CPPFLAGS)  $(INC)

clean:
	\rm -rf obj $(MANIP)glnemo.so
## DEPENDENCIES

#------------------------------------------------------------------------------
# Network Engine                                                               
#
$(OBJ)MessageBuffer.o : $(SRC)MessageBuffer.cc $(SRC)MessageBuffer.h $(SRC)ChunkBuffer.h
	$(CPP) -fPIC $(INC) -c $(SRC)MessageBuffer.cc -o $(OBJ)MessageBuffer.o

$(OBJ)ChunkBuffer.o : $(SRC)ChunkBuffer.cc $(SRC)ChunkBuffer.h
	$(CPP) -fPIC $(INC) -c $(SRC)ChunkBuffer.cc -o $(OBJ)ChunkBuffer.o

$(OBJ)Socket.o : $(SRC)Socket.h $(SRC)Socket.cc
	$(CPP) -fPIC $(INC) -c $(SRC)Socket.cc -o $(OBJ)Socket.o

$(OBJ)SocketServer.o : $(SRC)Socket.h $(SRC)SocketServer.h $(SRC)SocketServer.cc
	$(CPP) -fPIC $(INC) -c $(SRC)SocketServer.cc -o $(OBJ)SocketServer.o

#------------------------------------------------------------------------------
# Multi Threaded Server Engine                                                 
#
$(OBJ)generic_thread.o :                                                       \
		$(SRC)generic_thread.h        $(SRC)generic_thread.cc
	$(CPP) -fPIC $(INC)                                                    \
		-c $(SRC)generic_thread.cc -o  $(OBJ)generic_thread.o

$(OBJ)master_server_thread.o :                                                 \
		$(SRC)MessageBuffer.h         $(SRC)SocketServer.h             \
		$(SRC)generic_thread.h        $(SRC)server_thread.h            \
		$(SRC)master_server_thread.h  $(SRC)master_server_thread.cc
	$(CPP) -fPIC $(INC)  $(NBDYFLAGS) $(CPPFLAGS)                          \
		-c $(SRC)master_server_thread.cc -o $(OBJ)master_server_thread.o

$(OBJ)server_thread.o : 	  					       \
		$(SRC)MessageBuffer.h         $(SRC)generic_thread.h           \
		$(SRC)server_thread.h         $(SRC)server_thread.cc
	$(CPP) -fPIC $(INC) $(NBDYFLAGS) $(CPPFLAGS)                           \
		-c $(SRC)server_thread.cc -o $(OBJ)server_thread.o

#------------------------------------------------------------------------------
# dirs                                                                         

dirs:
	@mkdir -p $(MANIP) $(OBJ)

#------------------------------------------------------------------------------
# package                                                                      
tar:
	(cd ..; tar czhvf manipulators/$(MNAME).tar.gz `ls manipulators/*.{cc,h,txt} manipulators/Makefile`)

# install
install: ommanips  #$(MANIP)$(MNAME).so
	@echo "****"
	@echo "**** Installing " $(VERSION) " version of the manipulator"
	@echo "**** here [" $(GYRFALCON)/manip "] "
	@echo "****"
	mkdir -p $(GYRFALCON)/manip
	\cp -p $(MANIP)$(MNAME).so $(GYRFALCON)/manip

#
