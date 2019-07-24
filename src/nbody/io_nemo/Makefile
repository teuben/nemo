#--------------------------------------------------------
#
#	NEMO/src/nbody/io_nemo   - for NEMO V3
#
#--------------------------------------------------------

# Nemo's definitions
include $(NEMOLIB)/makedefs

#--------------------------------------------------------
# Library path
LIBPATH = -L$(NEMOLIB) -L.
L       = $(NEMOLIB)/libnemo.a 
Lg77    = $(NEMOLIB)/libnemomaing77.a 
Ldummyc = $(NEMOLIB)/libio_nemo.a
Ldummyf = $(NEMOLIB)/libio_nemo_f.a

LIBNEMO = -lnemo
#--------------------------------------------------------
# compilers options
FC     = gfortran

ifeq ($(FC),gfortran)
FFLAGS = -O -fno-second-underscore
DFC    = GFORT
else
FFLAGS = -O -fno-second-underscore -Wno-globals
DFC    = G77
endif
OPT    = -Wall -g
#--------------------------------------------------------
# Absoft f77's options
#FC=f77
#FFLAGS = -O -B108 -K

#--------------------------------------------------------
# Macros definitions
# DESTDIR
DESTDIR = compile
# io_nemo OBJ
IONO = ${DESTDIR}/obj
# io_nemo BIN
IONB = ${DESTDIR}/bin
# io_nemo Test Src
IONT = test_src
# io_nemo src
IONS = src

# LIB_IO_NEMO
LIB_IO_NEMO = ${IONO}/parameters.o    ${IONO}/check_file.o       ${IONO}/io_init.o \
	      ${IONO}/io_nemo_tools.o ${IONO}/get_data_wrapper.o                   \
              ${IONO}/io_nemo.o       ${IONO}/io_get_put.o                         \
              ${IONO}/io_nemo_f.o     ${IONO}/io_get_put_f.o

LIB_IO_NEMOF= ${IONO}/nemo_g77.o      $(LIB_IO_NEMO)

LIB_DUMMY_IO_NEMO = ${IONO}/dummy.o

# Includes
INC         = -I${IONT} -I${IONS}
#--------------------------------------------------------
# Main targets
#

clean:    allclean

allclean  : clean_obj clean_test

clean_obj :
	/bin/rm -f ${IONO}/*.o

clean_test :
	/bin/rm -rf compile/bin/*

run_Testfile : dummy
	@make -f Testfile clean testc testf

run_test : dummy
	(make test  ; script/test_io_nemo.csh )
	(make testf ; script/fortran_test_io_nemo.csh )

run_stress: dummy
	(make test  ; script/stress_io_nemo.pl)

test  : dirs                                                   \
	${IONB}/io_nemo_test_float ${IONB}/io_nemo_test_double \
	${IONB}/stress_io_nemo     ${IONB}/stress_io_nemo_f    \
	${IONB}/snapmask_s         ${IONB}/snapmask_d

testf : dirs 						       \
	${IONB}/nemo_fortran_f_3n  ${IONB}/nemo_fortran_d_3n   \
        ${IONB}/nemo_fortran_f_n3  ${IONB}/nemo_fortran_d_n3   \
	${IONB}/snapmask_s         ${IONB}/snapmask_d

install:  dirs	.install_lib .install_lib_maing77 .install_lib_dummy test testf man

.install_lib: $(LIB_IO_NEMO)
	ar ruv $(L) $?
	ranlib $(L)

.install_lib_maing77: ${IONO}/nemo_g77.o
	ar ruv $(Lg77) $?
	ranlib $(Lg77)

.install_lib_dummy : $(LIB_DUMMY_IO_NEMO)
	ar ruv $(Ldummyc) $?
	ar ruv $(Ldummyf) $?

lib: dirs $(LIB_IO_NEMO)
bins:
#	@echo No bins in nemo_io

tar :
	        name=io_nemo-`date "+%d%b-%H.%M"`; cd ..; tar czhvf io_nemo/store/$${name}.tar.gz `find io_nemo/ -type f \( -name "*.[f,F,c,h,3]"  -o -name "*.cc" -o -name "Makefile"  -o -name "*.ps"  -o -name "README*" -o -name "*.csh" -o -name "*.pl" -o -name "Testfile" \)`; cp io_nemo/store/$${name}.tar.gz ~/io_nemo.tar.gz
dirs : 
	@script/make_dirs.csh

doc  : docps docascii

docps : 
	@groff -Tps man/io_nemo.3   > doc/io_nemo.ps
	@groff -Tps man/io_nemo_f.3 > doc/io_nemo_f.ps
docascii : 
	@groff -Tascii man/io_nemo.3   > doc/io_nemo.ascii
	@groff -Tascii man/io_nemo_f.3 > doc/io_nemo_f.ascii

man: dummy
	/bin/cp man/io_nemo.3 man/io_nemo_f.3 $(NEMO)/man/man3

dummy :
#--------------------------------------------------------
# dependencies

# Specific io_nemo's files :
${IONO}/io_nemo.o          : ${IONS}/io_nemo.h      ${IONS}/maxio.h          \
			     ${IONS}/io_init.h      ${IONS}/parameters.h     \
			     ${IONS}/io_get_put.h   ${IONS}/check_file.h     \
			     ${IONS}/io_nemo.c
	$(CC) $(CFLAGS) $(OPT) $(INC) -c ${IONS}/io_nemo.c -o ${IONO}/io_nemo.o

${IONO}/io_get_put.o       : ${IONS}/io_get_put.h   ${IONS}/io_nemo_data.h   \
			     ${IONS}/flags_data.h   ${IONS}/history_data.h   \
			     ${IONS}/check_file.h   ${IONS}/io_nemo_tools.h  \
			     ${IONS}/get_data_wrapper.h                      \
			     ${IONS}/io_get_put.c
	$(CC) $(CFLAGS) $(OPT) $(INC) -c ${IONS}/io_get_put.c -o ${IONO}/io_get_put.o

# Specific io_nemo_f's files :
${IONO}/io_nemo_f.o        : ${IONS}/io_nemo_f.h    ${IONS}/maxio.h          \
			     ${IONS}/io_init.h      ${IONS}/parameters.h     \
			     ${IONS}/io_get_put_f.h ${IONS}/check_file.h     \
			     ${IONS}/io_nemo_f.c
	$(CC) $(CFLAGS) $(OPT) $(INC) -c ${IONS}/io_nemo_f.c -o ${IONO}/io_nemo_f.o

${IONO}/io_get_put_f.o     : ${IONS}/io_get_put_f.h ${IONS}/io_nemo_data_f.h \
			     ${IONS}/flags_data.h   ${IONS}/history_data_f.h \
		     	     ${IONS}/check_file.h   ${IONS}/io_nemo_tools.h  \
			     ${IONS}/get_data_wrapper.h                      \
			     ${IONS}/io_get_put_f.c
	$(CC) $(CFLAGS) $(OPT) $(INC) -c ${IONS}/io_get_put_f.c -o ${IONO}/io_get_put_f.o


# Common C files with "io_nemo" and "io_nemo_f" :
${IONO}/io_init.o          : ${IONS}/io_init.h      ${IONS}/io_nemo_tools.h  \
			     ${IONS}/io_init.c
	$(CC) $(CFLAGS) $(OPT) $(INC) -c ${IONS}/io_init.c -o ${IONO}/io_init.o 

${IONO}/parameters.o       : ${IONS}/parameters.h   ${IONS}/flags_data.h     \
			     ${IONS}/parameters.c
	$(CC) $(CFLAGS) $(OPT) $(INC) -c ${IONS}/parameters.c -o ${IONO}/parameters.o

${IONO}/io_nemo_tools.o    : ${IONS}/io_nemo_tools.h ${IONS}/io_nemo_tools.c
	$(CC) $(CFLAGS) $(OPT) $(INC) -c ${IONS}/io_nemo_tools.c -o ${IONO}/io_nemo_tools.o

${IONO}/get_data_wrapper.o : ${IONS}/get_data_wrapper.h ${IONS}/get_data_wrapper.c
	$(CC) $(CFLAGS) $(OPT) $(INC) -c ${IONS}/get_data_wrapper.c -o ${IONO}/get_data_wrapper.o

${IONO}/check_file.o       : ${IONS}/check_file.h ${IONS}/check_file.c
	$(CC) $(CFLAGS) $(OPT) $(INC) -c ${IONS}/check_file.c -o ${IONO}/check_file.o

# main program used with g77 compiler
${IONO}/nemo_g77.o         : ${IONS}/nemo_g77.c
	$(CC) $(CFLAGS) -D$(DFC) $(OPT) $(INC) -c ${IONS}/nemo_g77.c -o ${IONO}/nemo_g77.o 

#
# io_nemo test programs
#

${IONB}/io_nemo_test_float : ${IONO}/io_nemo_test_float.o $(LIB_IO_NEMO)
	$(CC) -o ${IONB}/io_nemo_test_float ${IONO}/io_nemo_test_float.o    \
	                          $(LIB_IO_NEMO) -L$(NEMOLIB) -lnemo -lm

${IONO}/io_nemo_test_float.o: ${IONT}/io_nemo_test.c
	$(CC) $(CFLAGS) $(OPT) $(INC) -DSINGLEPREC -c ${IONT}/io_nemo_test.c       \
			       -o ${IONO}/io_nemo_test_float.o

${IONB}/io_nemo_test_double : ${IONO}/io_nemo_test_double.o $(LIB_IO_NEMO)
	$(CC) -o ${IONB}/io_nemo_test_double ${IONO}/io_nemo_test_double.o  \
	                            $(LIB_IO_NEMO) -L$(NEMOLIB) -lnemo -lm

${IONO}/io_nemo_test_double.o: ${IONT}/io_nemo_test.c
	$(CC) $(CFLAGS) $(OPT) $(INC) -c ${IONT}/io_nemo_test.c                    \
			       -o ${IONO}/io_nemo_test_double.o


${IONO}/stress_io_nemo.o: ${IONT}/stress_io_nemo.cc
	$(CXX) $(CFLAGS) $(OPT) $(INC) -DSINGLEPREC -c ${IONT}/stress_io_nemo.cc       \
			       -o ${IONO}/stress_io_nemo.o

${IONB}/stress_io_nemo : ${IONO}/stress_io_nemo.o $(LIB_IO_NEMO)
	$(CXX) -o ${IONB}/stress_io_nemo ${IONO}/stress_io_nemo.o    \
	                          $(LIB_IO_NEMO) -L$(NEMOLIB) -lnemo -lm

${IONO}/stress_io_nemo_f.o: ${IONT}/stress_io_nemo_f.cc
	$(CXX) $(CFLAGS) $(OPT) $(INC) -DSINGLEPREC -c ${IONT}/stress_io_nemo_f.cc       \
			       -o ${IONO}/stress_io_nemo_f.o

${IONB}/stress_io_nemo_f : ${IONO}/stress_io_nemo_f.o $(LIB_IO_NEMO)
	$(CXX) -o ${IONB}/stress_io_nemo_f ${IONO}/stress_io_nemo_f.o    \
	                          $(LIB_IO_NEMO) -L$(NEMOLIB) -lnemo -lm

#
# io_nemo_f test programs
#
${IONB}/nemo_fortran_f_3n  : ${IONO}/nemo_fortran_f_3n.o $(LIB_IO_NEMOF)
	$(FC) -o ${IONB}/nemo_fortran_f_3n ${IONO}/nemo_fortran_f_3n.o      \
               $(LIBPATH)    $(LIB_IO_NEMOF) -lnemo  -lm

${IONO}/nemo_fortran_f_3n.o: ${IONT}/nemo_fortran_f_3n.F
	$(FC) $(FFLAGS) -c   ${IONT}/nemo_fortran_f_3n.F                    \
                        -o   ${IONO}/nemo_fortran_f_3n.o

${IONB}/nemo_fortran_d_3n  : ${IONO}/nemo_fortran_d_3n.o $(LIB_IO_NEMOF)  
	$(FC) -o ${IONB}/nemo_fortran_d_3n ${IONO}/nemo_fortran_d_3n.o      \
               $(LIBPATH)    $(LIB_IO_NEMOF) -lnemo  -lm

${IONO}/nemo_fortran_d_3n.o: ${IONT}/nemo_fortran_d_3n.F
	$(FC) $(FFLAGS) -c   ${IONT}/nemo_fortran_d_3n.F                    \
                        -o   ${IONO}/nemo_fortran_d_3n.o

${IONB}/nemo_fortran_f_n3  : ${IONO}/nemo_fortran_f_n3.o $(LIB_IO_NEMOF)  
	$(FC) -o ${IONB}/nemo_fortran_f_n3 ${IONO}/nemo_fortran_f_n3.o      \
               $(LIBPATH)    $(LIB_IO_NEMOF) -lnemo  -lm

${IONO}/nemo_fortran_f_n3.o: ${IONT}/nemo_fortran_f_n3.F
	$(FC) $(FFLAGS) -c   ${IONT}/nemo_fortran_f_n3.F                    \
                        -o   ${IONO}/nemo_fortran_f_n3.o

${IONB}/nemo_fortran_d_n3  : ${IONO}/nemo_fortran_d_n3.o $(LIB_IO_NEMOF)  
	$(FC) -o ${IONB}/nemo_fortran_d_n3 ${IONO}/nemo_fortran_d_n3.o      \
               $(LIBPATH)    $(LIB_IO_NEMOF) -lnemo  -lm

${IONO}/nemo_fortran_d_n3.o: ${IONT}/nemo_fortran_d_n3.F
	$(FC) $(FFLAGS) -c   ${IONT}/nemo_fortran_d_n3.F                    \
                        -o   ${IONO}/nemo_fortran_d_n3.o
#
# snapmask_[sd] program
#
${IONB}/snapmask_s   : ${IONO}/snapmask_s.o
	$(CC)  -o ${IONB}/snapmask_s ${IONO}/snapmask_s.o -L$(NEMOLIB) -lnemo -lm

${IONO}/snapmask_s.o : ${IONT}/snapmask2.c
	$(CC) $(CFLAGS) $(OPT) $(INC) -DSINGLEPREC -c ${IONT}/snapmask2.c          \
	                                    -o ${IONO}/snapmask_s.o
${IONB}/snapmask_d   : ${IONO}/snapmask_d.o
	$(CC)  -o  ${IONB}/snapmask_d ${IONO}/snapmask_d.o -L$(NEMOLIB) -lnemo -lm

${IONO}/snapmask_d.o : ${IONT}/snapmask2.c
	$(CC) $(CFLAGS) $(OPT) $(INC) -c ${IONT}/snapmask2.c                       \
                               -o ${IONO}/snapmask_d.o
#
# Dummy function
#
${IONO}/dummy.o : ${IONS}/dummy.c
	$(CC) $(CFLAGS) $(OPT) $(INC) -c ${IONS}/dummy.c -o ${IONO}/dummy.o
#--------------------------------------------------------
#
