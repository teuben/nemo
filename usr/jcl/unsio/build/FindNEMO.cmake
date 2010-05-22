SET(NEMO_INSTALLED FALSE)

FILE(GLOB GLOB_TEMP_VAR $ENV{NEMO})
IF(GLOB_TEMP_VAR)
  SET(NEMO_INSTALLED TRUE)

  # Detect fortran compiler from ${NEMOLIB}/makedefs
  # files (variable FC  = )

  FILE(READ $ENV{NEMOLIB}/makedefs MAKEDEFS)
  STRING(REGEX MATCH "FC  = gfortran" GFORTRAN_IS_SET ${MAKEDEFS}) 

  
  IF ("FC  = gfortran" STREQUAL "${GFORTRAN_IS_SET}")
     # gfortran stuffs
     MESSAGE("GFortran compiler Detected......")
     #SET (FC_COMPILER gfortran)
     execute_process(COMMAND gfortran -print-file-name=libgfortran.a 
	OUTPUT_VARIABLE G2CFPATH OUTPUT_STRIP_TRAILING_WHITESPACE)
     execute_process(COMMAND dirname ${G2CFPATH} 
	OUTPUT_VARIABLE FC_LIB_PATH OUTPUT_STRIP_TRAILING_WHITESPACE)
     SET (FC_LIB gfortran)
  ELSE("FC  = gfortran" STREQUAL "${GFORTRAN_IS_SET}")
     # g77 stuffs     
     MESSAGE("Assuming G77 compiler")
     #SET (FC_COMPILER g77)
     execute_process(COMMAND g77 -print-file-name=libg2c.a 
	OUTPUT_VARIABLE G2CFPATH OUTPUT_STRIP_TRAILING_WHITESPACE)
     execute_process(COMMAND dirname ${G2CFPATH} 
	OUTPUT_VARIABLE FC_LIB_PATH OUTPUT_STRIP_TRAILING_WHITESPACE)
     SET (FC_LIB g2c)
  ENDIF  ("FC  = gfortran" STREQUAL "${GFORTRAN_IS_SET}")

  MESSAGE( STATUS "FC LIB PATH : " ${FC_LIB_PATH})

ENDIF(GLOB_TEMP_VAR)

if (NOT NEMO_INSTALLED)
MESSAGE(STATUS "NEMO environement not loaded.... using nemo_light")
endif (NOT NEMO_INSTALLED)
