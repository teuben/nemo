#-*-cmake-*-
# ============================================================================
# Copyright Jean-Charles LAMBERT - 2008-2013
# e-mail:   Jean-Charles.Lambert@oamp.fr
# address:  Dynamique des galaxies
#           Centre de donneeS Astrophysique de Marseille (CeSAM)
#           Laboratoire d'Astrophysique de Marseille
#           Pole de l'Etoile, site de Chateau-Gombert
#           38, rue Frederic Joliot-Curie
#           13388 Marseille cedex 13 France
#           CNRS U.M.R 6110
# ============================================================================
# CMakeListst.txt file for nemo library
# ============================================================================
SET(NEMO_INSTALLED FALSE)

FILE(GLOB GLOB_TEMP_VAR $ENV{NEMO})
IF(GLOB_TEMP_VAR)
  SET(NEMO_INSTALLED TRUE)

  # Detect fortran compiler from ${NEMOLIB}/makedefs
  # files (variable FC  = )

  FILE(READ $ENV{NEMOLIB}/makedefs MAKEDEFS)
  STRING(REGEX MATCH "FC  = gfortran" GFORTRAN_IS_SET ${MAKEDEFS}) 


  IF ("FC  = gfortran" STREQUAL "${GFORTRAN_IS_SET}")
    MESSAGE("GFortran compiler Detected......")
    SET (FORT_COMP gfortran)
  ELSE("FC  = gfortran" STREQUAL "${GFORTRAN_IS_SET}")
    # g77 stuffs     
    MESSAGE("Assuming G77 compiler")
    SET (FORT_COMP g77)
  ENDIF ("FC  = gfortran" STREQUAL "${GFORTRAN_IS_SET}")

  # gfortran stuffs
  SET (FC_GFORT_LIB "")
  SET (FC_GFORT_PATH "")
  set (GFORTPATH "")
  execute_process(COMMAND gfortran -print-file-name=libgfortran.a 
    OUTPUT_VARIABLE GFORTPATH OUTPUT_STRIP_TRAILING_WHITESPACE)
  IF   (NOT ${GFORTPATH} STREQUAL "")
    MESSAGE("GFORTRAN compiler exist, adding libgfortran......")
    execute_process(COMMAND dirname ${GFORTPATH} 
      OUTPUT_VARIABLE FC_GFORT_PATH OUTPUT_STRIP_TRAILING_WHITESPACE)
    SET (FC_GFORT_LIB gfortran)
  ENDIF(NOT ${GFORTPATH} STREQUAL "")
  
  # g77 stuffs     
  SET (FC_G77_LIB "")
  SET (FC_G77_PATH "")
  set (G2CFPATH "")
  execute_process(COMMAND g77 -print-file-name=libg2c.a 
    OUTPUT_VARIABLE G2CFPATH OUTPUT_STRIP_TRAILING_WHITESPACE)
  IF   (NOT ${G2CFPATH} STREQUAL "")
    MESSAGE("G77 compiler exist, adding libg2c......")
    execute_process(COMMAND dirname ${G2CFPATH} 
      OUTPUT_VARIABLE FC_G77_PATH OUTPUT_STRIP_TRAILING_WHITESPACE)
    SET (FC_G77_LIB  g2c)
  ENDIF(NOT ${G2CFPATH} STREQUAL "")

ENDIF(GLOB_TEMP_VAR)

if (NOT NEMO_INSTALLED)
  MESSAGE(STATUS "NEMO environement not loaded.... using nemo_light")
endif (NOT NEMO_INSTALLED)
