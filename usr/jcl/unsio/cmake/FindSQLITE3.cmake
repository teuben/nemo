# ============================================================================
# Copyright Jean-Charles LAMBERT - 2009-2012
# e-mail:   Jean-Charles.Lambert@oamp.fr
# address:  Dynamique des galaxies
#           Centre de donneeS Astrophysique de Marseille (CeSAM)
#           Laboratoire d'Astrophysique de Marseille
#           Pole de l'Etoile, site de Chateau-Gombert
#           38, rue Frederic Joliot-Curie
#           13388 Marseille cedex 13 France
#           CNRS U.M.R 6110
# ============================================================================
# CMake module to detect SQLITE3 library
# ============================================================================


SET(SQLITE3_FOUND FALSE)
SET(DNOSQLITE3 "-DNOSQLITE3") # if sqlite3 does not exist
SET(SQLITE3_LIB_PATH "")
SET(SQLITE3_H_PATH "")
SET(SQLITE3_LIB "")

find_file(SQLITE3_H NAMES "sqlite3.h" PATH /usr/local/include)
IF (SQLITE3_H)
  get_filename_component(SQLITE3_H_PATH  ${SQLITE3_H} PATH)
  MESSAGE(STATUS "Found sqlite3.h:" ${SQLITE3_H})
  MESSAGE(STATUS "Found sqlite3.h path:" ${SQLITE3_H_PATH})
  find_library(SQLITE3 NAMES sqlite3 PATH /usr/local/lib)
  IF (SQLITE3)
    SET(SQLITE3_FOUND TRUE)
    SET(DNOSQLITE3 "") # SQLITE3 exist
    SET(SQLITE3_LIB sqlite3)
    MESSAGE(STATUS "Found library here :" ${SQLITE3})
    get_filename_component(SQLITE3_LIB_PATH  ${SQLITE3} PATH)
    MESSAGE(STATUS "Found library PATH :" ${SQLITE3_LIB_PATH})    
  ENDIF(SQLITE3)
ENDIF(SQLITE3_H) 
