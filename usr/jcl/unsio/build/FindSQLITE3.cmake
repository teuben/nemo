# --------------------------------------------------------
# CMake module to detect SQLITE3 library                   
# --------------------------------------------------------

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
  ENDIF()
ENDIF() 
