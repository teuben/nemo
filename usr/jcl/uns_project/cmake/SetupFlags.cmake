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
if ( NOT SETUP_FLAGS_INCLUDED )
  set( SETUP_FLAGS_INCLUDED 1 )

macro( typed_cache_set type doc var )
    set ( ${var} ${ARGN} CACHE ${type} ${doc} FORCE )
    set ( ${var} ${ARGN} CACHE ${type} ${doc} FORCE )
endmacro()

typed_cache_set ( STRING "setup flags"  SETUP_FLAGS_INCLUDED 1  )

# Set a default build type if none is given
if ( NOT CMAKE_BUILD_TYPE ) # Debug default
    typed_cache_set ( STRING "Build type: Release or Debug" CMAKE_BUILD_TYPE Debug   )
 # else ()                   # else Release
 #   typed_cache_set ( STRING "Build type: Release or Debug" CMAKE_BUILD_TYPE Release )
endif()

# set optimizer flag
if ( ${CMAKE_BUILD_TYPE} STREQUAL "Debug" )
  typed_cache_set (STRING "Optimizer" OPT "-g -ggdb")
else ()
  typed_cache_set (STRING "Optimizer" OPT "-O3")
endif ()

#typed_cache_set (STRING "compilation warnings" WARNCPP "-Waddress -Warray-bounds  -Wc++0x-compat -Wchar-subscripts -Wenum-compare -Wcomment -Wformat -Wmain -Wmissing-braces -Wparentheses -Wreorder -Wreturn-type -Wsequence-point -Wsign-compare -Wstrict-aliasing -Wstrict-overflow=1 -Wswitch -Wtrigraphs -Wuninitialized -Wunknown-pragmas -Wunused-function -Wunused-label -Wunused-value -Wunused-variable -Wvolatile-register-var")
typed_cache_set (STRING "compilation warnings" WARNCPP "-Wall")

# check C compiler flags
include(CheckCCompilerFlag)
set (XX  -Warray-bounds  -Wchar-subscripts -Wcomment -Wformat -Wmain -Wmissing-braces -Wsequence-point -Wstrict-aliasing -Wstrict-overflow=1 -Wswitch -Wtrigraphs -Wunknown-pragmas -Wunused-label -Wvolatile-register-var)
set (WARNC_TMP "")
FOREACH(CWC ${XX})
   MESSAGE(STATUS "Test C options==> " ${CWC})
   check_c_compiler_flag(${CWC} HAS_OPTION)
   if (HAS_OPTION)
     set(WARNC_TMP "${WARNC_TMP} ${CWC}")
   endif()   
endforeach()
typed_cache_set (STRING "compilation warnings" WARNC ${WARNC_TMP}) 
message(STATUS "Warnc = " ${WARNC} )

typed_cache_set (STRING "compilation warnings" WARNF "-Wall" )
# just for testing
option(RELEASE "RELEASE " ON)  # SHARED lib default
if ( NOT RELEASE ) # Debug TRUE
    typed_cache_set ( BOOL "RELEASE Version" RELEASE FALSE )
  else ()                   # else  FALSE
    typed_cache_set ( BOOL "RELEASE Version" RELEASE TRUE )
endif()

# detect MacOS
if    (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
   typed_cache_set ( BOOL "OSX detection" OSX TRUE  )
else()
   typed_cache_set ( BOOL "OSX detection" OSX FALSE )
endif()

# SHARED lib ?
option(BUILD_SHARED_LIBS "Shared lib " ON)  # SHARED lib default
if ( NOT BUILD_SHARED_LIBS ) # Debug TRUE
    typed_cache_set ( BOOL "Build SHARED" BUILD_SHARED_LIBS FALSE )
  else ()                   # else  FALSE
    typed_cache_set ( BOOL "Build SHARED" BUILD_SHARED_LIBS TRUE)
endif()



# set library type , SHARED or STATIC + LIB extension
if ( BUILD_SHARED_LIBS )
   typed_cache_set ( STRING "LIB TYPE" LIBTYPE SHARED )
   if (OSX)
      typed_cache_set ( STRING "LIB TYPE" LIBEXT  "dylib" )
   else()
      typed_cache_set ( STRING "LIB TYPE" LIBEXT  "so" )
   endif()

else ()
   typed_cache_set ( STRING "LIB TYPE" LIBTYPE STATIC )
   typed_cache_set ( STRING "LIB TYPE" LIBEXT  "a" )
endif()

if (UNSIO_INSTALLPATH)
   typed_cache_set ( STRING "UNSIO loation" UNSIOPATH  ${UNSIO_INSTALLPATH} )
endif()
#

endif( NOT SETUP_FLAGS_INCLUDED )
# ============================================================================
