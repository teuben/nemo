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

# require boost
find_package( Boost 1.36.0)
set(LIBBOOST boost_thread-mt.a boost_system-mt.a)
set(Boost_FOUND FALSE)
set(DNOBOOST "")
if(NOT Boost_FOUND)
   MESSAGE (STATUS " Boost not found, uns_2dplot will run slowly.....")
   set(DNOBOOST "-DNOBOOST")
   set(LIBBOOST "")
endif(NOT Boost_FOUND)
