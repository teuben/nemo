// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// print_debug.h                                                               |
//                                                                             |
// C++ code                                                                    |
//                                                                             |
// Copyright Jean-Charles LAMBERT - 2005                                       |
// e-mail:   Jean-Charles.Lambert@oamp.fr                                      |
// address:  Dynamique des galaxies                                            |
//           Laboratoire d'Astrophysique de Marseille                          |
//           2, place Le Verrier                                               |
//           13248 Marseille Cedex 4, France                                   |
//           CNRS U.M.R 6110                                                   |
//                                                                             |
//-----------------------------------------------------------------------------+
// Activate debug print statement on screen                                     
//------------------------------------------------------------------------------

#define ALL_DEBUG 0  // set to 1 if you want to debug all files

#if (ALL_DEBUG || LOCAL_DEBUG)
#  define PRINT_D if (1)
#else
#  define PRINT_D if (0)
#endif
//------------------------------------------------------------------------------

