// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// deft.h                                                                      |
//                                                                             |
// C++ code                                                                    |
//                                                                             |
// Copyright Walter Dehnen, 2000-2002                                          |
// e-mail:   wdehnen@aip.de                                                    |
// address:  Astrophysikalisches Institut Potsdam,                             |
//           An der Sternwarte 16, D-14482 Potsdam, Germany                    |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
// defines some constants as default parameter in namespace nbdy::Default      |
//                                                                             |
#ifndef included_deft_h
#define included_deft_h
#ifndef included_enum_h                     //                                 |
#  include <public/enum.h>                  // soft_type, kern_type, MAC_type  |
#endif                                      //                                 |
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// define some default values for nbdy code                                     
//                                                                              
////////////////////////////////////////////////////////////////////////////////
namespace nbdy {
  namespace Default {
#ifdef USE_SSE
#  define DEFAULT_THETA  0.64
#else
#  define DEFAULT_THETA  0.6
#endif
    const real theta       = DEFAULT_THETA;
#define DEFAULT_KERNEL 1
    const kern_type kernel = p1;
#define DEFAULT_MAC    1
    const MAC_type  mac    = theta_of_M;
#ifdef ALLOW_INDI
    const soft_type soften = global;
#endif
#ifdef USE_SSE
//-----------------------------------------------------------------------------+
//                                                                             |
// direct[0] = N_cb^pre:  a C-B interaction is done via direct summation       |
//                        if N_cell <= N_cb^pre                                |
//                                                                             |
// direct[1] = N_cb^post: a not well-separated C-B interaction is done via     |
//                        direct summation if N_cell <= N_cb^post              |
//                                                                             |
// direct[2] = N_cc^post: a not well-separated C-C interaction is done via     |
//                        direct summation if N_cell_1,N_cell_2 <= N_cc^post   |
//                                                                             |
// direct[3] = N_cs:      a cell self-interation is done via direct summation  |
//                        if N_cell <= N_cs                                    |
//                                                                             |
//-----------------------------------------------------------------------------+
# define DEFAULT_NCRIT  16
    const int  Ncrit     = DEFAULT_NCRIT;
# if NDIM==2
    const int  direct[4] = {4, 64, 8,32};
# else
    const int  direct[4] = {4,128,16,64};             
# endif
#else
# define DEFAULT_NCRIT  6
    const int  Ncrit     = DEFAULT_NCRIT;
# if NDIM==2
    const int  direct[4] = {4, 64, 8,32};
# else
    const int  direct[4] = {3,128, 6,64};             
# endif
#endif
  }
}
////////////////////////////////////////////////////////////////////////////////
#endif // included_deft_h
