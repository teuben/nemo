// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// deft.h                                                                      |
//                                                                             |
// C++ code                                                                    |
//                                                                             |
// Copyright Walter Dehnen, 2000-2003                                          |
// e-mail:   wdehnen@aip.de                                                    |
// address:  Astrophysikalisches Institut Potsdam,                             |
//           An der Sternwarte 16, D-14482 Potsdam, Germany                    |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
// defines some constants as default parameter in namespace nbdy::Default      |
//                                                                             |
//-----------------------------------------------------------------------------+
#ifndef falcON_included_deft_h
#define falcON_included_deft_h
#ifndef falcON_included_enum_h              //                                 |
#  include <public/enum.h>                  // soft_type, kern_type, MAC_type  |
#endif                                      //                                 |
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// define some default values for nbdy code                                     
//                                                                              
////////////////////////////////////////////////////////////////////////////////
namespace nbdy {
  namespace Default {
    //--------------------------------------------------------------------------
    // the default for the tolerance parameter theta  (depends on SSE)          
    //--------------------------------------------------------------------------
#ifdef falcON_SSE
#  define falcON_THETA_TEXT "0.64"
    const real theta       = 0.64;
#else
#  define falcON_THETA_TEXT "0.60"
    const real theta       = 0.6;
#endif
    //--------------------------------------------------------------------------
    // the default kernel                                                       
    //--------------------------------------------------------------------------
#define falcON_KERNEL_TEXT  "1"
    const kern_type kernel = p1;
    //--------------------------------------------------------------------------
    // the default multipole acceptance criterion (MAC)                         
    //--------------------------------------------------------------------------
#define falcON_MAC_TEXT     "1"
    const MAC_type  mac    = theta_of_M;
    //--------------------------------------------------------------------------
    // the default behaviour of softening (global/individual)                   
    //--------------------------------------------------------------------------
#ifdef falcON_INDI
    const soft_type soften = global;
#endif
    //--------------------------------------------------------------------------
    // the default N_crit, the max # number of bodies in a unsplit cell         
    //--------------------------------------------------------------------------
#ifdef falcON_SSE
#  define falcON_NCRIT_TEXT  "16"
    const int  Ncrit        = 16;
#else
#  define falcON_NCRIT_TEXT  "6"
    const int  Ncrit        = 6;
#endif
#ifdef falcON_SPH
#  define falcON_SPHNCRIT_TEXT "32"
    const int  SPHNcrit       = 32;
#endif
    //--------------------------------------------------------------------------
    // the default maximum tree depth                                           
    //--------------------------------------------------------------------------
#define falcON_MAXDEPTH_TEXT "100"
    const int  MaxDepth     = 100;
    //--------------------------------------------------------------------------
    // the default values for the parameter controlling direct summation:       
    //-------------------------------------------------------------------------+
    //                                                                         |
    // direct[0] = N_cb^pre:  a C-B interaction is done via direct summation   |
    //                        if N_cell <= N_cb^pre                            |
    //                                                                         |
    // direct[1] = N_cb^post: a not well-separated C-B interaction is done via |
    //                        direct summation if N_cell <= N_cb^post          |
    //                                                                         |
    // direct[2] = N_cc^post: a not well-separated C-C interaction is done via |
    //                        direct summation if N_cell_1,N_cell_2<=N_cc^post |
    //                                                                         |
    // direct[3] = N_cs:      a cell self-interation is done via direct sum    |
    //                        if N_cell <= N_cs                                |
    //                                                                         |
    //-------------------------------------------------------------------------+
#ifdef falcON_SSE
#  if falcON_NDIM==2
    const int  direct[4] = {4, 64, 8,32};
#  else
    const int  direct[4] = {4,128,16,64};             
#  endif
#else
#  if falcON_NDIM==2
    const int  direct[4] = {4, 64, 8,32};
#  else
    const int  direct[4] = {3,128, 6,64};             
#  endif
#endif
#ifdef falcON_SPH
    const int  SPHdirect[3] = {128,32,64};
#endif
  }
}
////////////////////////////////////////////////////////////////////////////////
#endif // falcON_included_deft_h
