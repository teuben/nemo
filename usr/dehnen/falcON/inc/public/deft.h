// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// deft.h                                                                      |
//                                                                             |
// C++ code                                                                    |
//                                                                             |
// Copyright Walter Dehnen, 2000-2003                                          |
// e-mail:   walter.dehnen@astro.le.ac.uk                                      |
// address:  Department of Physics and Astronomy, University of Leicester      |
//           University Road, Leicester LE1 7RH, United Kingdom                |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
// defines some constants as default parameter in namespace nbdy::Default      |
//                                                                             |
//-----------------------------------------------------------------------------+
#ifndef falcON_included_deft_h
#define falcON_included_deft_h
#ifndef falcON_included_auxx_h
#  include <public/auxx.h>
#endif
#ifndef falcON_included_enum_h              //                                 |
#  include <public/enum.h>                  // kern_type, MAC_type             |
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
#define falcON_KERNEL_NAME  "P1"
    const kern_type kernel = p1;
  }
  //----------------------------------------------------------------------------
  inline kern_type kern(const indx KERN) {
    switch(KERN % 10) { 
    case  0: return p0;
    case  1: return p1;
    case  2: return p2;
    case  3: return p3;
    case  9: return newton;
    default: warning("kernel unknown, defaulting to " falcON_KERNEL_NAME );
      return Default::kernel;
    }
  }
  //----------------------------------------------------------------------------
  namespace Default {
    //--------------------------------------------------------------------------
    // the default multipole acceptance criterion (MAC)                         
    //--------------------------------------------------------------------------
#define falcON_MAC_TEXT     "1"
    const MAC_type  mac    = theta_of_M;
    //--------------------------------------------------------------------------
    // the default behaviour of softening (global/individual)                   
    //--------------------------------------------------------------------------
    const bool soften      = false;
    //--------------------------------------------------------------------------
    // the default maximum tree depth                                           
    //--------------------------------------------------------------------------
#define falcON_MAXDEPTH_TEXT "100"
    const int  MaxDepth     = 100;
    //--------------------------------------------------------------------------
    // the default values for the parameter controlling direct summation:       
    //-------------------------------------------------------------------------+
    //                                                                         |
    // Ncrit:                 maximum number of bodies in a unsplit cell       |
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
#ifdef falcON_SSE_CODE
# if   falcON_ORDER == 3
#   define falcON_NCRIT_TEXT  "16"
    const int  Ncrit         = 16;
    const int  direct[4]     = {4,128,16,64};             
#  elif falcON_ORDER == 4
#   define falcON_NCRIT_TEXT  "32"
    const int  Ncrit         = 32;
    const int  direct[4]     = {32,256,32,128};             
# else
#   define falcON_NCRIT_TEXT  "48"
    const int  Ncrit         = 48;
    const int  direct[4]     = {48,512,48,256};             
# endif

#else  // ! falcON_SSE_CODE

# if   falcON_ORDER == 3
#  define falcON_NCRIT_TEXT  "6"
    const int  Ncrit        = 6;
    const int  direct[4]    = {3,128, 6,64};             
# else
#  define falcON_NCRIT_TEXT  "20"
    const int  Ncrit        = 20;
    const int  direct[4]    = {20,128,20,64};             
# endif
#endif //   falcON_SSE_CODE

    //-------------------------------------------------------------------------+
    // the default values for the parameters controlling SPH code              |
    //-------------------------------------------------------------------------+
#ifdef falcON_SPH
#  define falcON_SPHNCRIT_TEXT "32"
    const int  SPHNcrit     =   32;
    const int  SPHdirect[3] = {128,32,64};
#endif
  }
}
////////////////////////////////////////////////////////////////////////////////
#endif // falcON_included_deft_h
