// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// default.h                                                                   |
//                                                                             |
// Copyright (C) 2000, 2001, 2002, 2003, 2008  Walter Dehnen                   |
//                                                                             |
// This program is free software; you can redistribute it and/or modify        |
// it under the terms of the GNU General Public License as published by        |
// the Free Software Foundation; either version 2 of the License, or (at       |
// your option) any later version.                                             |
//                                                                             |
// This program is distributed in the hope that it will be useful, but         |
// WITHOUT ANY WARRANTY; without even the implied warranty of                  |
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU           |
// General Public License for more details.                                    |
//                                                                             |
// You should have received a copy of the GNU General Public License           |
// along with this program; if not, write to the Free Software                 |
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.                   |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
// defines some constants as default parameter in namespace falcON::Default    |
//                                                                             |
//-----------------------------------------------------------------------------+
#ifndef falcON_included_default_h
#define falcON_included_default_h

#ifndef falcON_included_basic_h
#  include <public/basic.h>
#endif
////////////////////////////////////////////////////////////////////////////////
namespace falcON {
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // enum falcON::kern_type                                                   //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  enum kern_type {                                 // type of kernel            
    p0     = 0,                                    // P0 = Plummer kernel       
    p1     = 1,                                    // P1,P2,P3: enhanced Plummer
    p2     = 2,                                    //                           
    p3     = 3,                                    //                           
    newton = 9                                     // no softening              
  };                                               //                           
  //----------------------------------------------------------------------------
  inline
  const char* describe(const kern_type KER)        // describe kernel           
  {
    switch(KER) {
    case p0: return "P0==Plummer";
    case p1: return "P1";
    case p2: return "P2";
    case p3: return "P3";
    case newton:
    default: return "Newtonian";
    }
  }
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // enum falcON::soft_type                                                   //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  enum soft_type {                                 // type of softening         
    global_fixed        = 0,                       // globally time-constant    
    individual_fixed    = 1,                       // individual but time-const 
#ifdef falcON_ADAP
    individual_adaptive = 2                        // individual and adaptive   
#endif
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // enum falcON::MAC_type                                                    //
  //                                                                          //
  // NOTE The latter two possibilities are experimental and should not be     //
  // used by the general user (he won't know what theta_0 to use with them).  //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  enum MAC_type {                                  // type of MAC               
    const_theta     =0,                            // theta = const             
    theta_of_M      =1,                            // theta = theta(M)          
    theta_of_M_ov_rq=2, /* experimental */         // theta = theta(M/rmax^2)   
    theta_of_M_ov_r =3  /* experimental */         // theta = theta(M/rmax)     
  };                                               //                           
  inline
  const char* describe(const MAC_type MAC)         // describe MAC              
  {
    switch(MAC) {
    case const_theta:      return "theta=const";
    case theta_of_M:       return "theta(M)";
    case theta_of_M_ov_r:  return "theta(M/rmax)";
    case theta_of_M_ov_rq: return "theta(M/rmax^2)";
    default:               return "unknown MAC";
    }
  }
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // define some default values for falcON code                               //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
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
    default:
      falcON_Warning("kernel unknown, defaulting to " falcON_KERNEL_NAME );
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
  } // namespace Default {
} // namespace falcON {
////////////////////////////////////////////////////////////////////////////////
#endif // falcON_included_default_
