// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// enum.h                                                                      |
//                                                                             |
// C++ code                                                                    |
//                                                                             |
// Copyright Walter Dehnen, 2002-2003                                          |
// e-mail:   walter.dehnen@astro.le.ac.uk                                      |
// address:  Department of Physics and Astronomy, University of Leicester      |
//           University Road, Leicester LE1 7RH, United Kingdom                |
//                                                                             |
//-----------------------------------------------------------------------------+
#ifndef falcON_included_enum_h
#define falcON_included_enum_h 1
namespace nbdy {
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // enum nbdy::kern_type                                                     //
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
  // enum nbdy::MAC_type                                                      //
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
}                                                  // END: namespace nbdy       
#endif                                             // falcON_included_enum_h    
