// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// coef.h                                                                      |
//                                                                             |
// C++ code                                                                    |
//                                                                             |
// Copyright Walter Dehnen, 2000-2003                                          |
// e-mail:   walter.dehnen@astro.le.ac.uk                                      |
// address:  Department of Physics and Astronomy, University of Leicester      |
//           University Road, Leicester LE1 7RH, United Kingdom                |
//                                                                             |
//-----------------------------------------------------------------------------+
#ifndef falcON_included_coef_h
#define falcON_included_coef_h

#ifndef falcON_included_grav_h
#  include <public/grav.h>
#endif
#ifndef falcON_included_tset_h
#  include <public/tset.h>
#endif


#ifdef falcON_SSE_CODE
#  define __ARG_D D,J
#else
#  define __ARG_D D
#endif

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// Cell-Leaf interactions                                                     //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#define CellLeaf(A,B,D,J,R) {						       \
  register grav::Cset F;                           /* to hold F^(n)        */  \
  if(is_active(A)) {                               /* IF A is active       */  \
    set_dPhi(F,R,__ARG_D);                         /*   F^(n) = d^nPhi/dR^n*/  \
    add_C_B2C(A->Coeffs(),F);                      /*   C_A   = ...        */  \
    if(is_active(B)) {                             /*   IF B is active, too*/  \
      F.flip_sign_odd();                           /*     flip sign:F^(odd)*/  \
      add_C_C2B(B->Coeffs(),F,A->poles());         /*     C_B   = ...      */  \
    }                                              /*   ENDIF              */  \
  } else if(is_active(B)) {                        /* ELIF B is active     */  \
    R.negate();                                    /*   flip sign: R       */  \
    set_dPhi(F,R,__ARG_D);                         /*   F^(n) = d^nPhi/dR^n*/  \
    add_C_C2B(B->Coeffs(),F,A->poles());           /*   C_B   = ...        */  \
  }                                                /* ENDIF                */  \
}

#define CellLeafAll(A,B,D,J,R) {					       \
  register grav::Cset F;                           /* to hold F^(n)        */  \
  set_dPhi(F,R,__ARG_D);                           /* F^(n) = d^nPhi/dR^n  */  \
  add_C_B2C(A->Coeffs(),F);                        /* C_A   = ...          */  \
  F.flip_sign_odd();                               /* F^(n) = d^nPhi/dR^n  */  \
  add_C_C2B(B->Coeffs(),F,A->poles());             /* C_B   = ...          */  \
}

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// Cell-Cell interactions                                                     //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#define CellCell(A,B,D,J,R) {						       \
  register grav::Cset F;                           /* to hold F^(n)        */  \
  if(is_active(A)) {                               /* IF A is active       */  \
    set_dPhi(F,R,__ARG_D);                         /*   F^(n) = d^nPhi/dR^n*/  \
    add_C_C2C(A->Coeffs(),F,B->poles());           /*   C_A   = ...        */  \
    if(is_active(B)) {                             /*   IF B is active, too*/  \
      F.flip_sign_odd();                           /*     flip sign:F^(odd)*/  \
      add_C_C2C(B->Coeffs(),F,A->poles());         /*     C_B   = ...      */  \
    }                                              /*   ENDIF              */  \
  } else if(is_active(B)) {                        /* ELIF B is active     */  \
    R.negate();                                    /*   flip sign: R       */  \
    set_dPhi(F,R,__ARG_D);                         /*   F^(n) = d^nPhi/dR^n*/  \
    add_C_C2C(B->Coeffs(),F,A->poles());           /*   C_B   = ...        */  \
  }                                                /* ENDIF                */  \
}

#define CellCellAll(A,B,D,J,R) {					       \
  register grav::Cset F;                           /* to hold F^(n)        */  \
  set_dPhi(F,R,__ARG_D);                           /* F^(n) = d^nPhi/dR^n  */  \
  add_C_C2C(A->Coeffs(),F,B->poles());             /* C_A   = ...          */  \
  F.flip_sign_odd();                               /* F^(n) = d^nPhi/dR^n  */  \
  add_C_C2C(B->Coeffs(),F,A->poles());             /* C_B   = ...          */  \
}

////////////////////////////////////////////////////////////////////////////////
#endif                                             // falcON_included_coef_h    
