// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// kmac.h                                                                      |
//                                                                             |
// C++ code                                                                    |
//                                                                             |
// Copyright Walter Dehnen, 2000-2003                                          |
// e-mail:   wdehnen@aip.de                                                    |
// address:  Astrophysikalisches Institut Potsdam,                             |
//           An der Sternwarte 16, D-14482 Potsdam, Germany                    |
//                                                                             |
//-----------------------------------------------------------------------------+
#ifndef falcON_included_kmac_h
#define falcON_included_kmac_h
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// Cell-Body interaction with approximated gravity                            //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
#define CELL_BODY(_A,_B,_D0,_D1,_D2,_D3,_R)                                    \
  if(is_active(_A)) {                              /* IF(interaction B->A) >*/ \
    c0(_A)+= _D0;                                  /*   add to C0 of A      */ \
    c1(_A).sub_mul(_R, _D1);                       /*   add to C1 of A      */ \
    c2(_A).add_quadrup(_R,-_D1,_D2);               /*   add to C2 of A      */ \
    c3(_A).add_quadrup(_R,_D2,-_D3);               /*   add to C3 of A      */ \
  }                                                /* <                     */ \
  if(is_active(_B)) {                              /* IF(interaction B->A) >*/ \
    register vect RQ;                              /*   auxiliary vector    */ \
    quad(_A).inner_prod(_R,RQ);                    /*   R*Q_A               */ \
    register real tQ=trace(quad(_A)),RQR=_R*RQ;    /*   tr(Q_A), R*Q*R      */ \
    _B->pot()-=          _D0-_D1*tQ+_D2*RQR;       /*   pot A->B            */ \
    _B->acc().add_mul(_R,_D1-_D2*tQ+_D3*RQR);      /*   acc A->B (in R)     */ \
    RQ       *= _D2;                               /*   R*Q_A*D2            */ \
    _B->acc().sub_twice(RQ);                       /*   acc A->B (in R*Q_A) */ \
  }                                                /* <                     */
//------------------------------------------------------------------------------
#define CELL_BODY_ALL(_A,_B,_D0,_D1,_D2,_D3,_R)                                \
  c0(_A)+= _D0;                                    /* add to C0 of A        */ \
  c1(_A).sub_mul(_R, _D1);                         /* add to C1 of A        */ \
  c2(_A).add_quadrup(_R,-_D1,_D2);                 /* add to C2 of A        */ \
  c3(_A).add_quadrup(_R,_D2,-_D3);                 /* add to C3 of A        */ \
  register vect RQ;                                /* auxiliary vector      */ \
  quad(_A).inner_prod(_R,RQ);                      /* R*Q_A                 */ \
  register real tQ=trace(quad(_A)),RQR=_R*RQ;      /* tr(Q_A), R*Q*R        */ \
  _B->pot()-=          _D0-_D1*tQ+_D2*RQR;         /* pot A->B              */ \
  _B->acc().add_mul(_R,_D1-_D2*tQ+_D3*RQR);        /* acc A->B (in R)       */ \
  RQ       *= _D2;                                 /* R*Q_A*D2              */ \
  _B->acc().sub_twice(RQ);                         /* acc A->B (in R*Q_A)   */
//------------------------------------------------------------------------------
#define CELL_BODY_FP(_A,_B,_D0,_D1,_D2,_D3,_R)				       \
  falcON_SYM3(C3); C3.quadrup(_R,_D2,_D3);         /* pre-compute C3        */ \
  if(is_active(_A)) {                              /* IF(interaction B->A) >*/ \
    c0(_A)+= _D0;                                  /*   add to C0 of A      */ \
    c1(_A).add_mul(_R, _D1);                       /*   add to C1 of A      */ \
    c2(_A).add_quadrup(_R,_D1,_D2);                /*   add to C2 of A      */ \
    c3(_A) +=C3;                                   /*   add to C3 of A      */ \
  }                                                /* <                     */ \
  if(is_active(_B)) {                              /* IF(interaction B->A) >*/ \
    register vect RQ;                              /*   auxiliary vector    */ \
    quad(_A).inner_prod(_R,_RQ);                   /*   R*Q_A               */ \
    register real tQ=trace(quad(_A)),RQR=_R*RQ;    /*   tr(Q_A), R*Q*R      */ \
    _B->pot()-= C3.inner_prod(octo(_A))            /*   add to pot of B     */ \
                      +  _D0+_D1*tQ+_D2*RQR;       /*   pot A->B            */ \
    _B->acc().sub_mul(_R,_D1+_D2*tQ+_D3*RQR);      /*   acc A->B (in R)     */ \
    RQ      *= _D2;                                /*   R*Q_A*D2            */ \
    _B->acc().sub_twice(RQ);                       /*   acc A->B (in R*Q_A) */ \
  }                                                /* <                     */
//------------------------------------------------------------------------------
#define CELL_BODY_FP_ALL(_A,_B,_D0,_D1,_D2,_D3,_R)			       \
  falcON_SYM3(C3); C3.quadrup(_R,_D2,_D3);         /* pre-compute C3        */ \
  c0(_A)+= _D0;                                    /* add to C0 of A        */ \
  c1(_A).add_mul(_R, _D1);                         /* add to C1 of A        */ \
  c2(_A).add_quadrup(_R,_D1,_D2);                  /* add to C2 of A        */ \
  c3(_A) +=C3;                                     /* add to C3 of A        */ \
  register vect RQ;                                /* auxiliary vector      */ \
  quad(_A).inner_prod(_R,_RQ);                     /* R*Q_A                 */ \
  register real tQ=trace(quad(_A)),RQR=_R*RQ;      /* tr(Q_A), R*Q*R        */ \
  _B->pot()-= C3.inner_prod(octo(_A))              /* add to pot of B       */ \
                    +  _D0+_D1*tQ+_D2*RQR;         /* pot A->B              */ \
  _B->acc().sub_mul(_R,_D1+_D2*tQ+_D3*RQR);        /* acc A->B (in R)       */ \
  RQ      *= _D2;                                  /* R*Q_A*D2              */ \
  _B->acc().sub_twice(RQ);                         /* acc A->B (in R*Q_A)   */
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// Cell-Cell interaction with approximated gravity                            //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
#define CELL_CELL(_A,_B,_D0,_D1,_D2,_D3,_R)                                    \
  falcON_SYM2(C2); C2.quadrup(_R,-_D1,_D2);        /* pre-compute C2        */ \
  falcON_SYM3(C3); C3.quadrup(_R,_D2,-_D3);        /* pre-compute C3        */ \
  register vect RQ;                                /* to hold R*Q           */ \
  register real tQ, RQR;                           /* to hold tr(Q), R*Q*R  */ \
  if(is_active(_A)) {                              /* IF(interaction B->A) >*/ \
    quad(_B).inner_prod(_R,RQ);                    /*   R*Q_B               */ \
    tQ    = trace(quad(_B));                       /*   tr(Q_B)             */ \
    RQR   = _R*RQ;                                 /*   R*Q_B*R             */ \
    c0(_A)+=           _D0 - _D1*tQ + _D2*RQR;     /*   add to C0 of A      */ \
    c1(_A).sub_mul(_R, _D1 - _D2*tQ + _D3*RQR);    /*   add to C1 of A (R)  */ \
    RQ   *= _D2;                                   /*   R*Q_B*D2            */ \
    c1(_A).add_twice(RQ);                          /*   add to C1 of A(R*QB)*/ \
    c2(_A)+= C2;                                   /*   add to C2 of A      */ \
    c3(_A)+= C3;                                   /*   add to C3 of A      */ \
  }                                                /* <                     */ \
  if(is_active(_B)) {                              /* IF(interaction B->A) >*/ \
    _R.negate();                                   /*   set R to -R         */ \
    quad(_A).inner_prod(_R,RQ);                    /*   R*Q_B               */ \
    tQ    = trace(quad(_A));                       /*   tr(Q_A)             */ \
    RQR   = _R*RQ;                                 /*   R*Q_A*R             */ \
    c0(_B)+=           _D0 - _D1*tQ + _D2*RQR;     /*   add to C0 of B      */ \
    c1(_B).sub_mul(_R, _D1 - _D2*tQ + _D3*RQR);    /*   add to C1 of B (R)  */ \
    RQ   *= _D2;                                   /*   R*Q_A*D2            */ \
    c1(_B).add_twice(RQ);                          /*   add to C1 of B(R*QA)*/ \
    c2(_B)+= C2;                                   /*   add to C2 of B      */ \
    c3(_B)-= C3;                                   /*   add to C3 of B      */ \
  }                                                /* <                     */
//------------------------------------------------------------------------------
#define CELL_CELL_ALL(_A,_B,_D0,_D1,_D2,_D3,_R)                                \
  falcON_SYM2(C2); C2.quadrup(_R,-_D1,_D2);        /* pre-compute C2        */ \
  falcON_SYM3(C3); C3.quadrup(_R,_D2,-_D3);        /* pre-compute C3        */ \
  register vect RQ;                                /* to hold R*Q           */ \
  register real tQ, RQR;                           /* to hold tr(Q), R*Q*R  */ \
  quad(_B).inner_prod(_R,RQ);                      /* R*Q_B                 */ \
  tQ    = trace(quad(_B));                         /* tr(Q_B)               */ \
  RQR   = _R*RQ;                                   /* R*Q_B*R               */ \
  c0(_A)+=           _D0 - _D1*tQ + _D2*RQR;       /* add to C0 of A        */ \
  c1(_A).sub_mul(_R, _D1 - _D2*tQ + _D3*RQR);      /* add to C1 of A (R)    */ \
  RQ   *= _D2;                                     /* R*Q_B*D2              */ \
  c1(_A).add_twice(RQ);                            /* add to C1 of A(R*QB)  */ \
  c2(_A)+= C2;                                     /* add to C2 of A        */ \
  c3(_A)+= C3;                                     /* add to C3 of A        */ \
  _R.negate();                                     /* set R to -R           */ \
  quad(_A).inner_prod(_R,RQ);                      /* R*Q_B                 */ \
  tQ    = trace(quad(_A));                         /* tr(Q_A)               */ \
  RQR   = _R*RQ;                                   /* R*Q_A*R               */ \
  c0(_B)+=           _D0 - _D1*tQ + _D2*RQR;       /* add to C0 of B        */ \
  c1(_B).sub_mul(_R, _D1 - _D2*tQ + _D3*RQR);      /* add to C1 of B (R)    */ \
  RQ   *= _D2;                                     /* R*Q_A*D2              */ \
  c1(_B).add_twice(RQ);                            /* add to C1 of B(R*QA)  */ \
  c2(_B)+= C2;                                     /* add to C2 of B        */ \
  c3(_B)-= C3;                                     /* add to C3 of B        */
//------------------------------------------------------------------------------
#define CELL_CELL_FP(_A,_B,_D0,_D1,_D2,_D3,_R)				       \
  falcON_SYM2(C2); C2.quadrup(_R,_D1,_D2);         /* pre-compute C2        */ \
  falcON_SYM3(C3); C3.quadrup(_R,_D2,_D3);         /* pre-compute C3        */ \
  register vect RQ;                                /* to hold R*Q           */ \
  register real tQ, RQR;                           /* to hold tr(Q), R*Q*R  */ \
  if(is_active(_A)) {                              /* IF(interaction B->A) >*/ \
    quad(_B).inner_prod(_R,RQ);                    /*   R*Q_B               */ \
    tQ    = trace(quad(_B));                       /*   tr(Q_B)             */ \
    RQR   = _R*RQ;                                 /*   R*Q_B*R             */ \
    c0(_A)+=           _D0 + _D1*tQ + _D2*RQR      /*   add to C0 of A      */ \
                     - C3.inner_prod(octo(_B));    /*   add to C0 of A      */ \
    c1(_A).add_mul(_R, _D1 + _D2*tQ + _D3*RQR);    /*   add to C1 of A (R)  */ \
    RQ   *= _D2;                                   /*   R*Q_B*D2            */ \
    c1(_A).add_twice(RQ);                          /*   add to C1 of A(R*QB)*/ \
    c2(_A)+=_ C2;                                  /*   add to C2 of A      */ \
    c3(_A)+= C3;                                   /*   add to C3 of A      */ \
  }                                                /* <                     */ \
  if(is_active(_B)) {                              /* IF(interaction B->A) >*/ \
    _R.negate();                                   /*   set R to -R         */ \
    quad(_A).inner_prod(_R,RQ);                    /*   R*Q_B               */ \
    tQ    = trace(quad(_A));                       /*   tr(Q_A)             */ \
    RQR   = _R*RQ;                                 /*   R*Q_A*R             */ \
    c0(_B)+=           _D0 + _D1*tQ + _D2*RQR;     /*   add to C0 of B      */ \
                     + C3.inner_prod(octo(_A));    /*   add to C0 of B      */ \
    c1(_B).add_mul(_R, _D1 + _D2*tQ + _D3*RQR);    /*   add to C1 of B (R)  */ \
    RQ   *= _D2;                                   /*   R*Q_A*D2            */ \
    c1(_B).add_twice(RQ);                          /*   add to C1 of B(R*QA)*/ \
    c2(_B)+= C2;                                   /*   add to C2 of B      */ \
    c3(_B)-= C3;                                   /*   add to C3 of B      */ \
  }                                                /* <                     */
//------------------------------------------------------------------------------
#define CELL_CELL_FP_ALL(_A,_B,_D0,_D1,_D2,_D3,_R)			       \
  falcON_SYM2(C2); C2.quadrup(_R,_D1,_D2);         /* pre-compute C2        */ \
  falcON_SYM3(C3); C3.quadrup(_R,_D2,_D3);         /* pre-compute C3        */ \
  register vect RQ;                                /* to hold R*Q           */ \
  register real tQ, RQR;                           /* to hold tr(Q), R*Q*R  */ \
  quad(_B).inner_prod(_R,RQ);                      /* R*Q_B                 */ \
  tQ    = trace(quad(_B));                         /* tr(Q_B)               */ \
  RQR   = _R*RQ;                                   /* R*Q_B*R               */ \
  c0(_A)+=           _D0 + _D1*tQ + _D2*RQR        /* add to C0 of A        */ \
                   - C3.inner_prod(octo(_B));      /* add to C0 of A        */ \
  c1(_A).add_mul(_R, _D1 + _D2*tQ + _D3*RQR);      /* add to C1 of A (R)    */ \
  RQ   *= _D2;                                     /* R*Q_B*D2              */ \
  c1(_A).add_twice(RQ);                            /* add to C1 of A(R*QB)  */ \
  c2(_A)+=_ C2;                                    /* add to C2 of A        */ \
  c3(_A)+= C3;                                     /* add to C3 of A        */ \
  _R.negate();                                     /* set R to -R           */ \
  quad(_A).inner_prod(_R,RQ);                      /* R*Q_B                 */ \
  tQ    = trace(quad(_A));                         /* tr(Q_A)               */ \
  RQR   = _R*RQ;                                   /* R*Q_A*R               */ \
  c0(_B)+=           _D0 + _D1*tQ + _D2*RQR;       /* add to C0 of B        */ \
                   + C3.inner_prod(octo(_A));      /* add to C0 of B        */ \
  c1(_B).add_mul(_R, _D1 + _D2*tQ + _D3*RQR);      /* add to C1 of B (R)    */ \
  RQ   *= _D2;                                     /* R*Q_A*D2              */ \
  c1(_B).add_twice(RQ);                            /* add to C1 of B(R*QA)  */ \
  c2(_B)+= C2;                                     /* add to C2 of B        */ \
  c3(_B)-= C3;                                     /* add to C3 of B        */
////////////////////////////////////////////////////////////////////////////////
#endif                                             // falcON_included_kmac_h    
