//-----------------------------------------------------------------------------+
//                                                                             |
// grav.cc                                                                     |
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
// implementing nbdy/grav.h                                                    |
//                                                                             |
//-----------------------------------------------------------------------------+
#include <public/grav.h>

using namespace nbdy;
using namespace nbdy::grav;
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// class nbdy::grav_iact_base                                                 //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
inline grav_iact_base::TaylorSeries::TaylorSeries(vect const&x)
  : X(x)
{
  for(register int i=0; i!=N_COEFF; ++i)
    C[i] = zero;
}

////////////////////////////////////////////////////////////////////////////////
inline grav_iact_base::TaylorSeries::TaylorSeries(TaylorSeries const&T)
  : X(T.X)
{
  for(register int i=0; i!=N_COEFF; ++i)
    C[i] = T.C[i];
}
////////////////////////////////////////////////////////////////////////////////
inline void grav_iact_base::shift_by(real*const&T, vect const&X) const {
  //----------------------------------------------------------------------------
  // computational effort:                                                      
  //                  p   |  1    2    3    4    5    6    7                    
  //                 -----+----------------------------------                   
  //                  *   |  3   18   51  114  222  393  648                    
  //                  +/- |  1   12   35   78  151  266  437                    
  //----------------------------------------------------------------------------
  if(X!=zero && c0(T) != zero ) {                  // 1st order                 
    real         &C0(*T);                          //   C0                      
    register ten1 C1(T+N_C1);                      //   C1                      
    C0 += X * C1;                                  //   1st order shift: C0     
#if falcON_ORDER > 1                               // 2nd order                 
    register vect X2 = half*X;                     //   X/2                     
    register vect F1;                              //   auxiliary Tensor_1      
    register ten2 C2(T+N_C2);                      //   C2                      
    C2.inner_prod(X,F1);                           //   F1 = X * C2             
    C1 += F1;                                      //   2nd order shift: C1     
    C0 += X2 * F1;                                 //   2nd order shift: C0     
#if falcON_ORDER > 2                               // 3rd order                 
    register vect X3 = third*X;                    //   X/3                     
    falcON_SYM2(F2);                               //   auxiliary Tensor_2      
    register ten3 C3(T+N_C3);                      //   C3                      
    C3.inner_prod(X,F2);                           //   F2 = X * C3             
    C2 += F2;                                      //   3rd order shift: C2     
    F2.inner_prod(X2,F1);                          //   F1 = X^2*C3 / 2         
    C1 += F1;                                      //   3rd order shift: C1     
    C0 += X3 * F1;                                 //   3rd order shift: C0     
#if falcON_ORDER > 3                               // 4th order                 
    register vect X4 = quarter*X;                  //   X/4                     
    falcON_SYM3(F3);                               //   auxiliary Tensor_3      
    register ten4 C4(T+N_C4);                      //   C4                      
    C4.inner_prod(X,F3);                           //   F3 = X * C4             
    C3 += F3;                                      //   3rd order shift: C3     
    F3.inner_prod(X2,F2);                          //   F2 = X^2*C4 / 2         
    C2 += F2;                                      //   3rd order shift: C2     
    F2.inner_prod(X3,F1);                          //   F1 = X^3*C4 / 6         
    C1 += F1;                                      //   3rd order shift: C1     
    C0 += X4 * F1;                                 //   3rd order shift: C0     
#if falcON_ORDER > 4                               // 5th order                 
    register vect X5 = fifth*X;                    //   X/5                     
    falcON_SYM4(F4);                               //   auxiliary Tensor_4      
    register ten5 C5(T+N_C5);                      //   C5                      
    C5.inner_prod(X,F4);                           //   F4 = X * C5             
    C4 += F4;                                      //   3rd order shift: C4     
    F4.inner_prod(X2,F3);                          //   F3 = X^2*C5 / 2         
    C3 += F3;                                      //   3rd order shift: C3     
    F3.inner_prod(X3,F2);                          //   F2 = X^3*C5 / 6         
    C2 += F2;                                      //   3rd order shift: C2     
    F2.inner_prod(X4,F1);                          //   F1 = X^4*C5 / 24        
    C1 += F1;                                      //   3rd order shift: C1     
    C0 += X5 * F1;                                 //   3rd order shift: C0     
#endif
#endif
#endif
#endif
  }
}
//------------------------------------------------------------------------------
inline
void grav_iact_base::extract_grav(TaylorSeries const&T, soul_iter const&S) const
{
  //----------------------------------------------------------------------------
  // computational effort:                                                      
  //                  p   |  1    2    3    4    5    6    7                    
  //                 -----+----------------------------------                   
  //                  *   |  3   15   43  114  222  393  648                    
  //                  +/- |  5   15   37   78  151  266  437                    
  //----------------------------------------------------------------------------
  register vect X = cofm(S)-T.X;                   // expansion argument        
  register vect A = c1(T);                         //   1st order acceleration  
  register real P = c0(T) + X * A;                 //   1st order potential     
#if falcON_ORDER > 1                               // 2nd order                 
    register vect F;                               //   auxiliary vector        
    c2(T).inner_prod(X,F);                         //   F = X * C2              
    A += F;                                        //   2nd order acceleration  
    P += half * (X*F);                             //   2nd order potential     
#if falcON_ORDER > 2                               // 3rd order                 
    const real i3= third*half;                     //   1/3!                    
    falcON_SYM2(X2); X2.outer_prod(X);             //   X2ij = Xi*Xj            
    c3(T).inner_prod(X2,F);                        //   F = X^2*C3              
    A += half * F;                                 //   3rd order acceleration  
    P += i3 * (X*F);                               //   3rd order potential     
#if falcON_ORDER > 3                               // 4th order                 
    const real i4= quarter*i3;                     //   1/4!                    
    falcON_SYM3(X3); X3.outer_prod(X2,X);          //   X3ijk = Xi*Xj*Xk        
    c4(T).inner_prod(X3,F);                        //   F = X^3*C4              
    A += i3 * F;                                   //   4th order acceleration  
    P += i4 * (X*F);                               //   4th order potential     
#if falcON_ORDER > 4                               // 5th order                 
    const real i5= fifth*i4;                       //   1/5!                    
    falcON_SYM4(X4); X4.outer_prod(X3,X);          //   X4ijkl = Xi*Xj*Xk*Xl    
    c5(T).inner_prod(X4,F);                        //   F = X^4*C5              
    A += i4 * F;                                   //   5th order acceleration  
    P += i5 * (X*F);                               //   5th order potential     
#endif
#endif
#endif
#endif
  S->acc() += A;                                   // add to souls acceleration 
  S->pot() -= P;                                   // add to souls potential    
}
//------------------------------------------------------------------------------
inline void grav_iact_base::shift_to(TaylorSeries&T, vect const&x) const {
  shift_by(T.C,x-T.X);
  T.X = x;
}
//------------------------------------------------------------------------------
inline void grav_iact_base::shift_and_add(TaylorSeries       &T,
					  vect          const&x,
					  const real*   const&t,
					  real          const&f) const {
  if(t && t[0]) {
    shift_to(T,x);
    for(register int i=0; i!=N_COEFF; ++i) T.C[i] += t[i] * f;
  }
}
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// class nbdy::grav_iact                                                      //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
inline void grav_iact::eval_grav(cell_iter    const&C,
				 TaylorSeries const&T) const
{
  TaylorSeries G(T);                               // G = copy of T             
  shift_and_add(G,cofm(C),coeffs(C),imass(C));     // shift G; G+=T_C           
  LoopSoulKids(cell_iter,C,s) if(is_active(s)) {   // LOOP C's active soul kids 
    s->normalize_grav();                           //   pot,acc/=mass           
    extract_grav(G,s);                             //   add pot,acc due to G    
  }                                                // END LOOP                  
  LoopCellKids(cell_iter,C,c) if(is_active(c))     // LOOP C's active cell kids 
    eval_grav(c,G);                                //   recursive call          
}
//------------------------------------------------------------------------------
void grav_iact::evaluate(cell_iter const&C) const
{
  flush_buffers();                                 // finish interactions       
  eval_grav(C,TaylorSeries(cofm(C)));              // start recursion           
}
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// class nbdy::grav_iact_s                                                    //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
inline void grav_iact_s::eval_grav(cell_iter    const&C,
				   TaylorSeries const&T) const
{
  TaylorSeries G(T);                               // G = copy of T             
  shift_and_add(G,cofm(C),coeffs(C),imass(C));     // shift G; G+=T_C           
  take_coeffs(C);                                  // free memory: C's coeffs   
  LoopSoulKids(cell_iter,C,s) if(is_active(s)) {   // LOOP C's active soul kids 
    s->normalize_grav();                           //   pot,acc/=mass           
    extract_grav(G,s);                             //   add pot,acc due to G    
  }                                                // END LOOP                  
  LoopCellKids(cell_iter,C,c) if(is_active(c))     // LOOP C's active cell kids 
    eval_grav(c,G);                                //   recursive call          
}
//------------------------------------------------------------------------------
void grav_iact_s::evaluate(cell_iter const&Ci) const
{
  flush_buffers();                                 // finish interactions       
  if(NC > MAXNC) MAXNC = NC;                       // record maximum N_coeffs   
  eval_grav(Ci,TaylorSeries(cofm(Ci)));            // start recursion           
}
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// class nbdy::grav_iact_all                                                  //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
inline void grav_iact_all::eval_grav(cell_iter    const&C,
				     TaylorSeries const&T) const
{
  TaylorSeries G(T);                               // G = copy of T             
  shift_and_add(G,cofm(C),coeffs(C),imass(C));     // shift G; G+=T_C           
  LoopSoulKids(cell_iter,C,s) {                    // LOOP C's soul kids        
    s->normalize_grav();                           //   pot,acc/=mass           
    extract_grav(G,s);                             //   add pot,acc due to G    
  }                                                // END LOOP                  
  LoopCellKids(cell_iter,C,c)                      // LOOP C's cell kids        
    eval_grav(c,G);                                //   recursive call          
}
//------------------------------------------------------------------------------
void grav_iact_all::evaluate(cell_iter const&C) const
{
  flush_buffers();                                 // finish interactions       
  eval_grav(C,TaylorSeries(cofm(C)));              // start recursion           
}
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// class nbdy::grav_iact_all_s                                                //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
inline void grav_iact_all_s::eval_grav(cell_iter    const&C,
				   TaylorSeries const&T) const
{
  TaylorSeries G(T);                               // G = copy of T             
  shift_and_add(G,cofm(C),coeffs(C),imass(C));     // shift G; G+=T_C           
  take_coeffs(C);                                  // free memory: C's coeffs   
  LoopSoulKids(cell_iter,C,s) {                    // LOOP C's soul kids        
    s->normalize_grav();                           //   pot,acc/=mass           
    extract_grav(G,s);                             //   add pot,acc due to G    
  }                                                // END LOOP                  
  LoopCellKids(cell_iter,C,c)                      // LOOP C's cell kids        
    eval_grav(c,G);                                //   recursive call          
}
//------------------------------------------------------------------------------
void grav_iact_all_s::evaluate(cell_iter const&Ci) const
{
  flush_buffers();                                 // finish interactions       
  if(NC > MAXNC) MAXNC = NC;                       // record maximum N_coeffs   
  eval_grav(Ci,TaylorSeries(cofm(Ci)));            // start recursion           
}
////////////////////////////////////////////////////////////////////////////////
