// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// tens.h                                                                      |
//                                                                             |
// C++ code                                                                    |
//                                                                             |
// Copyright Walter Dehnen, 1999-2003                                          |
// e-mail:   wdehnen@aip.de                                                    |
// address:  Astrophysikalisches Institut Potsdam,                             |
//           An der Sternwarte 16, D-14482 Potsdam, Germany                    |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
// defines                                                                     |
// class nbdy::sym1<SCALAR>                                                    |
// class nbdy::sym2<SCALAR>                                                    |
// class nbdy::sym3<SCALAR>                                                    |
// class nbdy::sym4<SCALAR>                                                    |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
// last change 28/07/03: happy gcc 3.3                                         |
//                                                                             |
//-----------------------------------------------------------------------------+
#ifndef falcON_included_tens_h
#define falcON_included_tens_h

#ifndef falcON_included_vect_h
#  include <public/vect.h>
#endif
////////////////////////////////////////////////////////////////////////////////
#ifndef falcON_NDIM
#  error falcON_NDIM not #defined in tens.h
#endif
////////////////////////////////////////////////////////////////////////////////
namespace nbdy {
#define OP      operator
#define TV      template<typename VECTOR>
#define c_VECT  const VECTOR&
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class nbdy::sym1                                                         //
  //                                                                          //
  //--------------------------------------------------------------------------//
  //                                                                          //
  // a vector with external memory                                            //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
#if falcON_NDIM == 2

#define CI(op,b,ox) A[0] op b[0] ox  A[1] op b[1]
#define CX(op,b,ox) A[0] op b    ox  A[1] op b
#define FI(op,b)    A[0] op b[0];    A[1] op b[1]
#define FX(op,b)    A[0] op b;       A[1] op b
#define DO(macro)   macro(0); macro(1)

#else

#define CI(op,b,ox) A[0] op b[0] ox  A[1] op b[1] ox A[2] op b[2]
#define CX(op,b,ox) A[0] op b    ox  A[1] op b    ox A[2] op b
#define FI(op,b)    A[0] op b[0];    A[1] op b[1];   A[2] op b[2]
#define FX(op,b)    A[0] op b;       A[1] op b;      A[2] op b
#define DO(macro)   macro(0); macro(1); macro(2)

#endif
#define falcON_NSYM1 falcON_NDIM
  //////////////////////////////////////////////////////////////////////////////
  template<typename SCALAR>                    // scalar type                   
  class sym1 {
  public:
    static const int NDAT = falcON_NSYM1;
    //--------------------------------------------------------------------------
    // types                                                                    
    //--------------------------------------------------------------------------
    typedef SCALAR       element_type;         // type for public use           
  private:
    typedef element_type         scal;         // type for const returns        
    typedef sym1                 tsy1;         // type for non-const returns    
    typedef vector<scal>         vect;         // type for const returns        
    typedef int          const&c_indx;         // type for const args           
    typedef scal         const&c_scal;         // type for const args           
    typedef tsy1         const&c_tsy1;         // type for const args           
    //--------------------------------------------------------------------------
    // friendships                                                              
    //--------------------------------------------------------------------------
    friend class vector<element_type>;
    friend class sym2<element_type>;
    friend class sym3<element_type>;
    friend class sym4<element_type>;
    //--------------------------------------------------------------------------
    // data                                                                     
    //--------------------------------------------------------------------------
  private:
    element_type *A;
    //--------------------------------------------------------------------------
    // construction                                                             
    //--------------------------------------------------------------------------
  public:
    sym1 (element_type* x) : A(x) {}
    //--------------------------------------------------------------------------
    // non-const methods                                                        
    //--------------------------------------------------------------------------
       scal* pointer()                   { return A; }
       scal &OP[]   (c_indx i)           { return A[i]; }
       tsy1 &OP=    (c_scal x)           { FX(=,x); return *this; }
       tsy1 &OP=    (const scal* x)      { FI(=,x); return *this; }
    TV tsy1 &OP=    (c_VECT x)           { FI(=,x.A); return *this; }
    TV tsy1 &OP+=   (c_VECT x)           { FI(+=,x.A); return *this; }
    TV tsy1 &OP-=   (c_VECT x)           { FI(-=,x.A); return *this; }
       tsy1 &OP*=   (c_scal x)           { FX(*=,x); return *this; }
       tsy1 &OP/=   (c_scal x)           { return OP*=(scal(1)/x); }
       tsy1 &negate ()                   { FI(=-,A); return *this; }
    TV tsy1 &add_mul(c_VECT y, c_scal x) { FI(+=,x*y.A); return *this; }
    TV tsy1 &sub_mul(c_VECT y, c_scal x) { FI(-=,x*y.A); return *this; }
    //--------------------------------------------------------------------------
#define JOB(i) A[i] += y.A[i]+y.A[i]
    TV tsy1 &add_twice(c_VECT y)         { DO(JOB); return *this; }
#undef  JOB
#define JOB(i) A[i] -= y.A[i]+y.A[i]
    TV tsy1 &sub_twice(c_VECT y)         { DO(JOB); return *this; }
#undef  JOB
    //--------------------------------------------------------------------------
#define JOB(i) A[i]=f(A[i])
    tsy1 &apply  (scal(*f)(c_scal))   { DO(JOB); return *this; }
#undef  JOB
    //--------------------------------------------------------------------------
#define JOB(i) A[i]=f(A[i],x.A[i])
    TV tsy1 &connect(c_VECT x, scal(*f)(c_scal, c_scal))
    { DO(JOB); return *this; }
#undef  JOB
    //--------------------------------------------------------------------------
    // const methods                                                            
    //--------------------------------------------------------------------------
    const scal* const_pointer()     const { return A; }
    scal const&OP[](c_indx i) const { return A[i]; }
       vect OP*    (c_scal x) const { register vect y=*this; y*=x; return y; }
    friend
       vect OP*    (c_scal x,
		    c_tsy1 y)       { return y*x; }
       vect OP/    (c_scal x) const { register vect y(*this); y*=1/x;return y; }
    TV vect OP+    (c_VECT x) const { register vect y(*this); y+=x; return y; }
    TV vect OP-    (c_VECT x) const { register vect y(*this); y-=x; return y; }
       vect OP-    ()         const { register vect y(0); y-= *this; return y; }
       scal norm   ()         const { return CI(*,A,+); }
    friend
       scal norm   (c_tsy1 x)       { return x.norm(); }
    TV scal OP*    (c_VECT x) const { return CI(*,x.A,+); }
    //--------------------------------------------------------------------------
#if falcON_NDIM == 2
    scal maxnorm   ()         const { max(abs(A[0]),abs(A[1])); }
    scal min       ()         const { return nbdy::min(A[0],A[1]); }
    scal max       ()         const { return nbdy::max(A[0],A[1]); }
    TV scal dist_sq(c_VECT x) const { return square(A[0]-x.A[0])
				           + square(A[1]-x.A[1]); }
    TV scal sum_sq (c_VECT x) const { return square(A[0]+x.A[0])
				           + square(A[1]+x.A[1]); }
    TV scal OP^    (c_VECT x) const { return A[0]*x.A[1] - A[1]*x.A[0]; }
#else
    //==========================================================================
    scal maxnorm   ()         const { return nbdy::max(abs(A[0]),abs(A[1]),
						       abs(A[2])); }
    scal min       ()         const { return nbdy::min(A[0],A[1],A[2]); }
    scal max       ()         const { return nbdy::max(A[0],A[1],A[2]); }
    TV scal dist_sq(c_VECT x) const { return square(A[0]-x.A[0])
				           + square(A[1]-x.A[1])
				           + square(A[2]-x.A[2]); }
    TV scal sum_sq (c_VECT x) const { return square(A[0]+x.A[0])
				           + square(A[1]+x.A[1])
				           + square(A[2]+x.A[2]); }
    TV vect OP^    (c_VECT x) const {
      register vect z;
      z.A[0] = A[1]*x.A[2] - A[2]*x.A[1];
      z.A[1] = A[2]*x.A[0] - A[0]*x.A[2];
      z.A[2] = A[0]*x.A[1] - A[1]*x.A[0];
      return z;
    }
#endif
    TV scal dist(c_VECT x) const { return sqrt(dist_sq(x)); }
    //--------------------------------------------------------------------------
    TV friend scal dist_sq(c_tsy1 x, c_VECT y) { return x.dist_sq(y); }
    TV friend scal sum_sq (c_tsy1 x, c_VECT y) { return x.sum_sq(y); }
    TV friend scal dist   (c_tsy1 x, c_VECT y) { return x.dist(y); }
    //--------------------------------------------------------------------------
       bool OP==   (c_scal x) const { return CX(==,x,&&); }
    TV bool OP==   (c_VECT x) const { return CI(==,x.A,&&); }
       bool OP!=   (c_scal x) const { return CX(!=,x,||); }
    TV bool OP!=   (c_VECT x) const { return CI(!=,x.A,||); }
    //--------------------------------------------------------------------------
    // ascii I/O                                                                
    //--------------------------------------------------------------------------
    friend std::ostream& OP<< (std::ostream&s, c_tsy1 x)
    {
      return write_array(s,x.A,NDAT);
    }
    //--------------------------------------------------------------------------
    friend std::istream& OP>> (std::istream&s, tsy1 &x)
    {
      return read_array(s,x.A,NDAT);
    }
  };
  //////////////////////////////////////////////////////////////////////////////
#undef DO
#undef CI
#undef CX
#undef FI
#undef FX
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class nbdy::sym2                                                         //
  //                                                                          //
  //--------------------------------------------------------------------------//
  // a symmetric 3x3 matrix                                                   //
  // element mapping is as follows:   0 -> 0,0                                //
  //                                  1 -> 0,1 & 1,0                          //
  //                                  2 -> 0,2 & 2,0                          //
  //                                  3 -> 1,1                                //
  //                                  4 -> 1,2 & 2,1                          //
  //                                  5 -> 2,2 & 2,2                          //
  //--------------------------------------------------------------------------//
  // or a symmetric 2x2 matrix                                                //
  // element mapping is as follows:   0 -> 0,0                                //
  //                                  1 -> 0,1 & 1,0                          //
  //                                  2 -> 1,1                                //
  //////////////////////////////////////////////////////////////////////////////
#if falcON_NDIM == 2

#define falcON_NSYM2 3
#define CI(op,b,ox) A[0] op b[0] ox A[1] op b[1] ox A[2] op b[2]
#define CX(op,b,ox) A[0] op b    ox A[1] op b    ox A[2] op b
#define FI(op,b)    A[0] op b[0];   A[1] op b[1];   A[2] op b[2]
#define FX(op,b)    A[0] op b;      A[1] op b;      A[2] op b

#else

#define falcON_NSYM2 6
#define CI(op,b,ox) A[0] op b[0] ox A[1] op b[1] ox A[2] op b[2] ox	\
                    A[3] op b[3] ox A[4] op b[4] ox A[5] op b[5]
#define CX(op,b,ox) A[0] op b    ox A[1] op b    ox A[2] op b    ox	\
                    A[3] op b    ox A[4] op b    ox A[5] op b
#define FI(op,b)    A[0] op b[0];   A[1] op b[1];   A[2] op b[2];	\
                    A[3] op b[3];   A[4] op b[4];   A[5] op b[5]
#define FX(op,b)    A[0] op b;      A[1] op b;      A[2] op b;		\
                    A[3] op b;      A[4] op b;      A[5] op b

#endif
  //////////////////////////////////////////////////////////////////////////////
  template<typename SCALAR>                    // scalar type                   
  class sym2 {
  public:
    static const int NDAT = falcON_NSYM2;
    //--------------------------------------------------------------------------
    // types                                                                    
    //--------------------------------------------------------------------------
    typedef SCALAR          element_type;      // type for public use           
  protected:
    typedef element_type            scal;      // type for const returns        
    typedef vector<scal>            vect;      // type for const returns        
    typedef sym1  <scal>            tsy1;      // type for non-const args       
    typedef sym2                    tsy2;      // type for non-const args       
    typedef int             const&c_indx;      // type for const args           
    typedef scal            const&c_scal;      // type for const args           
    typedef tsy1            const&c_tsy1;      // type for const args           
    typedef tsy2            const&c_tsy2;      // type for const args           
    //--------------------------------------------------------------------------
    // friendships                                                              
    //--------------------------------------------------------------------------
    friend class sym3<element_type>;
    friend class sym4<element_type>;
    //--------------------------------------------------------------------------
    // data                                                                     
    //--------------------------------------------------------------------------
  private:
    element_type *A;
    //--------------------------------------------------------------------------
    // construction                                                             
    //--------------------------------------------------------------------------
  private:
    sym2 ();                                   // not implemented               
  public:
    sym2 (element_type* x) : A(x) {}
    //--------------------------------------------------------------------------
    // non-const methods (operator=() to be defined in derived class)           
    //--------------------------------------------------------------------------
    // simple element wise operations, most are needed by falcON                
    scal &OP[]   (c_indx i)           { return A[i]; }
    tsy2 &OP=    (c_scal x)           { FX(=,x); return *this; }
    tsy2 &OP=    (c_tsy2 x)           { FI(=,x.A); return *this; }
    tsy2 &OP+=   (c_tsy2 x)           { FI(+=,x.A); return *this; }
    tsy2 &OP-=   (c_tsy2 x)           { FI(-=,x.A); return *this; }
    tsy2 &OP*=   (c_scal x)           { FX(*=,x); return *this; }
    tsy2 &OP/=   (c_scal x)           { return OP*=(scal(1)/x); }
    tsy2 &negate ()                   { FI(=-,A); return *this; }
    tsy2 &add_mul(c_tsy2 y, c_scal x) { FI(+=,x*y.A); return *this; }
    tsy2 &sub_mul(c_tsy2 y, c_scal x) { FI(-=,x*y.A); return *this; }
    //--------------------------------------------------------------------------
    // unity tensor                                                             
    //   not needed ?                                                           
#if falcON_NDIM==2
    tsy2 &unity ()         { A[0]=A[2]=1; A[1]=0; return *this; }
    tsy2 &unity (c_scal x) { A[0]=A[2]=x; A[1]=0; return *this; }
    tsy2 &add_unity(c_scal x) { A[0]+=x; A[2]+=x; return *this; }
#else
    tsy2 &unity ()         { A[0]=A[3]=A[5]=1; A[1]=A[2]=A[4]=0; return *this; }
    tsy2 &unity (c_scal x) { A[0]=A[3]=A[5]=x; A[1]=A[2]=A[4]=0; return *this; }
    tsy2 &add_unity(c_scal x) { A[0]+=x; A[3]+=x; A[5]+=x; return *this; }
#endif
    //--------------------------------------------------------------------------
    // A_ij = f(A_ij)                                                           
    //   not needed in falcON                                                   
    tsy2 &apply  (scal(*f)(c_scal))   { 
      A[0]=f(A[0]); A[1]=f(A[1]); A[2]=f(A[2]);
#if falcON_NDIM==3
      A[3]=f(A[3]); A[4]=f(A[4]); A[5]=f(A[5]);
#endif
      return *this;
    }
    //--------------------------------------------------------------------------
    // A_ij = f(A_ij, B_ij)                                                     
    //   not needed in falcON                                                   
    tsy2 &connect(c_tsy2 x, scal(*f)(c_scal, c_scal)) {
      A[0]=f(A[0],x.A[0]); A[1]=f(A[1],x.A[1]); A[2]=f(A[2],x.A[2]);
#if falcON_NDIM==3
      A[3]=f(A[3],x.A[3]); A[4]=f(A[4],x.A[4]); A[5]=f(A[5],x.A[5]);
#endif
      return *this;
    }
    //--------------------------------------------------------------------------
    // outer product and related                                                
    //   Aij = Xi*Xj                                                            
    TV tsy2 &outer_prod(c_VECT x) {
      A[0] = x.A[0]*x.A[0];
      A[1] = x.A[0]*x.A[1];
#if falcON_NDIM==2
      A[2] = x.A[1]*x.A[1];
#else
      A[2] = x.A[0]*x.A[2]; 
      A[3] = x.A[1]*x.A[1];
      A[4] = x.A[1]*x.A[2];
      A[5] = x.A[2]*x.A[2];
#endif
      return *this;
    }
    // Aij = Xi*Xj * f                                                          
    TV tsy2 &outer_prod(c_VECT x, c_scal f) {
      register VECTOR y(x); y*=f;
      A[0] = x.A[0]*y.A[0];
      A[1] = x.A[0]*y.A[1];
#if falcON_NDIM==2
      A[2] = x.A[1]*y.A[1];
#else
      A[2] = x.A[0]*y.A[2]; 
      A[3] = x.A[1]*y.A[1];
      A[4] = x.A[1]*y.A[2];
      A[5] = x.A[2]*y.A[2];
#endif
      return *this;
    }
    // Aij += Xi*Xj                                                             
    TV tsy2 &add_outer_prod(c_VECT x) {
      A[0]+= x.A[0]*x.A[0];
      A[1]+= x.A[0]*x.A[1];
#if falcON_NDIM==2
      A[2]+= x.A[1]*x.A[1];
#else
      A[2]+= x.A[0]*x.A[2]; 
      A[3]+= x.A[1]*x.A[1];
      A[4]+= x.A[1]*x.A[2];
      A[5]+= x.A[2]*x.A[2];
#endif
      return *this;
    }
    // Aij += Xi*Xj * f                                                         
    TV tsy2 &add_outer_prod(c_VECT x, c_scal f) {
      register VECTOR y(x); y*=f;
      A[0]+= x.A[0]*y.A[0];
      A[1]+= x.A[0]*y.A[1];
#if falcON_NDIM==2
      A[2]+= x.A[1]*y.A[1];
#else
      A[2]+= x.A[0]*y.A[2]; 
      A[3]+= x.A[1]*y.A[1];
      A[4]+= x.A[1]*y.A[2];
      A[5]+= x.A[2]*y.A[2];
#endif
      return *this;
    }
    //--------------------------------------------------------------------------
    // Aij = a*D_ij + b*Xi*Xj                                                   
    //   needed in cell-cell interaction, file grav.cc                          
    TV tsy2 &quadrup(c_VECT x, c_scal a, c_scal b) {
      register scal tmp=b*x.A[0];
      A[0]= a+tmp*x.A[0];
      A[1]=   tmp*x.A[1];
#if falcON_NDIM==2
      A[2]= a+b*x.A[1]*x.A[1];
#else
      A[2]=   tmp*x.A[2];
      tmp =     b*x.A[1];
      A[3]= a+tmp*x.A[1];
      A[4]=   tmp*x.A[2];
      A[5]= a+b*x.A[2]*x.A[2];
#endif
      return *this;
    }      
    //--------------------------------------------------------------------------
    // Aij += a*D_ij + b*Xi*Xj                                                  
    //   needed in cell-body interaction, file grav.cc                          
    TV tsy2 &add_quadrup(c_VECT x, c_scal a, c_scal b) {
      register scal tmp=b*x.A[0];
      A[0]+= a+tmp*x.A[0];
      A[1]+=   tmp*x.A[1];
#if falcON_NDIM==2
      A[2]+= a+b*x.A[1]*x.A[1];
#else
      A[2]+=   tmp*x.A[2];
      tmp  =     b*x.A[1];
      A[3]+= a+tmp*x.A[1];
      A[4]+=   tmp*x.A[2];
      A[5]+= a+b*x.A[2]*x.A[2];
#endif
      return *this;
    }      
    //--------------------------------------------------------------------------
    // const methods                                                            
    //--------------------------------------------------------------------------
    // simple element wise operations, most are not needed by falcON            
    scal const&OP[](c_indx i) const { return A[i]; }
    bool OP==      (c_scal x) const { return CX(==,x,&&); }
    bool OP==      (c_tsy2 x) const { return CI(==,x.A,&&); }
    bool OP!=      (c_scal x) const { return CX(!=,x,||); }
    bool OP!=      (c_tsy2 x) const { return CI(!=,x.A,||); }
    //--------------------------------------------------------------------------
    // trace = A_ii                                                             
    //   needed in cell-cell interaction in grav.cc                             
    scal trace() const {
#if falcON_NDIM==2
      return A[0]+A[2];
#else
      return A[0]+A[3]+A[5];
#endif
    }
    friend scal trace(c_tsy2 x) { return x.trace(); }
    //--------------------------------------------------------------------------
    // double inner product:  result = Xj*Xi*Aij                                
    //   not needed by falcON ?                                                 
    TV scal contract(c_VECT x) const {
#if falcON_NDIM==2
      return x.A[0]*(A[0]*x.A[0]+twice(A[1]*x.A[1]))
	    +x.A[1]* A[2]*x.A[1];
#else
      return x.A[0]*(x.A[0]*A[0]+twice(x.A[1]*A[1]+x.A[2]*A[2]))
	    +x.A[1]*(x.A[1]*A[3]+twice(x.A[2]*A[4]))
	    +x.A[2]* x.A[2]*A[5];
#endif
    }  
    //------------------------------------------------------------+-------------
    // inner product with vector:  result_i = Xj*Aij              | 9* 6+       
    //    needed several times in grav.cc & tayl.cc                             
    template<typename VECA, typename VECB>
    void inner_prod(const VECA&x, VECB&z) const {
#if falcON_NDIM==2
      z.A[0] = A[0]*x.A[0] + A[1]*x.A[1];
      z.A[1] = A[1]*x.A[0] + A[2]*x.A[1];
#else
      z.A[0] = A[0]*x.A[0] + A[1]*x.A[1] + A[2]*x.A[2];
      z.A[1] = A[1]*x.A[0] + A[3]*x.A[1] + A[4]*x.A[2];
      z.A[2] = A[2]*x.A[0] + A[4]*x.A[1] + A[5]*x.A[2];
#endif
    }
    //--------------------------------------------------------------------------
    // inner product with matrix:  result = Aij*Xij                             
    //   needed ?                                                               
    scal OP* (c_tsy2 x) const {
#if falcON_NDIM==2
      return twice(A[1]*x.A[1])+A[1]*x.A[1]; 
#else
      return twice(A[1]*x.A[1] + A[2]*x.A[2] + A[4]*x.A[4])
                 + A[0]*x.A[0] + A[3]*x.A[3] + A[5]*x.A[5];
#endif
    }
    //--------------------------------------------------------------------------
    // ascii I/O                                                                
    //--------------------------------------------------------------------------
    friend std::ostream& OP<< (std::ostream&s, c_tsy2 x)
    {
      return write_array(s,x.A,NDAT);
    }
    //--------------------------------------------------------------------------
    friend std::istream& OP>> (std::istream&s, tsy2 &x)
    {
      return read_array(s,x.A,NDAT);
    }
  };
  //////////////////////////////////////////////////////////////////////////////
#undef CI
#undef CX
#undef FI
#undef FX
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class nbdy::sym3                                                         //
  //                                                                          //
  //--------------------------------------------------------------------------//
  // a symmetric 3x3x3 tensor                                                 //
  // element mapping is as follows:   0 -> 0,0,0                              //
  //                                  1 -> 0,0,1 (3 permutations)             //
  //                                  2 -> 0,0,2 (3 permutations)             //
  //                                  3 -> 0,1,1 (3 permutations)             //
  //                                  4 -> 0,1,2 (6 permutations)             //
  //                                  5 -> 0,2,2 (3 permutations)             //
  //                                  6 -> 1,1,1                              //
  //                                  7 -> 1,1,2 (3 permutations)             //
  //                                  8 -> 1,2,2 (3 permutations)             //
  //                                  9 -> 2,2,2                              //
  //--------------------------------------------------------------------------//
  // or a symmetric 2x2x2 tensor                                              //
  // element mapping is as follows:   0 -> 0,0,0                              //
  //                                  1 -> 0,0,1 (3 permutations)             //
  //                                  2 -> 0,1,1 (3 permutations)             //
  //                                  3 -> 1,1,1                              //
  //////////////////////////////////////////////////////////////////////////////
#if falcON_NDIM == 2

#define falcON_NSYM3 4
#define CI(op,b,ox) A[0] op b[0] ox A[1] op b[1] ox A[2] op b[2] ox A[3] op b[3]
#define CX(op,b,ox) A[0] op b    ox A[1] op b    ox A[2] op b    ox A[3] op b
#define FI(op,b)    A[0] op b[0];   A[1] op b[1];   A[2] op b[2];   A[3] op b[3]
#define FX(op,b)    A[0] op b;      A[1] op b;      A[2] op b;      A[3] op b

#else

#define falcON_NSYM3 10
#define CI(op,b,ox)							\
     A[0] op b[0] ox A[1] op b[1] ox A[2] op b[2] ox A[3] op b[3]	\
  ox A[4] op b[4] ox A[5] op b[5] ox A[6] op b[6] ox A[7] op b[7]	\
  ox A[8] op b[8] ox A[9] op b[9]
#define CX(op,b,ox)							\
     A[0] op b    ox A[1] op b    ox A[2] op b    ox A[3] op b		\
  ox A[4] op b    ox A[5] op b    ox A[6] op b    ox A[7] op b		\
  ox A[8] op b    ox A[9] op b
#define FI(op,b)							\
     A[0] op b[0];   A[1] op b[1];   A[2] op b[2];   A[3] op b[3];	\
     A[4] op b[4];   A[5] op b[5];   A[6] op b[6];   A[7] op b[7];	\
     A[8] op b[8];   A[9] op b[9]
#define FX(op,b)							\
     A[0] op b;      A[1] op b;      A[2] op b;      A[3] op b;		\
     A[4] op b;      A[5] op b;      A[6] op b;      A[7] op b;		\
     A[8] op b;      A[9] op b

#endif
  //////////////////////////////////////////////////////////////////////////////
  template<typename SCALAR>                    // scalar type                   
  class sym3 {
  public:
    static const int NDAT = falcON_NSYM3;
    //--------------------------------------------------------------------------
    // types                                                                    
    //--------------------------------------------------------------------------
    typedef SCALAR          element_type;      // type for public use           
  protected:
    typedef element_type            scal;      // type for const returns        
    typedef vector<scal>            vect;      // type for const returns        
    typedef sym1  <scal>            tsy1;      // type for non-const args       
    typedef sym2  <scal>            tsy2;      // type for non-const args       
    typedef sym3                    tsy3;      // type for non-const args/return
    typedef int             const&c_indx;      // type for const args           
    typedef scal            const&c_scal;      // type for const args           
    typedef tsy2            const&c_tsy2;      // type for const args           
    typedef tsy3            const&c_tsy3;      // type for const args           
    //--------------------------------------------------------------------------
    // friendships                                                              
    //--------------------------------------------------------------------------
    friend class sym4<element_type>;
    //--------------------------------------------------------------------------
    // data                                                                     
    //--------------------------------------------------------------------------
  private:
    element_type *A;
    //--------------------------------------------------------------------------
    // construction                                                             
    //--------------------------------------------------------------------------
  private:
    sym3 ();                                   // not implemented               
  public:
    sym3 (element_type*x) : A(x) {}
    //--------------------------------------------------------------------------
    // non-const methods                                                        
    //--------------------------------------------------------------------------
  public:
    scal &OP[]   (c_indx i)           { return A[i]; }
    tsy3 &OP=    (c_scal x)           { FX(=,x); return *this; }
    tsy3 &OP=    (c_tsy3 x)           { FI(=,x.A); return *this; }
    tsy3 &OP+=   (c_tsy3 x)           { FI(+=,x.A); return *this; }
    tsy3 &OP-=   (c_tsy3 x)           { FI(-=,x.A); return *this; }
    tsy3 &OP*=   (c_scal x)           { FX(*=,x); return *this; }
    tsy3 &OP/=   (c_scal x)           { return OP*=(scal(1)/x); }
    tsy3 &negate ()                   { FI(=-,A); return *this; }
    tsy3 &add_mul(c_tsy3 y, c_scal x) { FI(+=,x*y.A); return *this; }
    tsy3 &sub_mul(c_tsy3 y, c_scal x) { FI(-=,x*y.A); return *this; }
    //--------------------------------------------------------------------------
    // unity tensor                                                             
    //   not needed ?                                                           
#if falcON_NDIM==2
    tsy3 &unity  ()           { A[0]=A[3]=1; A[1]=A[2]=0; return *this; }
    tsy3 &unity  (c_scal x)   { A[0]=A[3]=x; A[1]=A[2]=0; return *this; }
    tsy3 &add_unity(c_scal x) { A[0]+=x; A[3]+=x;         return *this; }
#else
    //--------------------------------------------------------------------------
    tsy3 &unity  () {
      A[1]=A[2]=A[3]=A[4]=A[5]=A[7]=A[8]=0;
      A[0]=A[6]=A[9]=1;
      return *this;
    }
    tsy3 &unity  (c_scal x) {
      A[1]=A[2]=A[3]=A[4]=A[5]=A[7]=A[8]=0;
      A[0]=A[6]=A[9]=x;
      return *this;
    }
    tsy3 &add_unity  (c_scal x) {
      A[0]+= x;
      A[6]+= x;
      A[9]+= x;
      return *this;
    }
#endif
    //--------------------------------------------------------------------------
    // A_ijk = f(A_ijk)                                                         
    //   not needed in falcON                                                   
    tsy3 &apply  (scal(*f)(c_scal))   { 
      A[0]=f(A[0]); A[1]=f(A[1]); A[2]=f(A[2]); A[3]=f(A[3]); 
#if falcON_NDIM==3
      A[4]=f(A[4]); A[5]=f(A[5]); A[6]=f(A[6]); A[7]=f(A[7]);
      A[8]=f(A[8]); A[9]=f(A[9]);
#endif
      return *this;
    }
    //--------------------------------------------------------------------------
    // A_ijk = f(A_ijk, B_ijk)                                                  
    //   not needed in falcON                                                   
    tsy3 &connect(c_tsy3 x, scal(*f)(c_scal, c_scal)) {
      A[0]=f(A[0],x.A[0]);
      A[1]=f(A[1],x.A[1]);
      A[2]=f(A[2],x.A[2]);
      A[3]=f(A[3],x.A[3]); 
#if falcON_NDIM==3
      A[4]=f(A[4],x.A[4]);
      A[6]=f(A[5],x.A[6]);
      A[7]=f(A[7],x.A[7]);
      A[8]=f(A[8],x.A[8]);
      A[9]=f(A[9],x.A[9]);
#endif
      return *this;
    }
    //--------------------------------------------------------------------------
    // outer product Aijk = Xi*Xj*Xk, where Xi*Xj is given as argument, too     
    TV tsy3 &outer_prod(c_tsy2 XX, c_VECT X) {
#if falcON_NDIM==2
      A[0] = XX.A[0] * X.A[0];
      A[1] = XX.A[0] * X.A[1];
      A[2] = XX.A[1] * X.A[1];
      A[3] = XX.A[2] * X.A[1];
#else
      A[0] = XX.A[0] * X.A[0];
      A[1] = XX.A[0] * X.A[1];
      A[2] = XX.A[0] * X.A[2];
      A[3] = XX.A[1] * X.A[1];
      A[4] = XX.A[1] * X.A[2];
      A[5] = XX.A[2] * X.A[2];
      A[6] = XX.A[3] * X.A[1];
      A[7] = XX.A[3] * X.A[2];
      A[8] = XX.A[4] * X.A[2];
      A[9] = XX.A[5] * X.A[2];
#endif
      return *this;
    }
    //--------------------------------------------------------------------------
    TV tsy3 &add_outer_prod(c_tsy2 XX, c_VECT X) {
#if falcON_NDIM==2
      A[0] += XX.A[0] * X.A[0];
      A[1] += XX.A[0] * X.A[1];
      A[2] += XX.A[1] * X.A[1];
      A[3] += XX.A[2] * X.A[1];
#else
      A[0] += XX.A[0] * X.A[0];
      A[1] += XX.A[0] * X.A[1];
      A[2] += XX.A[0] * X.A[2];
      A[3] += XX.A[1] * X.A[1];
      A[4] += XX.A[1] * X.A[2];
      A[5] += XX.A[2] * X.A[2];
      A[6] += XX.A[3] * X.A[1];
      A[7] += XX.A[3] * X.A[2];
      A[8] += XX.A[4] * X.A[2];
      A[9] += XX.A[5] * X.A[2];
#endif
      return *this;
    }
    //--------------------------------------------------------------------------
    // outer symmetric product Aij = Xi*Yjk + Xj*Yik + Xk*Yij                   
    TV tsy3 &symm_prod(c_tsy2 y, c_VECT x) {
      A[0] = trice(y.A[0]*x.A[0]);
      A[1] = twice(y.A[1]*x.A[0]) + y.A[0]*x.A[1];
#if falcON_NDIM==2
      A[2] = twice(y.A[1]*x.A[1]) + y.A[2]*x.A[0];
      A[3] = trice(y.A[2]*x.A[1]);
#else
      A[6] = trice(y.A[3]*x.A[1]);
      A[9] = trice(y.A[5]*x.A[2]);
      A[2] = twice(y.A[2]*x.A[0]) + y.A[0]*x.A[2];
      A[3] = twice(y.A[1]*x.A[1]) + y.A[3]*x.A[0];
      A[5] = twice(y.A[2]*x.A[2]) + y.A[5]*x.A[0];
      A[7] = twice(y.A[4]*x.A[1]) + y.A[3]*x.A[2];
      A[8] = twice(y.A[4]*x.A[2]) + y.A[5]*x.A[1];
      A[4] = y.A[1]*x.A[2] +y.A[2]*x.A[1] +y.A[4]*x.A[0];
#endif
      return *this;
    }
    //--------------------------------------------------------------------------
    // outer symmetric product Aij = Xi*Yjk + Xj*Yik + Xk*Yij                   
    TV tsy3 &add_symm_prod(c_tsy2 y, c_VECT x) {
      A[0]+= trice(y.A[0]*x.A[0]);
      A[1]+= twice(y.A[1]*x.A[0]) + y.A[0]*x.A[1];
#if falcON_NDIM==2
      A[2]+= twice(y.A[1]*x.A[1]) + y.A[2]*x.A[0];
      A[3]+= trice(y.A[2]*x.A[1]);
#else
      A[6]+= trice(y.A[3]*x.A[1]);
      A[9]+= trice(y.A[5]*x.A[2]);
      A[2]+= twice(y.A[2]*x.A[0]) + y.A[0]*x.A[2];
      A[3]+= twice(y.A[1]*x.A[1]) + y.A[3]*x.A[0];
      A[5]+= twice(y.A[2]*x.A[2]) + y.A[5]*x.A[0];
      A[7]+= twice(y.A[4]*x.A[1]) + y.A[3]*x.A[2];
      A[8]+= twice(y.A[4]*x.A[2]) + y.A[5]*x.A[1];
      A[4]+= y.A[1]*x.A[2] +y.A[2]*x.A[1] +y.A[4]*x.A[0];
#endif
      return *this;
    }
    //--------------------------------------------------------------------------
    // Aijk = a * (Dij*Xk + Dik*Xj + Djk*Xi) + b * Xi*Xj*Xk                     
    //   needed in cell-cell interaction, file grav.cc                          
    TV tsy3 &quadrup(c_VECT x, c_scal a, c_scal b) {
      register scal a3=a+a+a, tmp=b*x.A[0]*x.A[0];
      A[0] = (a3+tmp)*x.A[0];
      A[1] = (a +tmp)*x.A[1];
#if falcON_NDIM==2
      tmp  = b*x.A[1]*x.A[1];
      A[2] = (a +tmp)*x.A[0];
      A[3] = (a3+tmp)*x.A[1];
#else
      A[2] = (a +tmp)*x.A[2];
      tmp  = b*x.A[1]*x.A[1];
      A[3] = (a +tmp)*x.A[0];
      A[6] = (a3+tmp)*x.A[1];
      A[7] = (a +tmp)*x.A[2];
      tmp  = b*x.A[2]*x.A[2];
      A[5] = (a +tmp)*x.A[0];
      A[8] = (a +tmp)*x.A[1];
      A[9] = (a3+tmp)*x.A[2];
      A[4] = b*x.A[0]*x.A[1]*x.A[2];
#endif
      return *this;
    }      
    //--------------------------------------------------------------------------
    // Aijk += a * (Dij*Xk + Dik*Xj + Djk*Xi) + b * Xi*Xj*Xk                    
    //   needed in cell-cell interaction, file grav.cc                          
    TV tsy3 &add_quadrup(c_VECT x, c_scal a, c_scal b) {
      register scal a3=a+a+a, tmp=b*x.A[0]*x.A[0];
      A[0]+= (a3+tmp)*x.A[0];
      A[1]+= (a +tmp)*x.A[1];
#if falcON_NDIM==2
      tmp  = b*x.A[1]*x.A[1];
      A[2]+= (a +tmp)*x.A[0];
      A[3]+= (a3+tmp)*x.A[1];
#else
      A[2]+= (a +tmp)*x.A[2];
      tmp  = b*x.A[1]*x.A[1];
      A[3]+= (a +tmp)*x.A[0];
      A[6]+= (a3+tmp)*x.A[1];
      A[7]+= (a +tmp)*x.A[2];
      tmp  = b*x.A[2]*x.A[2];
      A[5]+= (a +tmp)*x.A[0];
      A[8]+= (a +tmp)*x.A[1];
      A[9]+= (a3+tmp)*x.A[2];
      A[4]+= b*x.A[0]*x.A[1]*x.A[2];
#endif
      return *this;
    }      
    //--------------------------------------------------------------------------
    // const methods                                                            
    //--------------------------------------------------------------------------
    // simple element wise operations, most are not needed by falcON            
    scal const&OP[](c_indx i) const { return A[i]; }
    bool OP==      (c_scal x) const { return CX(==,x,&&); }
    bool OP==      (c_tsy3 x) const { return CI(==,x.A,&&); }
    bool OP!=      (c_scal x) const { return CX(!=,x,||); }
    bool OP!=      (c_tsy3 x) const { return CI(!=,x.A,||); }
    //--------------------------------------------------------------------------
    // return_k = Aijk * D_ij                                                   
    //   not currently needed by falcON                                         
    vect trace() const {
      register vect w;
#if falcON_NDIM==2
      w.A[0] = A[0]+A[2];
      w.A[1] = A[1]+A[3];
#else
      w.A[0] = A[0]+A[3]+A[5];
      w.A[1] = A[1]+A[6]+A[8];
      w.A[2] = A[3]+A[7]+A[9];
#endif
      return w;
    }
    //--------------------------------------------------------------------------
    // trace(x*A) = Dij * xk * Aijk                                             
    //   not currently needed by falcON                                         
    TV scal trace(c_VECT x) const {
#if falcON_NDIM==2
      return   x.A[0]*(A[0]+A[2]) + x.A[1]*(A[1]+A[3]);
#else
      return   x.A[0]*(A[0]+A[3]+A[5])
	     + x.A[1]*(A[1]+A[6]+A[8])
	     + x.A[2]*(A[2]+A[7]+A[9]);
#endif
    }
    //------------------------------------------------------------+-------------
    // inner product W_k = A_ijk*Xij                              | 18*  18+    
    TV void inner_prod(c_tsy2 x, VECTOR&w) const {
#if falcON_NDIM==2
      w.A[0] = A[0]*x.A[0] + A[2]*x.A[2] + twice(A[1]*x.A[1]);
      w.A[1] = A[1]*x.A[0] + A[3]*x.A[2] + twice(A[2]*x.A[1]);
#else
      w.A[0] = A[0]*x.A[0] +A[3]*x.A[3] +A[5]*x.A[5]
	+twice(A[1]*x.A[1] +A[2]*x.A[2] +A[4]*x.A[4]);
      w.A[1] = A[1]*x.A[0] +A[6]*x.A[3] +A[8]*x.A[5]
	+twice(A[3]*x.A[1] +A[4]*x.A[2] +A[7]*x.A[4]);
      w.A[2] = A[2]*x.A[0] +A[7]*x.A[3] +A[9]*x.A[5]
	+twice(A[4]*x.A[1] +A[5]*x.A[2] +A[8]*x.A[4]);
#endif
    }
    //--------------------------------------------------------------------------
    // inner product W_ij = A_ijk*Xk                                            
    TV void inner_prod(c_VECT x, tsy2& w) const {
#if falcON_NDIM==2
      w.A[0] = A[0]*x.A[0] +A[1]*x.A[1];
      w.A[1] = A[1]*x.A[0] +A[2]*x.A[1];
      w.A[2] = A[2]*x.A[0] +A[3]*x.A[1];
#else
      w.A[0] = A[0]*x.A[0] + A[1]*x.A[1] + A[2]*x.A[2];
      w.A[1] = A[1]*x.A[0] + A[3]*x.A[1] + A[4]*x.A[2];
      w.A[2] = A[2]*x.A[0] + A[4]*x.A[1] + A[5]*x.A[2];
      w.A[3] = A[3]*x.A[0] + A[6]*x.A[1] + A[7]*x.A[2];
      w.A[4] = A[4]*x.A[0] + A[7]*x.A[1] + A[8]*x.A[2];
      w.A[5] = A[5]*x.A[0] + A[8]*x.A[1] + A[9]*x.A[2];
#endif
    }
    //--------------------------------------------------------------------------
    // inner product X= A_ijk*B_ijk                                             
    //--------------------------------------------------------------------------
    scal inner_prod(c_tsy3 B) const {
#if falcON_NDIM==2
      return       (B.A[ 0]*A[ 0] +B.A[ 3]*A[ 3])
      + trice(B.A[ 1]*A[ 1] +B.A[ 2]*A[ 2]);
#else
      return       (B.A[ 0]*A[ 0] +B.A[ 6]*A[ 6] +B.A[ 9]*A[ 9])
      + trice(B.A[ 1]*A[ 1] +B.A[ 2]*A[ 2] +B.A[ 3]*A[ 3] +B.A[ 5]*A[ 5] +
                    B.A[ 7]*A[ 7] +B.A[ 8]*A[ 8])
      +        6 * (B.A[ 4]*A[ 4]);
#endif
    }
    //--------------------------------------------------------------------------
    // ascii I/O                                                                
    //--------------------------------------------------------------------------
    friend std::ostream& OP<< (std::ostream&s, c_tsy3 x)
    {
      return write_array(s,x.A,NDAT);
    }
    //--------------------------------------------------------------------------
    friend std::istream& OP>> (std::istream&s, tsy3 &x)
    {
      return read_array(s,x.A,NDAT);
    }
  };
  //////////////////////////////////////////////////////////////////////////////
#undef CI
#undef CX
#undef FI
#undef FX
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class nbdy::sym4                                                         //
  //                                                                          //
  //--------------------------------------------------------------------------//
  // a symmetric 3x3x3x3 tensor                                               //
  // element mapping is as follows:   0 -> 0,0,0,0                            //
  //                                  1 -> 0,0,0,1 ( 4 permutations)          //
  //                                  2 -> 0,0,0,2 ( 4 permutations)          //
  //                                  3 -> 0,0,1,1 ( 6 permutations)          //
  //                                  4 -> 0,0,1,2 (12 permutations)          //
  //                                  5 -> 0,0,2,2 ( 6 permutations)          //
  //                                  6 -> 0,1,1,1 ( 4 permutations)          //
  //                                  7 -> 0,1,1,2 (12 permutations)          //
  //                                  8 -> 0,1,2,2 (12 permutations)          //
  //                                  9 -> 0,2,2,2 ( 4 permutations)          //
  //                                 10 -> 1,1,1,1                            //
  //                                 11 -> 1,1,1,2 ( 4 permutations)          //
  //                                 12 -> 1,1,2,2 ( 6 permutations)          //
  //                                 13 -> 1,2,2,2 ( 4 permutations)          //
  //                                 14 -> 2,2,2,2                            //
  //--------------------------------------------------------------------------//
  // or a symmetric 2x2x2x2 tensor                                            //
  // element mapping is as follows:   0 -> 0,0,0,0                            //
  //                                  1 -> 0,0,0,1 & permutations             //
  //                                  2 -> 0,0,1,1 & permutations             //
  //                                  3 -> 0,1,1,1 & permutations             //
  //                                  4 -> 1,1,1,1                            //
  //////////////////////////////////////////////////////////////////////////////
#if falcON_NDIM == 2

#define falcON_NSYM4 5
#define CI(op,b,ox) A[ 0] op b[ 0] ox A[ 1] op b[ 1] ox A[ 2] op b[ 2]	\
                 ox A[ 3] op b[ 3] ox A[ 4] op b[ 4]
#define CX(op,b,ox) A[ 0] op b     ox A[ 1] op b     ox A[ 2] op b	\
                 ox A[ 3] op b     ox A[ 4] op b
#define FI(op,b)    A[ 0] op b[ 0];   A[ 1] op b[ 1];   A[ 2] op b[ 2];	\
                    A[ 3] op b[ 3];   A[ 4] op b[ 4]
#define FX(op,b)    A[ 0] op b;       A[ 1] op b;       A[ 2] op b;	\
                    A[ 3] op b;       A[ 4] op b

#else

#define falcON_NSYM4 15
#define CI(op,b,ox) A[ 0] op b[ 0] ox A[ 1] op b[ 1] ox A[ 2] op b[ 2]	\
                 ox A[ 3] op b[ 3] ox A[ 4] op b[ 4] ox A[ 5] op b[ 5]	\
		 ox A[ 6] op b[ 6] ox A[ 7] op b[ 7] ox A[ 8] op b[ 8]	\
		 ox A[ 9] op b[ 9] ox A[10] op b[10] ox A[11] op b[11]	\
		 ox A[12] op b[12] ox A[13] op b[13] ox A[14] op b[14]
#define CX(op,b,ox) A[ 0] op b     ox A[ 1] op b     ox A[ 2] op b	\
                 ox A[ 3] op b     ox A[ 4] op b     ox A[ 5] op b	\
		 ox A[ 6] op b     ox A[ 7] op b     ox A[ 8] op b	\
		 ox A[ 9] op b     ox A[10] op b     ox A[11] op b	\
		 ox A[12] op b     ox A[13] op b     ox A[14] op b
#define FI(op,b)    A[ 0] op b[ 0];   A[ 1] op b[ 1];   A[ 2] op b[ 2];	\
                    A[ 3] op b[ 3];   A[ 4] op b[ 4];   A[ 5] op b[ 5];	\
                    A[ 6] op b[ 6];   A[ 7] op b[ 7];   A[ 8] op b[ 8];	\
		    A[ 9] op b[ 9];   A[10] op b[10];   A[11] op b[11];	\
		    A[12] op b[12];   A[13] op b[13];   A[14] op b[14]
#define FX(op,b)    A[ 0] op b;       A[ 1] op b;       A[ 2] op b;	\
                    A[ 3] op b;       A[ 4] op b;       A[ 5] op b;	\
		    A[ 6] op b;       A[ 7] op b;       A[ 8] op b;	\
		    A[ 9] op b;       A[10] op b;       A[11] op b;	\
		    A[12] op b;       A[13] op b;       A[14] op b
#endif
  //////////////////////////////////////////////////////////////////////////////
  template<typename SCALAR>                    // scalar type                   
  class sym4 {
  public:
    static const int NDAT = falcON_NSYM4;
    //--------------------------------------------------------------------------
    // types                                                                    
    //--------------------------------------------------------------------------
    typedef SCALAR          element_type;      // type for public use           
  protected:
    typedef element_type            scal;      // type for const returns        
    typedef vector<scal>            vect;      // type for const returns        
    typedef sym1  <scal>            tsy1;      // type for non-const args       
    typedef sym2  <scal>            tsy2;      // type for non-const args       
    typedef sym3  <scal>            tsy3;      // type for non-const args       
    typedef sym4                    tsy4;      // type for non-const args/return
    typedef int             const&c_indx;      // type for const args           
    typedef scal            const&c_scal;      // type for const args           
    typedef vect            const&c_vect;      // type for const args           
    typedef tsy1            const&c_tsy1;      // type for const args           
    typedef tsy2            const&c_tsy2;      // type for const args           
    typedef tsy3            const&c_tsy3;      // type for const args           
    typedef tsy4            const&c_tsy4;      // type for const args           
    //--------------------------------------------------------------------------
    // data                                                                     
    //--------------------------------------------------------------------------
  private:
    element_type *A;
    //--------------------------------------------------------------------------
    // construction                                                             
    //--------------------------------------------------------------------------
  private:
    sym4 ();                                   // not implemented               
  public:
    sym4 (element_type* x) : A(x) {}
    //--------------------------------------------------------------------------
    // non-const methods                                                        
    //--------------------------------------------------------------------------
  public:
    scal &OP[]   (c_indx i)           { return A[i]; }
    tsy4 &OP=    (c_scal x)           { FX(=,x); return *this; }
    tsy4 &OP=    (c_tsy4 x)           { FI(=,x.A); return *this; }
    tsy4 &OP+=   (c_tsy4 x)           { FI(+=,x.A); return *this; }
    tsy4 &OP-=   (c_tsy4 x)           { FI(-=,x.A); return *this; }
    tsy4 &OP*=   (c_scal x)           { FX(*=,x); return *this; }
    tsy4 &OP/=   (c_scal x)           { return OP*=(scal(1)/x); }
    tsy4 &negate ()                   { FI(=-,A); return *this; }
    tsy4 &add_mul(c_tsy4 y, c_scal x) { FI(+=,x*y.A); return *this; }
    tsy4 &sub_mul(c_tsy4 y, c_scal x) { FI(-=,x*y.A); return *this; }
    //--------------------------------------------------------------------------
    // Aijkl =   a (  Dij*Dkl + Djk*Dil + Dki*Djl)                              
    //         + b (  Dij*Xkl + Djk*Xil + Dki*Xjl                               
    //              + Xij*Dkl + Xjk*Dil + Xki*Djl)                              
    //         + c    Xij*Xkl                                                   
    // where Xij = Xi Xj is given                                               
    tsy4 &quadrup(c_tsy2 x, c_scal a, c_scal b, c_scal c) {
      register scal a3=a+a+a, b3=b+b+b, b6=b3+b3, tmp=c*x.A[0];
      A[ 0] = a3 + x.A[0] * (b6 + tmp); 
      A[ 1] =      x.A[1] * (b3 + tmp);
#if falcON_NDIM==2
      tmp   = c*x.A[2];
      A[ 2] = a  + x.A[0] * (b  + tmp) + b * x.A[2];
      A[ 3] =      x.A[1] * (b3 + tmp);
      A[ 4] =      x.A[2] * (b6 + tmp);
#else
      A[ 2] =      x.A[2] * (b3 + tmp);
      A[ 4] =      x.A[4] * (b  + tmp);
      A[ 5] = a  + x.A[5] * (b  + tmp) + b * x.A[0];
      tmp   = c*x.A[3];
      A[ 3] = a  + x.A[0] * (b  + tmp) + b * x.A[3];
      A[ 6] =      x.A[1] * (b3 + tmp);
      A[ 7] =      x.A[2] * (b  + tmp);
      A[10] =      x.A[3] * (b6 + tmp);
      A[11] =      x.A[4] * (b3 + tmp);
      tmp   = c*x.A[5];
      A[ 8] =      x.A[1] * (b  + tmp);
      A[ 9] = a3 + x.A[2] * (b3 + tmp);
      A[12] = a  + x.A[3] * (b  + tmp) + b * x.A[5];
      A[13] =      x.A[4] * (b3 + tmp);
      A[14] = a3 + x.A[5] * (b6 + tmp);
#endif
      return *this;
    }
    //--------------------------------------------------------------------------
    // Aijkl +=   a (  Dij*Dkl + Djk*Dil + Dki*Djl)                             
    //          + b (  Dij*Xkl + Djk*Xil + Dki*Xjl                              
    //               + Xij*Dkl + Xjk*Dil + Xki*Djl)                             
    //          + c    Xij*Xkl                                                  
    // where Xij = Xi Xj is given                                               
    tsy4 &add_quadrup(c_tsy2 x, c_scal a, c_scal b, c_scal c) {
      register scal a3=a+a+a, b3=b+b+b, b6=b3+b3, tmp=c*x.A[0];
      A[ 0]+= a3 + x.A[0] * (b6 + tmp); 
      A[ 1]+=      x.A[1] * (b3 + tmp);
#if falcON_NDIM==2
      tmp   = c*x.A[2];
      A[ 2]+= a  + x.A[0] * (b  + tmp) + b * x.A[2];
      A[ 3]+=      x.A[1] * (b3 + tmp);
      A[ 4]+=      x.A[2] * (b6 + tmp);
#else
      A[ 2]+=      x.A[2] * (b3 + tmp);
      A[ 4]+=      x.A[4] * (b  + tmp);
      A[ 5]+= a  + x.A[5] * (b  + tmp) + b * x.A[0];
      tmp   = c*x.A[3];
      A[ 3]+= a  + x.A[0] * (b  + tmp) + b * x.A[3];
      A[ 6]+=      x.A[1] * (b3 + tmp);
      A[ 7]+=      x.A[2] * (b  + tmp);
      A[10]+=      x.A[3] * (b6 + tmp);
      A[11]+=      x.A[4] * (b3 + tmp);
      tmp   = c*x.A[5];
      A[ 8]+=      x.A[1] * (b  + tmp);
      A[ 9]+= a3 + x.A[2] * (b3 + tmp);
      A[12]+= a  + x.A[3] * (b  + tmp) + b * x.A[5];
      A[13]+=      x.A[4] * (b3 + tmp);
      A[14]+= a3 + x.A[5] * (b6 + tmp);
#endif
      return *this;
    }
    //--------------------------------------------------------------------------
    // Aijkl +=   a (  Dij*Dkl + Djk*Dil + Dki*Djl)                             
    //          + b (  Dij*Xkl + Djk*Xil + Dki*Xjl                              
    //               + Xij*Dkl + Xjk*Dil + Xki*Djl)                             
    //          + c    Xij*Xkl                                                  
    // where Xij = Xi Xj and X is given                                         
    TV tsy4 &add_quadrup(c_VECT x, c_scal a, c_scal b, c_scal c) {
      SYM2(X);
      return add_quadrup(X.outer_prod(x),a,b,c);
    }
    //--------------------------------------------------------------------------
    // outer product Aijkl = Xi*Xj*Xk*Xl, where Xi*Xj*Xk is given, too          
    TV tsy4 &outer_prod(c_tsy3 XXX, c_VECT X) {
#if falcON_NDIM==2
      A[ 0] = XXX.A[0] * X.A[0];
      A[ 1] = XXX.A[0] * X.A[1];
      A[ 2] = XXX.A[1] * X.A[1];
      A[ 3] = XXX.A[2] * X.A[1];
      A[ 4] = XXX.A[3] * X.A[1];
#else
      A[ 0] = XXX.A[0] * X.A[0];
      A[ 1] = XXX.A[0] * X.A[1];
      A[ 2] = XXX.A[0] * X.A[2];
      A[ 3] = XXX.A[1] * X.A[1];
      A[ 4] = XXX.A[1] * X.A[2];
      A[ 5] = XXX.A[2] * X.A[2];
      A[ 6] = XXX.A[3] * X.A[1];
      A[ 7] = XXX.A[3] * X.A[2];
      A[ 8] = XXX.A[4] * X.A[2];
      A[ 9] = XXX.A[5] * X.A[2];
      A[10] = XXX.A[6] * X.A[1];
      A[11] = XXX.A[6] * X.A[2];
      A[12] = XXX.A[7] * X.A[2];
      A[13] = XXX.A[8] * X.A[2];
      A[14] = XXX.A[9] * X.A[2];
#endif
      return *this;
    }
    //--------------------------------------------------------------------------
    TV tsy4 &add_outer_prod(c_tsy3 XXX, c_VECT X) {
#if falcON_NDIM==2
      A[ 0]+= XXX.A[0] * X.A[0];
      A[ 1]+= XXX.A[0] * X.A[1];
      A[ 2]+= XXX.A[1] * X.A[1];
      A[ 3]+= XXX.A[2] * X.A[1];
      A[ 4]+= XXX.A[3] * X.A[1];
#else
      A[ 0]+= XXX.A[0] * X.A[0];
      A[ 1]+= XXX.A[0] * X.A[1];
      A[ 2]+= XXX.A[0] * X.A[2];
      A[ 3]+= XXX.A[1] * X.A[1];
      A[ 4]+= XXX.A[1] * X.A[2];
      A[ 5]+= XXX.A[2] * X.A[2];
      A[ 6]+= XXX.A[3] * X.A[1];
      A[ 7]+= XXX.A[3] * X.A[2];
      A[ 8]+= XXX.A[4] * X.A[2];
      A[ 9]+= XXX.A[5] * X.A[2];
      A[10]+= XXX.A[6] * X.A[1];
      A[11]+= XXX.A[6] * X.A[2];
      A[12]+= XXX.A[7] * X.A[2];
      A[13]+= XXX.A[8] * X.A[2];
      A[14]+= XXX.A[9] * X.A[2];
#endif
      return *this;
    }
    //--------------------------------------------------------------------------
    // const methods                                                            
    //--------------------------------------------------------------------------
    // simple element wise operations, most are not needed by falcON            
    scal const&OP[](c_indx i) const { return A[i]; }
    bool OP==      (c_scal x) const { return CX(==,x,&&); }
    bool OP==      (c_tsy4 x) const { return CI(==,x.A,&&); }
    bool OP!=      (c_scal x) const { return CX(!=,x,||); }
    bool OP!=      (c_tsy4 x) const { return CI(!=,x.A,||); }
    //-----------------------------------------------------------+--------------
    // inner product W_l = A_ijkl*Xijk                           | 33* 11+      
    TV void inner_prod(c_tsy3 x, VECTOR&w) const {
#if falcON_NDIM==2
      w.A[0] = A[ 0]*x.A[0] + A[ 3]*x.A[3]
	+trice(A[ 1]*x.A[1] + A[ 2]*x.A[2]);
      w.A[1] = A[ 1]*x.A[0] + A[ 4]*x.A[3]
	+trice(A[ 2]*x.A[1] + A[ 3]*x.A[2]);
#else
      w.A[0] = A[ 0]*x.A[0] + A[ 6]*x.A[6] + A[ 9]*x.A[9]
	+trice(A[ 1]*x.A[1] + A[ 2]*x.A[2] + A[ 3]*x.A[3] + 
	       A[ 5]*x.A[5] + A[ 7]*x.A[7] + A[ 8]*x.A[8])
	+six * A[ 4]*x.A[4];
      w.A[1] = A[ 1]*x.A[0] + A[10]*x.A[6] + A[13]*x.A[9]
	+trice(A[ 3]*x.A[1] + A[ 4]*x.A[2] + A[ 6]*x.A[3] + 
	       A[ 8]*x.A[5] + A[11]*x.A[7] + A[12]*x.A[8])
	+six * A[ 7]*x.A[4];
      w.A[2] = A[ 2]*x.A[0] + A[11]*x.A[6] + A[14]*x.A[9]
	+trice(A[ 4]*x.A[1] + A[ 5]*x.A[2] + A[ 7]*x.A[3] + 
	       A[ 9]*x.A[5] + A[12]*x.A[7] + A[13]*x.A[8])
	+six * A[ 8]*x.A[4];
#endif
    }
    //--------------------------------------------------------------------------
    // inner product W_ijk = A_ijkl*Xl                                          
    TV void inner_prod(c_VECT x, tsy3& w) const {
#if falcON_NDIM==2
      w.A[0] = A[ 0]*x.A[0] + A[ 1]*x.A[1];
      w.A[1] = A[ 1]*x.A[0] + A[ 2]*x.A[1];
      w.A[2] = A[ 2]*x.A[0] + A[ 3]*x.A[1];
      w.A[3] = A[ 3]*x.A[0] + A[ 4]*x.A[1];
#else
      w.A[0] = A[ 0]*x.A[0] + A[ 1]*x.A[1] + A[ 2]*x.A[2];
      w.A[1] = A[ 1]*x.A[0] + A[ 3]*x.A[1] + A[ 4]*x.A[2];
      w.A[2] = A[ 2]*x.A[0] + A[ 4]*x.A[1] + A[ 5]*x.A[2];
      w.A[3] = A[ 3]*x.A[0] + A[ 6]*x.A[1] + A[ 7]*x.A[2];
      w.A[4] = A[ 4]*x.A[0] + A[ 7]*x.A[1] + A[ 8]*x.A[2];
      w.A[5] = A[ 5]*x.A[0] + A[ 8]*x.A[1] + A[ 9]*x.A[2];
      w.A[6] = A[ 6]*x.A[0] + A[10]*x.A[1] + A[11]*x.A[2];
      w.A[7] = A[ 7]*x.A[0] + A[11]*x.A[1] + A[12]*x.A[2];
      w.A[8] = A[ 8]*x.A[0] + A[12]*x.A[1] + A[13]*x.A[2];
      w.A[9] = A[ 9]*x.A[0] + A[13]*x.A[1] + A[14]*x.A[2];
#endif
    }
    //--------------------------------------------------------------------------
    // ascii I/O                                                                
    //--------------------------------------------------------------------------
    friend std::ostream& OP<< (std::ostream&s, c_tsy4 x)
    {
      return write_array(s,x.A,NDAT);
    }
    //--------------------------------------------------------------------------
    friend std::istream& OP>> (std::istream&s, tsy4 &x)
    {
      return read_array(s,x.A,NDAT);
    }
  };
  //////////////////////////////////////////////////////////////////////////////
#undef CI
#undef CX
#undef FI
#undef FX
}
////////////////////////////////////////////////////////////////////////////////
#define falcON_SYM(P,NAME)				\
  register real SYM##P##_##NAME[sym##P<real>::NDAT];	\
  register sym##P<real> NAME(SYM##P##_##NAME)
#define falcON_SYM1(NAME) falcON_SYM(1,NAME)
#define falcON_SYM2(NAME) falcON_SYM(2,NAME)
#define falcON_SYM3(NAME) falcON_SYM(3,NAME)
#define falcON_SYM4(NAME) falcON_SYM(4,NAME)
////////////////////////////////////////////////////////////////////////////////
#undef OP
#undef TV
#undef c_VECT
////////////////////////////////////////////////////////////////////////////////
#endif // falcON_included_tens_h
