// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// tupl.h                                                                      |
//                                                                             |
// C++ code                                                                    |
//                                                                             |
// Copyright Walter Dehnen, 1996-2001                                          |
// e-mail:   wdehnen@aip.de                                                    |
// address:  Astrophysikalisches Institut Potsdam,                             |
//           An der Sternwarte 16, D-14482 Potsdam, Germany                    |
//                                                                             |
//-----------------------------------------------------------------------------+

#ifndef included_tupl_h
#define included_tupl_h

#ifndef included_iostream
#  include <iostream>
#  define included_iostream
#endif

#ifndef included_inln_h
#  include <public/inln.h>
#endif

////////////////////////////////////////////////////////////////////////////////
namespace nbdy {
  //////////////////////////////////////////////////////////////////////////////
  // template class nbdy::tupel                                                 
  //////////////////////////////////////////////////////////////////////////////
#define c_REAL  const REAL
#define c_tupel const tupel
#define op      operator
#define RR      register REAL
#define FI(X)   for(register int i=X; i!=N; ++i)
  template<class REAL, int N>
  class tupel {
  private:
    REAL a[N];
  public:
    tupel          ()                   {}
   ~tupel          ()                   {}
    tupel          (c_REAL  x)          { FI(0) a[i] = x; }
    tupel          (c_REAL *b)          { FI(0) a[i] = b[i]; }
    tupel          (c_tupel&b)          { FI(0) a[i] = b.a[i]; }
    tupel& op=     (c_REAL  x)          { FI(0) a[i] = x; return *this; }
    tupel& op+=    (c_REAL  x)          { FI(0) a[i]+= x; return *this; }
    tupel& op-=    (c_REAL  x)          { FI(0) a[i]-= x; return *this; }
    tupel& op*=    (c_REAL  x)          { FI(0) a[i]*= x; return *this; }
    tupel& op/=    (c_REAL  x)          { FI(0) a[i]/= x; return *this; }
    tupel  op+     (c_REAL  x) const    { tupel y(*this); y+=x; return y; }
    tupel  op-     (c_REAL  x) const    { tupel y(*this); y-=x; return y; }
    tupel  op*     (c_REAL  x) const    { tupel y(*this); y*=x; return y; }
    tupel  op/     (c_REAL  x) const    { tupel y(*this); y/=x; return y; }
    tupel& op=     (c_REAL *b)          { FI(0) a[i] = b[i]; return *this; }
    tupel& op=     (c_tupel&b)          { FI(0) a[i] = b[i]; return *this; }
    tupel& op+=    (c_tupel&b)          { FI(0) a[i]+= b[i]; return *this; }
    tupel& op-=    (c_tupel&b)          { FI(0) a[i]-= b[i]; return *this; }
    tupel& op*=    (c_tupel&b)          { FI(0) a[i]*= b[i]; return *this; }
    tupel& op/=    (c_tupel&b)          { FI(0) a[i]/= b[i]; return *this; }
    tupel  op+     (c_tupel&b) const    { tupel y(*this); y+=b; return y; }
    tupel  op-     (c_tupel&b) const    { tupel y(*this); y-=b; return y; }
    tupel  op-	   ()		        { tupel y(REAL(0)); return y-= *this; }
    REAL&  op[]    (const int i)        { return a[i]; }
    c_REAL&op[]    (const int i) const  { return a[i]; }
    op REAL*       ()                   { return a; }
    op c_REAL*     () const             { return a; }
    bool   op==    (c_tupel &b) const   { FI(0) if(a[i] != b[i]) return false; 
                                          return true; }
    bool   op!=    (c_tupel &b) const   { FI(0) if(a[i] != b[i]) return true;
                                          return false; }
    bool   op<     (c_tupel &b) const   { FI(0) if(a[i] >= b[i]) return false; 
                                          return true; }
    bool   op<=    (c_tupel &b) const   { FI(0) if(a[i] >  b[i]) return false; 
                                          return true; }
    bool   op>     (c_tupel &b) const   { FI(0) if(a[i] <= b[i]) return false; 
                                          return true; }
    bool   op>=    (c_tupel &b) const   { FI(0) if(a[i] <  b[i]) return false; 
                                          return true; }
    void   apply   (REAL(*Y)(c_REAL&))  { FI(0) a[i] = Y(a[i]); }
    void   connect (c_tupel &b, REAL(*Y)(c_REAL&,c_REAL&))
                                        { FI(0) a[i] = Y(a[i],b[i]); }
    REAL   norm    () const             { RR x=0; FI(0) x+=a[i]*a[i]; return x; }
    REAL   dist_sq (c_tupel &b) const   { RR x=0; FI(0) x+=square(a[i]-b[i]);
                                          return x; }
    REAL   op*     (c_tupel &b) const   { RR x=0;   FI(0) x+=a[i]*b[i];
                                          return x; }
    REAL   maxnorm () const             { RR b,x=0; FI(0) if(x<(b=abs(a[i]))) x=b;
                                          return x; }
  };
#undef c_tupel
  //////////////////////////////////////////////////////////////////////////////
  // related functions                                                          
  //////////////////////////////////////////////////////////////////////////////
#define c_tupel const tupel<REAL,N>
#define TN      template<class REAL, int N> inline
#define TUPEL   tupel<REAL,N>
  TN TUPEL op+     (c_REAL x, c_tupel &y)  { return y+x; }
  TN TUPEL op-     (c_REAL x, c_tupel &y)  { return -(y-x); }
  TN TUPEL op*     (c_REAL x, c_tupel &y)  { return y*x; }
  TN REAL norm     (c_tupel x)             { return x.norm(); }
  TN REAL dist_sq  (c_tupel x, c_tupel y)  { return x.dist_sq(y); }
  TN REAL maxnorm  (c_tupel x)             { return x.maxnorm(); }
  //----------------------------------------------------------------------------
  TN std::ostream& op<< (std::ostream& s, c_tupel& a) {
    return write_array(s,a,N);
  }
  //----------------------------------------------------------------------------
  TN std::istream& op>> (std::istream& s, TUPEL& a) {
    return read_array(s,a,N);
  }
  //////////////////////////////////////////////////////////////////////////////
#undef c_REAL
#undef c_tupel
#undef op
#undef RR
#undef FI
#undef TN
#undef TUPEL
}
////////////////////////////////////////////////////////////////////////////////
#endif	// included_tupel_h
