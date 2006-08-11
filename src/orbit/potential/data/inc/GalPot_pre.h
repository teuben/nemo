// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// GalPot_pre.h                                                                |
//                                                                             |
// C++ code written by Walter Dehnen, 1996-2002,                               |
//                                                                             |
// present address:                                                            |
// Astrophysikalisches Institut Potsdam                                        |
// An der Sternwarte 16, D-14482 Potsdam, Germany                              |
// e@mail: wdehnen@aip.de                                                      |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
// TO BE INCLUDED BY GalPot.h WHICH IN TURN IS INCLUDED BY THE ENDUSER         |
//                                                                             |
// Version 0.0    15. July      1997                                           |
// Version 0.1    24. March     1998                                           |
// Version 0.2    05. October   1998                                           |
// Version 0.3    07. June      2001                                           |
// Version 0.4    22. April     2002                                           |
//                                                                             |
//-----------------------------------------------------------------------------+
#ifndef GalPot_pre
#define GalPot_pre 1
#include <cmath>                                 // v0.4                        
#include <cstdlib>                               // v0.4                        
namespace GalPot {                               // v0.4                        
////////////////////////////////////////////////////////////////////////////////
// a reduced version of class template tupel                                    
////////////////////////////////////////////////////////////////////////////////
  template<class T, int N>
  class tupel {
  protected:
    T a[N];
    static void division_by_zero_error();
  public:
    tupel() {}
    tupel(const T);
    tupel(const T*);
    tupel(const tupel&);
    ~tupel() {}

    tupel&  operator=  (const tupel&);
    tupel&  operator+= (const tupel&);
    tupel&  operator-= (const tupel&);
    tupel&  operator=  (const T);
    tupel&  operator=  (const T*);
    tupel&  operator+= (const T);
    tupel&  operator-= (const T);
    tupel&  operator*= (const T);
    tupel&  operator/= (const T);
    tupel&  apply      (T(*)(T));

    tupel   operator-  () const;
    tupel   operator+  (const tupel&) const;
    tupel   operator-  (const tupel&) const;
    T       operator*  (const tupel&) const;
    bool    operator== (const tupel&) const;
    bool    operator!= (const tupel&) const;
    tupel   operator+  (const T) const;
    tupel   operator-  (const T) const;
    tupel   operator*  (const T) const;
    tupel   operator/  (const T) const;

    T       operator() (const int n) const { return  a[n]; }
    T&      operator[] (const int n)       { return  a[n]; }
    int     NumberofTerms() const { return N; }
    T       norm      () const;

    operator T*       ()	 { return a; }
    operator const T* () const	 { return a; }
  };
#define TI template<class T, int N> inline 

  TI tupel<T,N> operator+ (const T x, const tupel<T,N>& V) {
    return V+x;
  }

  TI tupel<T,N> operator- (const T x, const tupel<T,N>& V) {
    register tupel<T,N> P(x);
    return P-=V;
  }

  TI tupel<T,N> operator* (const T x, const tupel<T,N>& V) {
    return V*x;
  }

  TI void tupel<T,N>::division_by_zero_error() {
    std::cerr << " tupel: division by zero \n";
    std::exit(1);
  }

  TI tupel<T,N>::tupel(const T fill_value) {
    for(register int i=0; i<N; i++) a[i] = fill_value;
  }

  TI tupel<T,N>::tupel(const T *array) {
    for(register int i=0; i<N; i++) a[i] = array[i];
  }

  TI tupel<T,N>::tupel(const tupel<T,N>& V) {
    for(register int i=0; i<N; i++) a[i] = V.a[i];
  }

  TI tupel<T,N>& tupel<T,N>::operator= (const tupel<T,N>& V) {
    for(register int i=0; i<N; i++) a[i] = V.a[i];
    return *this;
  }

  TI tupel<T,N>& tupel<T,N>::operator= (const T fill_value) {
    for(register int i=0; i<N; i++) a[i] = fill_value;
    return *this;
  }

  TI tupel<T,N>& tupel<T,N>::operator= (const T* array) {
    for(register int i=0; i<N; i++) a[i] = array[i];
    return *this;
  }

  TI tupel<T,N>& tupel<T,N>::operator+= (const tupel<T,N>& V) {
    for(register int i=0; i<N; i++) a[i] += V.a[i];
    return *this;
  }

  TI tupel<T,N>& tupel<T,N>::operator-= (const tupel<T,N>& V) {
    for(register int i=0; i<N; i++) a[i] -= V.a[i];
    return *this;
  }

  TI tupel<T,N>& tupel<T,N>::operator+= (const T m) {
    for(register int i=0; i<N; i++) a[i] += m;
    return *this;
  }

  TI tupel<T,N>& tupel<T,N>::operator-= (const T m) {
    for(register int i=0; i<N; i++) a[i] -= m;
    return *this;
  }

  TI tupel<T,N>& tupel<T,N>::operator*= (const T m) {
    for(register int i=0; i<N; i++) a[i] *= m;
    return *this;
  }

  TI tupel<T,N>& tupel<T,N>::operator/= (const T m) {
    if(m==T(0.)) division_by_zero_error();
    for(register int i=0; i<N; i++) a[i] /= m;
    return *this;
  }

  TI tupel<T,N>& tupel<T,N>::apply    ( T(*f)(T) ) {
    for(register int i=0; i<N; i++) a[i] = f(a[i]);
    return *this;
  }

  TI tupel<T,N> tupel<T,N>::operator- () const {
    register tupel<T,N> P(0);
    return P-=*this;
  }

  TI tupel<T,N> tupel<T,N>::operator+ (const tupel<T,N>& V) const {
    register tupel<T,N> P(*this);
    return P+=V;
  }

  TI tupel<T,N> tupel<T,N>::operator- (const tupel<T,N>& V) const {
    register tupel<T,N> P(*this);
    return P-=V;
  }

  TI T tupel<T,N>::operator* (const tupel<T,N>& V) const {
    register T x=a[0] * V.a[0];
    for(register int i=1; i<N; i++) x += a[i] * V.a[i];
    return x;
  }

  TI bool tupel<T,N>::operator== (const tupel<T,N>& V) const {
    for(register int i=0; i<N; i++) if(a[i] != V.a[i]) return false;
    return true;
  }

  TI bool tupel<T,N>::operator!= (const tupel<T,N>& V) const {
    for(register int i=0; i<N; i++) if(a[i] != V.a[i]) return true;
    return false;
  }

  TI tupel<T,N> tupel<T,N>::operator+ (const T x) const {
    register tupel<T,N> P(*this);
    return P+=x;
  }

  TI tupel<T,N> tupel<T,N>::operator- (const T x) const {
    register tupel<T,N> P(*this);
    return P-=x;
  }

  TI tupel<T,N> tupel<T,N>::operator* (const T x) const {
    register tupel<T,N> P(*this);
    return P*=x;
  }

  TI tupel<T,N> tupel<T,N>::operator/ (const T x) const {
    register tupel<T,N> P(*this);
    return P/=x;
  }

  TI T tupel<T,N>::norm() const {
    register T x = a[0] * a[0];
    for(register int i=1; i<N; i++) x += a[i] * a[i];
    return x;
  }

  TI T norm(tupel<T,N> const&t) { return t.norm(); }

  TI std::ostream& operator<< (std::ostream& s, const tupel<T,N>& V) {
    s << V(0);
    for(register int i=1; i<N; i++) s << ' ' << V(i);
    return s;
  }

  TI std::istream& operator>> (std::istream& s, tupel<T,N>& V) { 
    T x[N];
    char c=0;
    register int i;
    c = s.get();
    if(c == '(') {
      for(i=0; i<N; i++) s >> x[i];
      s >> c;
      if(c != ')') s.clear(std::ios::badbit);
    } else {
      s.unget();
      for(i=0; i<N; i++) s >> x[i];
    }
    if(s) for(i=0; i<N; i++) V[i] = x[i];
    return s;
  }
#undef TI
  //////////////////////////////////////////////////////////////////////////////
  typedef tupel<double,3> Frequs;
  typedef tupel<double,5> DiskPar;
  typedef tupel<double,6> SphrPar;
  //////////////////////////////////////////////////////////////////////////////
  class PotResidual {
  public:
    virtual double Density  (const double, const double)                const=0;
    virtual double Residual (const double, const double, const double)  const=0;
  };
  //////////////////////////////////////////////////////////////////////////////
  class DiskAnsatz : public PotResidual {
  private:
    double               S0, Rd, zd, R0, eps;           // defining  variables  
    int                  thin, hollow, isothermal;      // auxiliary variable   
    double               Rd2, zdoRd, fac, R0oRd;        // auxiliary variables  
    double               mass_integrand(const double) const;
  public:
    void                 setup(const DiskPar&);
    DiskAnsatz           () {}
    DiskAnsatz           (const DiskPar&d) { setup(d); }
    double operator()    (const double, const double, const double, double* =0)
                                                                        const;
    double      Laplace       (const double, const double) const;
    Frequs      kapnuom       (const double) const;
    bool        is_thin       () const { return thin; }
    bool        is_hollow     () const { return hollow; }
    double      mass          (const double=0.) const;
    double      SurfaceDensity(const double) const;
    double      Density       (const double, const double) const;
    double      Residual      (const double, const double, const double) const;
    void        DescribePot   (std::ostream&) const;
    DiskPar     parameter     () const;
  };
  inline DiskPar DiskAnsatz::parameter() const
  {
    DiskPar p;
    p[0]=S0; p[1]=Rd; p[2]=(isothermal)? -zd : zd; p[3]=R0; p[4]=eps;
    return p;
  }
  //////////////////////////////////////////////////////////////////////////////
  class Disks {
  protected:
    int                   nd;
    DiskAnsatz            *D, *Dup;
    void reset            (const int, const DiskPar*);
  public:
    Disks                 (std::istream&);
    Disks                 (const Disks&); 
    Disks                 (const int, const DiskPar*);
   ~Disks                 () { delete[] D; }
    bool    all_hollow    () const;
    bool    none_hollow   () const;
    double  Mass          (const double=0.) const;
    double  SurfaceDensity(const double) const;
    double  operator()    (const double, const double) const;
    double  operator()    (const double, const double, double&, double&) const;
    double  Laplace       (const double, const double) const;
    double  Density       (const double, const double) const;
    double  Residual      (const double, const double, const double) const;
    int     NumberofDisks () const { return nd; }
    DiskPar Parameter     (const int i) const { return (D+i)->parameter(); }
  };
  inline bool Disks::none_hollow() const { 
    if(nd==0) return 1;
    for(register DiskAnsatz *p=D; p<Dup; p++) if( (p->is_hollow()) ) return 0;
    return 1;
  }
  inline bool Disks::all_hollow() const {
    if(nd==0) return 0;
    for(register DiskAnsatz *p=D; p<Dup; p++) if( !(p->is_hollow()) ) return 0;
    return 1;
  }
  inline double Disks::Mass(const double r) const { 
    if(nd==0) return 0.; 
    register double R=0.;
    for(register DiskAnsatz *p=D; p<Dup; p++) R += p->mass(r);
    return R;
  }
  inline double Disks::SurfaceDensity(const double a) const {
    if(nd==0) return 0.;
    register double R=0.;
    for(register DiskAnsatz *p=D; p<Dup; p++) R += p->SurfaceDensity(a);
    return R;
  }
  inline double Disks::operator() (const double R, const double z) const {
    if(nd==0) return 0.;
    register double r=hypot(R,z), pot=0.;
    for(register DiskAnsatz *p=D; p<Dup; p++) pot+= (*p)(R,z,r);
    return pot;
  }
  inline double Disks::operator() (const double R, const double z,
				   double& dR, double& dz) const {
    if(nd==0) {
      dR = 0.;
      dz = 0.;
      return 0.;
    }
    register double d[2], r=hypot(R,z), pot=0.;
    register DiskAnsatz *p=D;
    for(dR=dz=0.; p<Dup; p++)
      { pot += (*p)(R,z,r,d); dR+=d[0]; dz+=d[1]; }
    return pot;
  }
  inline double Disks::Laplace(const double a, const double b) const { 
    if(nd==0) return 0.;
    register double L=0.;
    for(register DiskAnsatz *p=D; p<Dup; p++) L += p->Laplace(a,b);
    return L;
  }
  inline double Disks::Density(const double a, const double b) const {
    if(nd==0) return 0.;
    register double R=0.;
    for(register DiskAnsatz *p=D; p<Dup; p++) R += p->Density(a,b);
    return R;
  }
  inline double Disks::Residual(const double a, const double b, const double c)
    const {
    if(nd==0) return 0.;
    register double R=0.;
    for(register DiskAnsatz *p=D; p<Dup; p++) R += p->Residual(a,b,c);
    return R;
  }
  //////////////////////////////////////////////////////////////////////////////
  class SpheroidDensity : public PotResidual {
  private:
    double rh0, q, gam, bet, r0, rcut;              // defining variables       
    double beg, qi, r0i, rci;                       // auxiliary variables      
    double mass_integrand(const double) const;      // auxiliary function       
  public:
    void    setup(const SphrPar&);
    SpheroidDensity        () {}
    SpheroidDensity        (const SphrPar &s)       { setup(s); }
    bool    cut_off        () const { return (rci>0.);  }
    double  scale_density  () const { return rh0;  }
    double  inner_power    () const { return gam;  }
    double  outer_power    () const { return (rci)? 1.e3 : bet; }
    double  mass           (const double) const;
    double  Density        (const double, const double) const;
    double  Residual       (const double, const double, const double) const;
    SphrPar parameter      () const;
  };
  inline double SpheroidDensity::Residual(const double r, const double st,
					  const double ct) const {
    return fPiG * Density(r*st,r*ct);
  }
  inline SphrPar SpheroidDensity::parameter() const {
    SphrPar p;
    p[0]=rh0;
    p[1]=q;
    p[2]=gam;
    p[3]=bet;
    p[4]=r0;
    p[5]=rcut;
    return p;
  }
  //////////////////////////////////////////////////////////////////////////////
  class Spheroids {
    int                      ns;
    SpheroidDensity          *S, *Sup;
  protected:
    void      reset          (const int, const SphrPar*);
  public:
    Spheroids                (std::istream&);
    Spheroids                (const int, const SphrPar*);
    Spheroids                (const Spheroids&);
    Spheroids& operator=     (const Spheroids&);
    ~Spheroids                () { delete[] S; }
    bool    massive          () const;
    double  beta             () const;
    double  gamma            () const;
    double  Mass             (const double) const;
    double  Density          (const double, const double) const;
    double  Residual         (const double, const double, const double)const;
    int     NumberofSpheroids() const { return ns; }
    SphrPar Parameter        (const int i) const { return (S+i)->parameter(); }
  };
  inline bool Spheroids::massive() const {
    if(ns==0) return 0;
    for(register SpheroidDensity *p=S; p<Sup;p++)
      if(p->scale_density()) return true;
    return false;
  }
  inline double Spheroids::Mass(const double a) const {
    if(ns==0.) return 0.;
    register double R=0.;
    for(register SpheroidDensity *p=S; p<Sup; p++) R += p->mass(a);
    return R;
  }
  inline double Spheroids::Density(const double a, const double b) const {
    if(ns==0.) return 0.;
    register double R=0.;
    for(register SpheroidDensity *p=S; p<Sup; p++) R += p->Density(a,b);
    return R;
  }
  inline double Spheroids::Residual(const double a, const double b,
				    const double c) const {
    if(ns==0) return 0.;
    register double R=0.;
    for(register SpheroidDensity *p=S; p<Sup; p++) R += p->Residual(a,b,c);
    return R;
  }
  //////////////////////////////////////////////////////////////////////////////
  class Multipole {
  private:
    int         LR,K[2];
    double      Rmin, Rmax, gamma, beta, Phi0;
    double      lRmin, lRmax, g2;
    double      lzmin, lzmax, tg3, g3h;
    double      *logr, *lLc, *d2R, *d2L; 
    double      *X[2], **Y[3], **Z[4];
    void        AllocArrays();
    void        setup(const double, const double,       // r_min, r_max         
                      const double, const double,       // gamma, beta          
                      PotResidual const*);              // providing rho(x)     
  public:
    Multipole (const int,                               // points on log grid   
               const double, const double,              // r_min, r_max         
               const double, const double,              // gamma, beta          
               PotResidual const*,                      // providing rho(x)     
               const int =1);                           // routines LfromRc...  
    void reset(const double, const double,              // r_min, r_max         
               const double, const double,              // gamma, beta          
               PotResidual const*,                      // providing rho(x)     
               const int =1);                           // routines LfromRc...  
   ~Multipole();
    double      operator()(const double,const double,const double, double* =0) 
                                                        const;
    double      vcsquare  (const double)		const;
    double      vcsquare  (const double, double&)	const;
    double      Laplace   (const double, const double)  const;
    Frequs      kapnuom   (const double)		const;
  };
  inline double Multipole::vcsquare(const double R) const
  { 
    double d[2];
    operator() (R,0.,1.,d);
    return R*d[0];
  }

}                                                       // namespace GalPot     
#endif

// end of file GalPot_pre.h ////////////////////////////////////////////////////
