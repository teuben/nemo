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
// Version 0.5    04. February  2005                                           |
//                                                                             |
//-----------------------------------------------------------------------------+
#ifndef GalPot_pre
#define GalPot_pre 1
#include <cmath>                                 // v0.4                        
#include <cstdlib>                               // v0.4                        
#define WD_TUPEL_FUNCOP                          // v0.5                        
#include <tupel.h>                               // v0.5                        
namespace GalPot {                               // v0.4                        
  using WDutils::tupel;                          // v0.5                        
  typedef tupel<3,double> Frequs;
  typedef tupel<5,double> DiskPar;
  typedef tupel<6,double> SphrPar;
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // abstract base class PotResidual                                          //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class PotResidual {
  public:
    // density at given (R,z)
    virtual double Density  (const double, const double)                const=0;
    // residual density (input for multipole expansion) at given (r,sin/cos(th))
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
    // return potential (and its gradient) of part not in residual
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
