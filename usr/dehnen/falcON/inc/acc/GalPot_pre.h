// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// GalPot_pre.h                                                                |
//                                                                             |
// Copyright (C) 1996-2005 Walter Dehnen                                       |
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
// TO BE INCLUDED BY GalPot.h WHICH IN TURN IS INCLUDED BY THE ENDUSER         |
//                                                                             |
// Version 0.0    15. July      1997                                           |
// Version 0.1    24. March     1998                                           |
// Version 0.2    05. October   1998                                           |
// Version 0.3    07. June      2001                                           |
// Version 0.4    22. April     2002                                           |
// Version 0.5    04. February  2005                                           |
// Version 0.6    06. November  2007                                           |
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
    virtual double Density  (double, double) const=0;
    // residual density (input for multipole expansion) at given (r,sin/cos(th))
    virtual double Residual (double, double, double)  const=0;
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class DiskAnsatz                                                         //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class DiskAnsatz : public PotResidual {
  private:
    double               S0, Rd, zd, R0, eps;           // defining  variables  
    int                  thin, hollow, isothermal;      // auxiliary variable   
    double               Rd2, zdoRd, fac, R0oRd;        // auxiliary variables  
    double               mass_integrand(double) const;
  public:
    void                 setup(const DiskPar&);
    DiskAnsatz           () {}
    DiskAnsatz           (const DiskPar&d) { setup(d); }
    // return potential (and its gradient) of part not in residual
    double operator()    (double, double, double, double* =0)
                                                                        const;
    double      Laplace       (double, double) const;
    Frequs      kapnuom       (double) const;
    bool        is_thin       () const { return thin; }
    bool        is_hollow     () const { return hollow; }
    double      mass          (double=0.) const;
    double      SurfaceDensity(double) const;
    double      Density       (double, double) const;
    double      Residual      (double, double, double) const;
    void        DescribePot   (std::ostream&) const;
    DiskPar     parameter     () const;
  };
  inline DiskPar DiskAnsatz::parameter() const {
    return DiskPar(S0, Rd, isothermal? -zd : zd, R0, eps);
  }
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class Disks                                                              //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class Disks {
  protected:
    int                   nd;
    DiskAnsatz            *D, *Dup;
    void reset            (int, const DiskPar*);
  public:
    Disks                 (std::istream&);
    Disks                 (const Disks&); 
    Disks                 (int, const DiskPar*);
    ~Disks                () { if(D) delete[] D; }
    bool    all_hollow    () const;
    bool    none_hollow   () const;
    double  Mass          (double=0.) const;
    double  SurfaceDensity(double) const;
    double  operator()    (double, double) const;
    double  operator()    (double, double, double&, double&) const;
    double  Laplace       (double, double) const;
    double  Density       (double, double) const;
    double  Residual      (double, double, double) const;
    int     NumberofDisks () const { return nd; }
    DiskPar Parameter     (int i) const { return (D+i)->parameter(); }
  };
  inline bool Disks::none_hollow() const { 
    for(DiskAnsatz *p=D; p!=Dup; ++p)
      if( p->is_hollow() ) return false;
    return true;
  }
  inline bool Disks::all_hollow() const {
    for(DiskAnsatz *p=D; p!=Dup; ++p)
      if( !(p->is_hollow()) ) return false;
    return true;
  }
  inline double Disks::Mass(double r) const { 
    double M=0.;
    for(DiskAnsatz *p=D; p!=Dup; ++p)
      M += p->mass(r);
    return M;
  }
  inline double Disks::SurfaceDensity(double a) const {
    double R=0.;
    for(DiskAnsatz *p=D; p!=Dup; ++p)
      R += p->SurfaceDensity(a);
    return R;
  }
  inline double Disks::operator() (double R, double z) const {
    double r=hypot(R,z), pot=0.;
    for(DiskAnsatz *p=D; p!=Dup; ++p)
      pot += (*p)(R,z,r);
    return pot;
  }
  inline double Disks::operator() (double R, double z,
				   double&dR, double&dz) const {
    dR = 0.;
    dz = 0.;
    double r=hypot(R,z), d[2], pot=0.;
    for(DiskAnsatz *p=D; p!=Dup; ++p) {
      pot += (*p)(R,z,r,d);
      dR  += d[0];
      dz  += d[1];
    }
    return pot;
  }
  inline double Disks::Laplace(double a, double b) const { 
    double L=0.;
    for(DiskAnsatz *p=D; p!=Dup; ++p)
      L += p->Laplace(a,b);
    return L;
  }
  inline double Disks::Density(double a, double b) const {
    double R=0.;
    for(DiskAnsatz *p=D; p!=Dup; ++p)
      R += p->Density(a,b);
    return R;
  }
  inline double Disks::Residual(double a, double b, double c) const {
    double R=0.;
    for(DiskAnsatz *p=D; p!=Dup; ++p)
      R += p->Residual(a,b,c);
    return R;
  }
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class SpheroidDensity                                                    //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class SpheroidDensity : public PotResidual {
  private:
    double rh0, q, gam, bet, r0, rcut;              // defining variables       
    double beg, qi, r0i, rci;                       // auxiliary variables      
    double mass_integrand(double) const;            // auxiliary function       
  public:
    void    setup(const SphrPar&);
    SpheroidDensity        () {}
    SpheroidDensity        (const SphrPar &s)       { setup(s); }
    bool    cut_off        () const { return (rci>0.);  }
    double  scale_density  () const { return rh0;  }
    double  inner_power    () const { return gam;  }
    double  outer_power    () const { return (rci)? 1.e3 : bet; }
    double  mass           (double) const;
    double  Density        (double, double) const;
    double  Residual       (double, double, double) const;
    SphrPar parameter      () const;
  };
  inline double SpheroidDensity::Residual(double r, double st, double ct) const
  {
    return fPiG * Density(r*st,r*ct);
  }
  inline SphrPar SpheroidDensity::parameter() const {
    return SphrPar(rh0,q,gam,bet,r0,rcut);
  }
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class Spheroids                                                          //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class Spheroids {
    int                      ns;
    SpheroidDensity          *S, *Sup;
  protected:
    void      reset          (int, const SphrPar*);
  public:
    Spheroids                (std::istream&);
    Spheroids                (int, const SphrPar*);
    Spheroids                (const Spheroids&);
    Spheroids& operator=     (const Spheroids&);
    ~Spheroids               () { if(S) delete[] S; }
    bool    massive          () const;
    double  beta             () const;
    double  gamma            () const;
    double  Mass             (double) const;
    double  Density          (double, double) const;
    double  Residual         (double, double, double)const;
    int     NumberofSpheroids() const { return ns; }
    SphrPar Parameter        (int i) const { return (S+i)->parameter(); }
  };
  inline bool Spheroids::massive() const {
    for(register SpheroidDensity *p=S; p!=Sup; ++p)
      if(p->scale_density()) return true;
    return false;
  }
  inline double Spheroids::Mass(double a) const {
    double R=0.;
    for(SpheroidDensity *p=S; p!=Sup; ++p)
      R += p->mass(a);
    return R;
  }
  inline double Spheroids::Density(double a, double b) const {
    double R=0.;
    for(SpheroidDensity *p=S; p!=Sup; ++p)
      R += p->Density(a,b);
    return R;
  }
  inline double Spheroids::Residual(double a, double b, double c) const {
    double R=0.;
    for(SpheroidDensity *p=S; p!=Sup; ++p)
      R += p->Residual(a,b,c);
    return R;
  }
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class Multipole                                                          //
  //                                                                          //
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
    void        setup(double, double,                   // r_min, r_max         
                      double, double,                   // gamma, beta          
                      PotResidual*);                    // providing rho(x)     
  public:
    Multipole (int,                                     // points on log grid   
               double, double,                          // r_min, r_max         
               double, double,                          // gamma, beta          
               PotResidual*,                            // providing rho(x)     
               int =1);                                 // routines LfromRc...  
    void reset(double, double,                          // r_min, r_max         
               double, double,                          // gamma, beta          
               PotResidual*,                            // providing rho(x)     
               int =1);                                 // routines LfromRc...  
   ~Multipole();
    double      operator()(double, double, double, double* =0) const;
    double      vcsquare  (double) const;
    double      vcsquare  (double, double&) const;
    double      Laplace   (double, double) const;
    Frequs      kapnuom   (double) const;
  };
  inline double Multipole::vcsquare(double R) const
  { 
    double d[2];
    operator() (R,0.,1.,d);
    return R*d[0];
  }

} // namespace GalPot     
#endif

// end of file GalPot_pre.h ////////////////////////////////////////////////////
