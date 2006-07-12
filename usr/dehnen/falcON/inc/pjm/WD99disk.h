// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// WD99disk.h                                                                  |
//                                                                             |
// C++ code                                                                    |
//                                                                             |
// Copyright Walter Dehnen, 2002-2005                                          |
//           Paul McMillan, 2005-2006                                          |
// e-mail:   walter.dehnen@astro.le.ac.uk                                      |
//           paul.mcmillan@astro.le.ac.uk                                      |
// address:  Department of Physics and Astronomy, University of Leicester      |
//           University Road, Leicester LE1 7RH, United Kingdom                |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
// class for a disk's particle distribution based on Walter's paper in 99      |
//                                                                             |
//-----------------------------------------------------------------------------+

#ifndef falcON_included_WD99disk_h
#define falcON_included_WD99disk_h
#ifndef falcON_included_rand_h
#include <public/random.h>
#endif
#ifndef falcON_included_externacc_h
#  include<externacc.h>
#endif

namespace falcON {
typedef tupel<3,double> vectd;
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class nbdy::WD99disk                                                     //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class WD99disk 
  {
  private:
    int          n,n1;                   // Table sizes, other
    double       *lr,*pot,*dpdr,sig0;
    const int    No;                            // parameter
    const double Rd,iRd,Dens0,Rsig,Qmin,Z0,Eps,Mt,rmin,rmax;
    const        acceleration *Pe;              // External combi-potential
    const        ExpDisk Disk;                  // disc for particles

    //--------------------------------------------------------------------------
    struct PlanarOrbit {
      bool   range;
      int    i,Npoints;
      double re,xi,rs,s0,ire,POT,CENACC,Q,sdens;
      double Ec,Lc,Omc,kap,gam,sigre,Lorb,e,g2,omr,Tr,dT;
      double *ttable,*Rtable,*vRtable;
      vectd  *W;

      PlanarOrbit(double R, 
		  double Xi, 
		  double Rs, 
		  double S0, 
		  double Qmin, 
		  double Sdens,
		  double Sigcorr);

      void sample( // samples a point on the orbit
		  Random const&,           // I: random
		  bool   const&,           // I: quasi?
		  double      &,           // O: radius
		  double      &,           // O: radial velocity
		  double      &,           // O: azimuth
		  double      &) const;    // 0: azimuthal velocity

      vectd Dt(vectd const&);
      vectd DR(vectd const&);
      void step_back_to_R_equal(const double , vectd& );
      double CashKarp(vectd& , const double , const double );
      void CashKarpStep(vectd& , const double);

      ~PlanarOrbit();
    };
    //--------------------------------------------------------------------------
  public:
    WD99disk ( int ,                        // Number of orbits
	       double ,                     // Disk scale radius
	       double ,                     // Disk scale radius
	       double ,                     // Vdisp scale radius
	       double ,                     // Sigma 0
	       double ,                     // scale height
	       double ,                     // smoothing length
	       const acceleration * = 0);   // External potential

    void sample( bodies const&,            // I/O: bodies to sample
		 int    const&,            // I: no. iterations
		 bool   const&,            // I: quasi random?          
		 Random const&) const;     // I: pseudo & quasi RNG
    
    
    // when creating a cold disk, the programme uses: 
    void coldsample( bodies const&,            // I/O: bodies to sample
		     bool   const&,            // I: quasi random?          
		     Random const&) const;     // I: pseudo & quasi RNG

    // to iterate the density distribution closer to 
    // the target, the programme uses:
    void iterate( int const& ,
		  int const&,
		  bool   const&,
		  Random const&,
		  double * =0,
		  double * =0,
		  double * =0,
		  double * =0) const;

    ~WD99disk();                              // destruction               
  };
}

////////////////////////////////////////////////////////////////////////////////
#endif                                             // falcON_included_WD99disk_h    

