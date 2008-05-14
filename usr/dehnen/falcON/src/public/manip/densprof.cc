// -*- C++ -*-                                                                  
////////////////////////////////////////////////////////////////////////////////
///                                                                             
/// \file   src/public/manip/densprof.cc                                        
///                                                                             
/// \author Walter Dehnen                                                       
/// \date   2006,2008                                                           
///                                                                             
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// Copyright (C) 2006,2008 Walter Dehnen                                        
//                                                                              
// This program is free software; you can redistribute it and/or modify         
// it under the terms of the GNU General Public License as published by         
// the Free Software Foundation; either version 2 of the License, or (at        
// your option) any later version.                                              
//                                                                              
// This program is distributed in the hope that it will be useful, but          
// WITHOUT ANY WARRANTY; without even the implied warranty of                   
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU            
// General Public License for more details.                                     
//                                                                              
// You should have received a copy of the GNU General Public License            
// along with this program; if not, write to the Free Software                  
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.                    
//                                                                              
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// history:                                                                     
//                                                                              
// v 0.0    02/05/2006  WD created                                              
// v 0.1    07/07/2006  WD using bodies in_subset()                             
// v 0.2    27/07/2006  WD made public                                          
// v 0.3    13/05/2008  WD debugged error with minor & major axes               
////////////////////////////////////////////////////////////////////////////////
#include <public/defman.h>
#include <public/io.h>
#include <utils/numerics.h>
#include <ctime>
#include <cstring>

namespace falcON { namespace Manipulate {
  using std::setprecision;
  using std::setw;
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // auxiliaries                                                              //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  inline real neg_density(body const&B) {
    return -rho(B);
  }
  //////////////////////////////////////////////////////////////////////////////
  const char*Line ="----------------------------------------";
  //////////////////////////////////////////////////////////////////////////////
  class PrintSmall {
    int N,P;
  public:
    PrintSmall(int n) : N(n), P(1) {
      for(int i=0; i!=n; ++i) P *= 10;
    }
    template<typename X>
    std::ostream&print_pos(std::ostream&o, X const&x) const {
      int I = int(P*x+0.5);
      if(I >= P) return o << "1." << std::setw(N-1) << std::setfill('0') << 0
			  << std::setfill(' ');
      else       return o << '.' << std::setw(N) << std::setfill('0') << I
			  << std::setfill(' ');
    }
    template<typename X>
    std::ostream&print(std::ostream&o, X const&x) const {
      if(x < 0) { o << '-'; return print_pos(o,-x); }
      else      { o << ' '; return print_pos(o, x); }
    }
    template<typename X>
    std::ostream&print_dir(std::ostream&o, tupel<3,X> const&x) const {
      return print(print(print(o,x[0])<<' ',x[1])<<' ',x[2]);
    }
    template<typename X>
    std::ostream&print_dir(std::ostream&o, const X x[3]) const {
      return print(print(print(o,x[0])<<' ',x[1])<<' ',x[2]);
    }
    const char*line_pos() const { return Line+(39-N); }
    const char*line    () const { return Line+(38-N); }
    const char*line_dir() const { return Line+(38-3*(N+2)); }
  };
  //----------------------------------------------------------------------------
  template<typename X> inline
  void add_outer_product(X p[3][3], tupel<3,X> const&x, X m) {
    p[0][0] += m * x[0] * x[0];
    p[0][1] += m * x[0] * x[1];
    p[0][2] += m * x[0] * x[2];
    p[1][1] += m * x[1] * x[1];
    p[1][2] += m * x[1] * x[2];
    p[2][2] += m * x[2] * x[2];
  }
  //----------------------------------------------------------------------------
  template<typename X> inline
  void symmetrize(X p[3][3]) {
    p[1][0] = p[0][1];
    p[2][0] = p[0][2];
    p[2][1] = p[1][2];
  }
  //////////////////////////////////////////////////////////////////////////////
  const int    W_default = 1000;
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class densprof                                                             
  //                                                                            
  /// manipulator: measures profiles of density bins                            
  ///                                                                           
  /// This manipulator uses a pre-computed mass-density (see manipulator        
  /// Manipulate::density) to compute the centre and profile of shells with     
  /// a small range in rho. Only bodies in_subset() are used.                   
  ///                                                                           
  /// Meaning of the parameters:\n                                              
  /// par[0]: # bodies per density shell (window size, default: 1000)\n         
  /// par[1]: delta time between analyses (default: 0)\n                        
  /// file  : format string for output table files\n                            
  ///                                                                           
  /// Usage of pointers: none\n                                                 
  /// Usage of flags:    uses in_subset()\n                                     
  ///                                                                           
  // ///////////////////////////////////////////////////////////////////////////
  class densprof : public manipulator {
  private:
    const int        W;         ///< window size
    double           STEP;      ///< delta time between manipulations
    mutable double   TMAN;      ///< time for next manipulation
    mutable bool     FST;       ///< first call to manipulate() ?
    mutable int      I;         ///< index of manipulation
    mutable output   OUT;       ///< file for output table
    char*  const     FILE;      ///< format for file name
    const PrintSmall PS;        ///< for printing out number in [0,1]
    //--------------------------------------------------------------------------
    void print_line() const;
    //--------------------------------------------------------------------------
  public:
    const char* name    () const { return "densprof"; }
    const char* describe() const {
      return
	"given density, compute centre and radial profile for "
	"bodies in_subset()";
    }
    //--------------------------------------------------------------------------
    fieldset need   () const { return fieldset::basic | fieldset::r; }
    fieldset provide() const { return fieldset::o; }
    fieldset change () const { return fieldset::o; }
    //--------------------------------------------------------------------------
    bool manipulate(const snapshot*) const;
    //--------------------------------------------------------------------------
    densprof(const double*, int npar, const char*) falcON_THROWING;
    //--------------------------------------------------------------------------
    ~densprof() { if(FILE) falcON_DEL_A(FILE); }
  };
  //////////////////////////////////////////////////////////////////////////////
  densprof::densprof(const double*pars,
		     int          npar,
		     const char  *file) falcON_THROWING
  : W    ( npar>0?     int(pars[0])    : W_default ),
    STEP ( npar>1?         pars[1]     : 0. ),
    I    ( 0 ),
    FST  ( true ),
    FILE ( (file && file[0])? falcON_NEW(char,strlen(file)+1) : 0 ),
    PS   ( 3 )
  {
    if(debug(2) || file==0 || file[0]==0 || npar>2)
      std::cerr<<
	" Manipulator \""<<name()<<"\":\n"
	" uses estimated density to compute centre and radial profiles\n"
	" par[0]: # bodies per density shell (window size, def: "
	       <<W_default<<")\n"
	" par[1]: delta time between analyses (def: 0)\n"
	" file  : format string for output table files\n";
    if(FILE) strcpy(FILE,file);
    if(file==0 || file[0]==0)
      falcON_THROW("Manipulator \"densprof\": no file given");
    if(W<=2) falcON_THROW("Manipulator \"%s\": W = %d <2\n",name(),W);
    if(npar > 2)
      warning("Manipulator \"%s\": skipping parameters beyond 2\n",name());
  }
  //////////////////////////////////////////////////////////////////////////////
  inline void densprof::print_line() const
  {
    OUT <<
      "#--------------------------------------------------------------"
      "---------------------------------------------------------------"
      "-------------"
      "-------------"
	<< PS.line_dir() << PS.line_dir() << PS.line_dir() <<'\n';
  }
  //////////////////////////////////////////////////////////////////////////////
  bool densprof::manipulate(const snapshot*SHOT) const
  {
    // 0. preliminaries
    // 0.1 first call ever:
    if(FST) {
      TMAN = SHOT->initial_time();
      FST  = false;
    }
    // 0.2 is it time for a manipulation?
    if(SHOT->time() < TMAN) return false;
    TMAN += STEP;
    // 0.3 are data sufficient?
    if(!SHOT->have_all(need()))
      error("densprof::manipulate(): need %s, but got %s\n",
	    word(need()), word(SHOT->all_data()));
    // 1. sort bodies in descending density
    Array<bodies::index> T;
    SHOT->sorted(T,&neg_density);
    const unsigned Nb = T.size();
    if(Nb < W)
      falcON_THROW("densprof::manipulate(): "
		   "fewer (%d) bodies than window (%d)\n",Nb,W);
    // 2. open output file and write header
    if(OUT.reopen(FILE,I++,1)) {
      print_line();
      OUT << "#\n"
	  << "# output from Manipulator \""<<name()<<"\"\n#\n";
      if(RunInfo::cmd_known ()) OUT<<"# command: \""<<RunInfo::cmd ()<<"\"\n";
      OUT  <<"# run at "<<RunInfo::time()<<'\n';
      if(RunInfo::user_known())
	OUT<<"#     by \""<<RunInfo::user()<<"\"\n";
      if(RunInfo::host_known())
	OUT<<"#     on \""<<RunInfo::host()<<"\"\n";
      if(RunInfo::pid_known())
	OUT<<"#     pid "<<RunInfo::pid()<<'\n';
      OUT <<"#\n";
    }
    OUT <<"# time = "<<SHOT->time()<<": "<<Nb
	<<" bodies (of "<<SHOT->N_bodies()
        <<")\n#\n"
	<<"#              xcen                          vcen           "
	<<"   radius "
	<<"      rho "
	<<"  <v_rad> "
	<<"  <v_rot> "
	<<"sigma_rad "
	<<"sigma_mer "
	<<"sigma_rot "
	<<"      c/a "
	<<"      b/a "
	<<"       major axis        minor axis     rotation axis\n";
    print_line();
    Array<real> Rq(W+W);
    Array<int>  Ir(W+W);
    // 3. loop windows of W
    for(unsigned ib=0,kb=W+W>Nb? Nb:W; ib!=Nb;
	ib=kb,kb=kb+W+W > Nb? Nb : kb+W) {
      const unsigned N = kb-ib;
      // 1st body loop: measure total mass, mean density, and centre
      double M(0.);
      vect_d X0(0.), V0(0.);
      double Rho(0.);
      for(unsigned i=ib; i!=kb; ++i) {
	const real mi = SHOT->mass(T[i]);
	M  += mi;
	X0 += mi * SHOT->pos(T[i]);
	V0 += mi * SHOT->vel(T[i]);
	Rho+= mi * SHOT->rho(T[i]);
      }
      double iM = 1./M;
      X0  *= iM;
      V0  *= iM;
      Rho *= iM;
      // 2nd body loop: measure moment of inertia, rotation and median radius
      vect_d Mvp(0.);
      double Mxx[3][3] = {{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}};
      for(unsigned i=ib,j=0; i!=kb; ++i,++j) {
	double mi = SHOT->mass(T[i]);
	vect_d ri = SHOT->pos(T[i]) - X0;
	vect_d vi = SHOT->vel(T[i]) - V0;
	vect_d er = normalized(ri);
	Mvp      += mi * (er^vi);
	Rq[j]     = norm(ri);
	add_outer_product(Mxx,ri,mi);
      }
      symmetrize(Mxx);
      double IV[3][3], ID[3];
      int    IR;
      EigenSymJacobiSorted<3,double>(Mxx,IV,ID,IR);
      double vp = abs(Mvp)*iM;
      HeapIndex(Rq.array(),N,Ir.array());
      double Mh = 0.5*M, Mc(0.), mR;
      for(unsigned j=0; j!=N; ++j) {
	double mi = SHOT->mass(T[ib+Ir[j]]);
	Mc += mi;
	if(Mc >= Mh) {
	  mR = sqrt((Rq[Ir[j-1]]*(Mc-Mh) + Rq[Ir[j]]*(Mh-Mc+mi))/mi);
	  break;
	}
      }
      // 3rd body loop: measure mean velocities and dispersion
      vect_d erot = norm(Mvp)>0.? normalized(Mvp) : vect_d(0.,0.,1.);
      double Mvr(0.), Mvrq(0.), Mvpq(0.), Mvtq(0.);
      for(unsigned i=ib; i!=kb; ++i) {
	double mi = SHOT->mass(T[i]);
	vect_d ri = SHOT->pos(T[i]) - X0;
	vect_d vi = SHOT->vel(T[i]) - V0;
	vect_d er = normalized(ri);
	vect_d ep = normalized(er^erot);
	vect_d et = normalized(er^ep);
	double ui = er*vi;
	Mvr      += mi * ui;
	Mvrq     += mi * ui*ui;
	Mvpq     += mi * square(vi*ep);
	Mvtq     += mi * square(vi*et);
      }
      double vr = Mvr*iM;
      double sr = sqrt(M*Mvrq-Mvr*Mvr)*iM;
      double st = sqrt(Mvtq*iM);
      double sp = sqrt(M*Mvpq-Mvp*Mvp)*iM;
      // finally print out the data
      OUT << print(X0 ,9,3)               <<' '  // centre position
	  << print(V0 ,9,3)               <<' '  // centre velocity
	  << print(mR ,9,3)               <<' '  // median radius
	  << print(Rho,9,3)               <<' '  // density
	  << print(vr ,9,3)               <<' '  // <v_r>
	  << print(vp ,9,3)               <<' '  // <v_rot>
	  << print(sr ,9,3)               <<' '  // sigma_r
	  << print(st ,9,3)               <<' '  // sigma_mer
	  << print(sp ,9,3)               <<' '  // sigma_rot
	  << print(sqrt(ID[2]/ID[0]),9,3) <<' '  // axis ratio c/a
	  << print(sqrt(ID[1]/ID[0]),9,3) <<' '; // axis ratio b/a
      Transpose<3>(IV);
      PS.print_dir(OUT,IV[0]) << ' ';
      PS.print_dir(OUT,IV[2]) << ' ';
      PS.print_dir(OUT,erot ) << std::endl;
    }
    return false;
  }
  //////////////////////////////////////////////////////////////////////////////
} }
////////////////////////////////////////////////////////////////////////////////
__DEF__MAN(falcON::Manipulate::densprof)

