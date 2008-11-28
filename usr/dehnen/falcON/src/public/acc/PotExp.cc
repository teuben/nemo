//-----------------------------------------------------------------------------+
//                                                                             |
/// \file src/public/acc/PotExp.cc                                             |
//                                                                             |
// Copyright (C) 2004-2008 Walter Dehnen                                       |
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
#include <defacc.h>
#include <ctime>
#define falcON_NEMO
#include <public/PotExp.h>
#include <public/nemo++.h>
////////////////////////////////////////////////////////////////////////////////
namespace {
  using namespace falcON;
  typedef falcON::tupel<3,float>  vectf;
  typedef falcON::tupel<3,double> vectd;
  //----------------------------------------------------------------------------
  const int    Nexp = 10;                          // max # of pot allowed
  const double a_def=1, r_def=1;                   // default parameter
  const int    n_def=8, l_def=8, s_def=1;          //   values
  //----------------------------------------------------------------------------
  struct PwithC : public PotExp {
    PotExp::Anlm Coef;
    PwithC(scalar a, scalar r, int n, int l, PotExp::symmetry s) :
      PotExp(a,r,n,l,s), Coef(*this) {}
  };
  //----------------------------------------------------------------------------
  inline PotExp::symmetry sym(int s) {
    return
      s==4? PotExp::spherical :
      s==3? PotExp::cylindrical : 
      s==2? PotExp::triaxial :
      s==1? PotExp::reflexion : 
            PotExp::none;
  }
  //----------------------------------------------------------------------------
  class PotExpansion {
    PwithC *P;
  public:
    static const char* name() {
      return "PotExp";
    }
    //--------------------------------------------------------------------------
    PotExpansion() : P(0)
    {}
    //--------------------------------------------------------------------------
    void init(const double*, int, const char*);
    //--------------------------------------------------------------------------
    bool is_init() const {
      return P!=0;
    }
    //--------------------------------------------------------------------------
    ~PotExpansion() {
      delete P;
    }
    //--------------------------------------------------------------------------
    void acc(int        n,
	     const void*x,
	     const int *f,
	     void      *p,
	     void      *a,
	     int        d,
	     char       t) const falcON_THROWING
    {
      clock_t cpu0 = clock();
      switch(t) {
      case 'f':
	P->SetGravity(P->Coef, n,
		      static_cast<const vectf*>(x),
		      static_cast<float*>(p),
		      static_cast<vectf*>(a),
		      f,d);
	if(P->has_error  ())
	  falcON_THROWN (const_cast<char*>(P->error_msg()));
	if(P->has_warning())
	  falcON_Warning(const_cast<char*>(P->warning_msg()));
	break;
      case 'd':
	P->SetGravity(P->Coef, n,
		      static_cast<const vectd*>(x),
		      static_cast<double*>(p),
		      static_cast<vectd*>(a),
		      f,d);
	if(P->has_error  ())
	  falcON_THROWN (const_cast<char*>(P->error_msg()));
	if(P->has_warning())
	  falcON_Warning(const_cast<char*>(P->warning_msg()));
	break;
      default:
	falcON_THROWN ("%s unknown type : '%c'",t);
      }
      clock_t cpu1 = clock();
      DebugInfo(2,"PotExp: gravity computed in %f sec CPU time\n",
		(cpu1 - cpu0)/double(CLOCKS_PER_SEC));
    }
  } // class PotExpansion
  Pexp[Nexp];                                      // array of Nexp PotExpansion
  //----------------------------------------------------------------------------
  void PotExpansion::init(const double*pars,
			  int          npar,
			  const char  *file)
  {
    // 0 checking consistency of arguments
    if(npar < 7)
      falcON_Warning(
      "%s: recognizing 7 parameters and one data file.\n"
      "Parameters:\n"
      " omega (real)           pattern speed (ignored)              [0]\n"
      " alpha (real)           shape parameter of expansion basis   [%f]\n"
      " r0    (real)           scale radius of expansion basis      [%f]\n"
      " nmax  (integer > 0)    max n in radial expansion            [%d]\n"
      " lmax  (integer, even)  max l in angular expansion           [%d]\n"
      " symm  (integer)        symmetry assumed (see below)         [%d]\n"
      " G     (real)           constant of gravity                  [1]\n\n"
      "The potential is given by the expansion\n\n"
      "    Phi(x) =  Sum  C_nlm Phi     (x)\n"
      "             n,l,m          n,l,m\n\n"
      "with the basis functions\n\n"
      "    Phi_nlm = - Psi_nl(r) * Y_lm(theta,phi).\n\n"
      "The lowest order radial basis function is given by\n\n"
      "                      (1/a)     -a\n"
      "    Psi_000 = ( [r/r0]      + 1)\n\n"
      "which gives a Hernquist sphere for a=alpha=1 and a Plummer sphere for a=1/2.\n"
      "The coefficients are such that potential approximates that of the first\n"
      "snapshot found in the data file.\n"
      "The last parameter, symm, allows to symmetrize the potential by constraining\n"
      "the coefficients:\n"
      " symm=0:   no symmetry: all coefficients used\n"
      " symm=1:   reflexion wrt origin: C_nlm=0 for odd (l,m)\n"
      " symm=2:   triaxial wrt xyz axes: C_nlm=0 for odd (l,m) and C_nlm = C_nl[-m]\n"
      " symm=3:   cylindrical: C_nlm=0 for odd l or m!=0\n"
      " symm=4:   spherical: C_nlm=0 for (l,m) != 0\n",
      name(),a_def,r_def,n_def,l_def,s_def);
    if(file==0 || file[0]==0)
      falcON_THROWN("%s: data file required\n","PotExp");
    // 1 reading in parameters and initializing potential expansion
    double
//       o = npar>0? pars[0] : 0,
      a = npar>1? pars[1] : a_def,
      r = npar>2? pars[2] : r_def;
    int
      n = npar>3? int(pars[3]) : n_def,
      l = npar>4? int(pars[4]) : l_def,
      s = npar>5? int(pars[5]) : s_def;
    double
      G = npar>6? pars[6] : 1.;
    if(s < 0 || s > 4) {
      falcON_WarningN("%s: symm out of range, defaulting to %d (%s symmetry)\n",
		      name(),s_def,PotExp::name_of_sym(sym(s_def)));
      s = s_def;
    }
    if(npar>7) falcON_WarningN("%s: skipped parameters beyond 6",name());
    P = new PwithC(a,r,n,l,sym(s));
    if(P->has_error  ()) falcON_THROWN (const_cast<char*>(P->error_msg()));
    if(P->has_warning()) falcON_Warning(const_cast<char*>(P->warning_msg()));
    DebugInfo(2,
	      "PotExp: initialized expansion with\n"
	      " alpha = %f\n"
	      " r0    = %f\n"
	      " nmax  = %d\n"
	      " lmax  = %d\n"
	      " assuming %s symmetry\n",
	      P->alpha(), P->scale(), P->Nmax(), P->Lmax(), 
	      P->symmetry_name());
    // 2 reading in positions and masses from snapshot
    nemo_in innemo(file);
    if(!innemo.is_open())
      falcON_ErrorN("PotExp: cannot open file %s\n",file);
    if(!innemo.has_snapshot())
      falcON_ErrorN("PotExp: no snapshot in file %s\n",file);
    snap_in insnap(innemo);
    int N = insnap.Ntot();
    // 2.1 read positions:
    vectf *x = falcON_NEW(vectf,N);
    if(insnap.has(nemo_io::pos)) {
      data_in indata(insnap,nemo_io::pos);
      if(indata.type() == nemo_io::Single) 
	indata.read(x);
      else if(indata.type() == nemo_io::Double) {
	vectd*_x = falcON_NEW(vectd,N);
	indata.read(_x);
        for(int i=0; i!=N; ++i) x[i] = _x[i];
        falcON_DEL_A(_x);
      } else
	falcON_ErrorN("PotExp: position not in float nor double\n");
    } else if(insnap.has(nemo_io::posvel)) {
      data_in indata(insnap,nemo_io::posvel);
      if(indata.type() == nemo_io::Single) {
	vectf*_x = falcON_NEW(vectf,2*N);
	indata.read(_x);
        for(int i=0,j=0; i!=N; ++i,j+=2) x[i] = _x[j];
        falcON_DEL_A(_x);
      } else if(indata.type() == nemo_io::Double) {
	vectd*_x = falcON_NEW(vectd,2*N);
	indata.read(_x);
        for(int i=0,j=0; i!=N; ++i,j+=2) x[i] = _x[j];
        falcON_DEL_A(_x);
      } else
	falcON_ErrorN("PotExp: phases not in float nor double\n");
    } else
      falcON_ErrorN("PotExp: no positions found in snapshot\n");
    // 2.2 read masses:
    float *m = falcON_NEW(float,N);
    if(insnap.has(nemo_io::mass)) {
      data_in indata(insnap,nemo_io::mass);
      if(indata.type() == nemo_io::Single) 
	indata.read(m);
      else if(indata.type() == nemo_io::Double) {
	double*_m = falcON_NEW(double,N);
	indata.read(_m);
        for(int i=0; i!=N; ++i) m[i] = _m[i];
        falcON_DEL_A(_m);
      } else
	falcON_ErrorN("PotExp: masses not in float nor double\n");
    } else
      falcON_ErrorN("PotExp: no masses found in snapshot\n");
    // 3 initializing coefficients
    clock_t cpu0 = clock();
    P->Coef.reset();
    P->AddCoeffs(P->Coef,N,m,x,0);
    if(P->has_error  ()) falcON_THROWN (const_cast<char*>(P->error_msg()));
    if(P->has_warning()) falcON_Warning(const_cast<char*>(P->warning_msg()));
    P->Normalize(P->Coef,G);
    if(P->has_error  ()) falcON_THROWN (const_cast<char*>(P->error_msg()));
    if(P->has_warning()) falcON_Warning(const_cast<char*>(P->warning_msg()));
    clock_t cpu1 = clock();
    DebugInfo(2,"PotExp: coefficients computed in %f sec CPU time\n",
	      (cpu1 - cpu0)/double(CLOCKS_PER_SEC));
    if(nemo_debug(2)) {
      std::cerr<<"PotExp: coefficients:\n";
      P->Coef.table_print(P->Symmetry(),std::cerr);
    }
    // 4 clearing up
    delete[] x;
    delete[] m;
  } // PotExpansion::init()
  //----------------------------------------------------------------------------
#define ACCFUNC(NUM)							\
  void accel##NUM(int        ndim,					\
		  double     time,					\
		  int        n,						\
		  const void*m,						\
		  const void*x,						\
		  const void*v,						\
		  const int *f,						\
		  void      *p,						\
		  void      *a,						\
		  int        d,						\
		  char       type)					\
  {									\
    if(ndim != 3)							\
      falcON_Error("%s: ndim=%d not supported\n",			\
		   PotExpansion::name(),ndim);				\
    Pexp[NUM].acc(n,x,f,p,a,d,type);					\
  }
  ACCFUNC(0);                           // accel0() to accel(9)
  ACCFUNC(1);                           // use Pexp[i].acc()
  ACCFUNC(2);  
  ACCFUNC(3);  
  ACCFUNC(4);  
  ACCFUNC(5);  
  ACCFUNC(6);  
  ACCFUNC(7);  
  ACCFUNC(8);  
  ACCFUNC(9);  

  // array of Nexp acc_pter: accel0 to accel9
  acc_pter ACCS[Nexp] = {accel0, accel1, accel2, accel3, accel4,
			 accel5, accel6, accel7, accel8, accel9};
  int Iexp = 0;
} // namespace MyPotExp {

void iniacceleration(
		     const double*pars,
		     int          npar,
		     const char  *file,
		     acc_pter    *accf,
		     bool        *needm,
		     bool        *needv)
{
  if(Iexp == Nexp)
    falcON_Error("iniacceleration(): "
		 "cannot have more than %d instances of '%s'\n",
		 Nexp,PotExpansion::name());
  if(needm) *needm = 0;
  if(needv) *needv = 0;
  Pexp[Iexp].init(pars,npar,file);
  *accf = ACCS[Iexp++];
}
//------------------------------------------------------------------------------
